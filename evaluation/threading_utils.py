
import threading
from .calculator import jobmaker, invoke_script
from . import submitter as subm
from . import datareader
import logging
import time
logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-10s) %(message)s',
                    )

qsta_out = dict()

debug = False

def runjobs_threading(mols_tocal, run, i):
    global qsta_out

    calc = run.calcs[i]

    def call(function, calc, *args, **kwargs):
        if isinstance(calc, list) or isinstance(calc, tuple):
            for j, cal in enumerate(calc):
                invoke_script(cal, locals(), (i + 1) * 100 + (j + 1))
                function(calc=cal, *args, **kwargs)
        else:
            invoke_script(calc, locals(), i + 1)
            function(calc=calc, *args, **kwargs)

    # 0. filter off ignored molecules
    mols_calc = [x for x in mols_tocal if not x.ignoremol]

    # 1. Make the jobs and add them to the molecules:
    call(jobmaker, mols=mols_calc, run=run, calc=calc)


    # From here can be threaded:
    jobs = get_jobs(run, mols_calc)

    # 1. make a list of Threads. one for each Job
    T = [ MyThread(name=jobindex, runparam=run, job=job) for jobindex, job in jobs ]

    # 2. run each thread
    for t in T: t.start()
    print(T)

    # 3. check the state of each job:
    waittime = 0.0
    timestep = 2
    while True:
        n=0
        for t in T:
            if t.isAlive():n+=1

        # renew qsta_out. This is a global variable accessed by all threads!
        qsta_out = get_qsta_out()

        if n==0:
            print("all threads are ready!")
            break

        # wait
        waittime += timestep
        time.sleep(timestep)
        if waittime < 300:
            timestep *= 1.2
            logging.info("{:d} threads are alive. (waittime: {:8.2f}s)".format(n, waittime))
        else:
            logging.info("{:d} threads are alive. (waittime: {:8.2f}h)".format(n, waittime/3600.))

    # now combine all information from the jobs:
    for molecule in mols_calc:
        molecule.predicted = False
        if any(job.ignorejob for job in molecule.jobs):
            molecule.discard()
            continue
        for job in molecule.jobs:
            molecule.props.update(job.readings)
        # only relevant for stab calculations
        if True:
            datareader.setEAHs(molecule)

    print(T)
    print("jobs:", jobs)

    # 4. test normal termination and read jobs (only run variable used is actually debug)
    # datareader.datareader(mols_calc, run)

    # until here

    return


def get_jobs(run, mols):
    jobs = []
    for i, molecule in enumerate(mols):
        for j, job in enumerate(molecule.jobs):
            jobindex = "{:d}.{:d}".format(i,j)
            # 1. try ready part
            if run.try_ready:
                name1 = job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                name2 = job.filepath + '.o[0-9][0-9][0-9][0-9]*'
                nameG16log= job.filepath[:-4] + '.log' #@@@david: added so it works with .log Gaussian output files (so no o files needed)
                if glob.glob(name1) or glob.glob(name2) or glob.glob(nameG16log):
                    print("already calculated:", job.filepath[:-4])#@@
                    continue
            jobs.append((jobindex, job))
    return jobs

class MyThread(threading.Thread):
    def __init__(self, group=None, target=None, name=None,
             job=None,
             runparam=None,
             verbose=None):
        threading.Thread.__init__(self, group=group, target=target, name=name,
                              verbose=verbose)
        self.job = job
        self.runparam = runparam
        return

    def __repr__(self):
        return "<MyTread: {}>".format(self.name)

    def __str__(self):
       return self.__repr__()

    def run(self):
        logging.debug('Starting')

        # run job:
        self.job.submit()

        # check if ready()
        self.handle_termination()

        # read data
        self.extract_data()

        return True

    def extract_data(self):
        readings = datareader.read_file(self.job)

        # 2.2 if there are multiple variants of the job, give each variant a different index _P#
        if hasattr(self.job, 'pos'):
            for old_key in list(readings.keys()):  # the .keys is very important here. iterkeys or for just readings do not work!
                readings["{}_P{}".format(old_key, str(job.pos))] = readings.pop(old_key)

        self.job.readings = readings

        logging.debug('Exiting')
        return

    def handle_termination(self):
        while True:
            time.sleep(5)
            if self.test_ready():
                break
            else:
                if debug: print('not ready {}'.format(self.job.filename))

        # check if normal termination. If not a zzz.log is tried to be made and submitted
        self.job.normaltermination(debug=self.runparam.debug)
        if self.job.IsReady:
            print("normal termination for:", self.job.filename)
        else:
            print("no normal termination?!", self.job.filename)
            # test normal termination and errorjob are ready or molecule is ignored
            while True:
                if self.job.IsReady:
                    break
                self.job.ready(ignore=self.runparam.ignore)
                time.sleep(5)
        return



    def test_ready(self, error=False):
        try:
            if error:
                filename = self.job.errorfile
            else:
                filename = self.job.filename
            if not qsta_out[filename] in ['H', 'E', 'R', 'Q']:
                return True
        except KeyError:
            print("KeyError:", self)
            return True


def get_qsta_out():
    """
    this function tests if the jobs are still queing or running based on the output of the 'qsta' command
    """
    output = {}
    tijdje = 0
    while True:
        count = 0
        qsta_raw = subm.qsta()
        if qsta_raw == False:
            print("qsta not working!")
            time.sleep(fileparameters['timestep'])
            continue
        qsta_out = [item.split() for item in qsta_raw.split('\n')]
        states = []
        jobs = []
        for item in qsta_out:
            if len(item) == 4:
                states.append(item[1])
                jobs.append(item[3])
            elif len(item) == 5:
                states.append(item[2])
                jobs.append(item[4])
        for job, state in zip(jobs, states):
            output[job]=state
        break
    return output

