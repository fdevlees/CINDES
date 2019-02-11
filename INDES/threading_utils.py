
import threading
from calculator import jobmaker, invoke_script
import submitter as subm
import datareader
import logging
import time
logging.basicConfig(level=logging.DEBUG,
                    format='(%(threadName)-10s) %(message)s',
                    )

qsta_out = dict()

debug = True

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
    mols_calc = filter(lambda x: not x.ignoremol, mols_tocal)

    # 1. Make the jobs and add them to the molecules:
    call(jobmaker, mols=mols_calc, run=run, calc=calc)

    jobs = get_jobs(run, mols_calc)

    # From here can be threaded:

    # 1. make a list of Threads. one for each Job
    T = [ MyThread(runparam=run, job=job) for job in jobs ]
    print T

    # 2. run each thread
    for t in T: t.start()
    print T

    # 3. check the state of each job:
    while True:
        n=0
        print "current threads:", threading.enumerate()
        for t in T:
            if t.isAlive():n+=1
        logging.debug("{:d} threads are alive".format(n))

        qsta_out = get_qsta_out()

        if n==0:
            print "all threads are ready!"
            break
        time.sleep(2)

    # now combine all information from the jobs:
    for molecule in mols_calc:
        molecule.predicted = False
        for job in molecule.jobs:
            molecule.props.update(job.readings)

    print T
    print "jobs:", jobs

    # 4. test normal termination and read jobs (only run variable used is actually debug)
    datareader.datareader(mols_calc, run)

    # until here

    return


def get_jobs(run, mols):
    jobs = []
    for molecule in mols:
        for job in molecule.jobs:

            # 1. try ready part
            if run.try_ready:
                name1 = job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                name2 = job.filepath + '.o[0-9][0-9][0-9][0-9]*'
                if glob.glob(name1) or glob.glob(name2):
                    print "already calculated:", name2
                    continue
            jobs.append(job)
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

    def run(self):
        logging.debug('Starting')

        # run job:
        self.job.submit()

        # check if ready()
        while True:
            time.sleep(5)
            if self.test_ready():
                print "job = ready"
                break
            else:
                if debug: print 'not ready {}'.format(self.job.filename)

        # check if normal termination. If not a zzz.log is tried to be made and submitted
        self.job.normaltermination(debug=self.runparam.debug)
        if self.job.IsReady:
            print "normal termination for:", self.job.filename
        else:
            print "no normal termination?!", self.job.filename
            # test normal termination and errorjob are ready or molecule is ignored
            while True:
                self.job.ready(ignore=self.runparam.ignore)
                if job.IsReady:
                    break
                time.sleep(5)

        # read data
        readings = datareader.read_file(self.job)

        # 2.2 if there are multiple variants of the job, give each variant a different index _P#
        if hasattr(self.job, 'pos'):
            for old_key in readings.keys():  # the .keys is very important here. iterkeys or for just readings do not work!
                readings["{}_P{}".format(old_key, str(job.pos))] = readings.pop(old_key)

        self.job.readings = readings

        logging.debug('Exiting')
        return True

    def test_ready(self):
        try:
            if not qsta_out[self.job.filename] in ['H', 'E', 'R', 'Q']:
                return True
        except KeyError:
            print "KeyError:", self
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
            print "qsta not working!"
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


    return mols_toread


def zzztester(mols):
    """
    this function tests if the jobs are still queing or running based on the output of the 'qsta' command
    function needs:
    mols:
    -jobs
    -run:
    -timestep
    """
    files = []
    for mol in mols:
        jobnames = [job.errorfile for job in mol.jobs if not job.errorfile is None]
        files.extend(jobnames)
    if files:
        #print "zzz-files:", files
        print "n zzz files:", len(files),
    else:  # here return so we don't need the qsta
        return False

    # get qstat
    while True:
        qsta_raw = qsta()
        if qsta_raw == False:
            print "qsta not working! trying again after 1 minute"
            time.sleep(60)
        else:
            break

    # get list of jobs and list of their states
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

    # inefficient loop
    #print "files:", files
    #print "jobs:", jobs
    #print "states:", states
    count = 0
    for filetje in files:
        for state, job in zip(states, jobs):
            if job in filetje:  # so there is a zzzjob in the queue!
                # filetje is whole path so job in filetje or filetje.split('/')[-1]==job
                if state in ['Q', 'R', 'H', 'E']:  # so if job still in queue and not has state=='C'
                    count += 1  # so count all the jobs still in queue
                else:
                    assert state == 'C'

    # return answer
    if count > 0:
        return True
    else:
        print "no zzzs (anymore) in queue",
        return False
