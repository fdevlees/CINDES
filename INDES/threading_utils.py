
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

def runjobs_threading(mols_tocal, myrun, i):
    global qsta_out

    calc = myrun.calcs[i]

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
    call(jobmaker, mols=mols_calc, myrun=myrun, calc=calc)

    jobs = get_jobs(myrun, mols_calc)

    # From here can be threaded:

    # 1. make a list of Threads. one for each Job
    T = [ MyThread(job=job) for job in jobs ]
    print T

    # 2. run each thread
    for t in T: t.start()
    print T

    # 3. check the state of each job:
    while True:
        n=0
        for t in T:
            if t.isAlive():n+=1
        logging.debug("{:d} threads are alive".format(n))

        qsta_out = get_qsta_out()

        if n==0:
            print "all threads are ready!"
            break
        time.sleep(2)

    print T
    print "jobs:", jobs

    # 3. test of all jobs are ready
    if jobids:
        subm.jobtester(mols_calc, myrun, jobids)
    else:
        print "no jobids so assume no jobs submitted"

    # 4. test normal termination and read jobs (only myrun variable used is actually debug)
    datareader.datareader(mols_calc, myrun)

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
             job=None, verbose=None):
        threading.Thread.__init__(self, group=group, target=target, name=name,
                              verbose=verbose)
        self.job = job
        return

    def run(self):
        logging.debug('Starting')

        # run job:
        self.job.submit()

        # check if ready()
        while True:
            if self.test_ready():
                print "job = ready"
                break
            else:
                time.sleep(5)

        # check if normal termination.
        

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
