''' this module submits the jobs and monitors them '''

import os
import random
import time
import string
import tempfile
import subprocess
import logging
import glob

from CINDES.utils.writings import log_io
once = 0
debug = False
safe = False
counter = 0

def submit(job, script='ID_gauss'):
    #print('filename:'),
    #print(job.filename)
    command = './' + script
    try:
        jobid = subprocess.check_output([command, job.filename], cwd=job.path)
        if not jobid:
            raise RuntimeError('no jobid')
    except subprocess.CalledProcessError as e:
        print "submitting error:", repr(e)
    except Exception as e:
        print "unforeseen submission error:", repr(e)
        raise
    time.sleep(1)
    return jobid

def submitworker(nprocs):
    # command to submit on 2 procs: wsub --batch myworker9.pbs --data loglist.csv -threaded 2 -master
    if nprocs==1:
        jobids = subprocess.check_output(['wsub', '-batch', 'CINDES_worker.pbs', '-data', 'loglist.csv'])
    elif nprocs==2:
        jobids = subprocess.check_output(['wsub', '-batch', 'CINDES_worker.pbs', '-data', 'loglist.csv', '-threaded', '2', '-master'])
    elif nprocs==4:
        jobids = subprocess.check_output(['wsub', '-batch', 'CINDES_worker.pbs', '-data', 'loglist.csv', '-threaded', '4', '-master'])
    else:
        raise NotImplementedError('other than 1,2,4 procs is currently not allowed via worker submission')
    return jobids

def nosubmit_orca(path, index, identify):
    inputname = identify + index
    outname = path + '/' + identify + index + '.out'
    fakename = path + '/' + identify + index + '.o123456'
    with open(outname, 'w') as out, open(fakename, 'w') as err:
        p = subprocess.Popen(['orca', inputname], stdout=out, stderr=err, cwd=path)
        p.wait()
    return 123456


def nosubmit(job, extension='.com', cmd=None):  # not tested
    if not cmd:
        #cmd='g09'
        cmd = 'g16'
    path = job.path
    inputname = job.filepath
    outname = job.logpath
    fakename = job.filepath + '.o12345'  # NOTE: on hydra this is without the [:-4] (.com)
    print "fakename:", fakename
    with open(inputname, 'r') as inp, open(outname, 'w') as out, open(fakename, 'w') as err:
        try:
            p = subprocess.Popen(cmd, stdin=inp, stdout=out, stderr=err, cwd=path)
            p.wait()
        except OSError as e:
            print "module probably not loaded."
            raise
    return 123456


def jobstatus(jobid):
    """ this function is not used anymore because it occured once that jobids changed during the calculation
    or that jobs need to be resubmitted because of torque problems. In that case jobids change
    Jobs are now declared ready when there exists a <name>.com.<6digits> file. (which can be empty)"""
    p1 = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep', '-w', '%s' % jobid], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0]
    return output


def qsta():
    try:
        p1 = subprocess.check_output(['qsta'])
    except subprocess.CalledProcessError as e:
        print "subprocess.CalledProcessError"
        print repr(e)
        p1 = False
    except LookupError as e:
        # This error occurred on VSC in subprocess module calling pickle using 'string_escape'
        print "LookupError"
        print repr(e)
        p1 = False
    return p1


# 2. submission
@log_io()
def submission(mols_tocal, myrun):
    global once
    if once == 1 and myrun.no1sub == 1:
        print "submit skipped"
        once = 2
        jobids = None
    else:
        jobids = submit_normal(mols_tocal, myrun)  # In here is decided to run on shell or to really submit!
    logging.info("----- END all jobs are submitted ----------")
    if safe:
        time.sleep(15)  # wait 15 seconds. to be sure that the jobs appear in the qstat command
    return jobids


def submit_normal(mols_tocal, myrun):
    ''' this function submits all the jobs of every mol.jobs list
    when worker=True the worker framework will be used:
        first a list of jobs is gathered in jobids and this list is used
        to submit the jobs simultaneously

        worker/1.6.8-intel-2018a
        '''

    jobids = []
    worker = myrun.worker
    for molecule in mols_tocal:
        for job in molecule.jobs:

            # 1. try ready part
            if myrun.try_ready:
                if worker:
                    name = job.logpath
                    if glob.glob(name):
                        print "already calculated:", name
                        continue
                else:
                    name1 = job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                    name2 = job.filepath + '.o[0-9][0-9][0-9][0-9]*'
                    if glob.glob(name1) or glob.glob(name2):
                        print "already calculated:", name2
                        continue

            # 2. submit part
            if worker:
                jobids.append(job)
            else:
                jobid = job.submit()
                jobids.append(jobid)
                time.sleep(1)

    minimum_number_worker_jobs = 10
    if worker and len(jobids)>minimum_number_worker_jobs:
        nprocs  = jobids[0].calc['nprocs'] # should be the same for all jobs!
        # 1. write csv file with job paths
        with open('loglist.csv', 'w') as f:
            f.write('job,log\n')
            for job in jobids:
                f.write('{},{}\n'.format(job.filepath,job.logpath))

        # 2. write pbs file with custom number of nodes:
        # always use at least one node and otherwise use njobs/28 or njobs/14 if nprocs==2
        nnodes= max(1, len(jobids)/(28/nprocs))
        # to prevent using similar names in queue
        global counter
        identifier = "{}{}".format(myrun.identify,str(counter))
        counter += 1
        with open('CINDES_worker.template','r') as f:
            template = f.read()
        with open('CINDES_worker.pbs','w') as f:
            f.write(template.format(nnodes=nnodes, identifier=identifier))

        submitworker(nprocs)
        jobfiles = [ identifier ]
        return jobfiles
    elif worker and len(jobids)>0:
        print "there are too few jobs to use the worker efficiently so jobs will be submitted to full nodes"
        jobfiles = []
        for job in jobids:
            #rewrite job to have %nprocshared=28
            job.rewrite()
            job.submit()
            jobfiles.append(job.filename)
            time.sleep(1)
        return jobfiles
    else:
        return jobids

# 3. testing


@log_io(signator='=')
def jobtester(mols_tocal, myrun, jobids=None):
    """ this tester tests if the jobs are ready by looking for a file <name><.extension>.o<6digits>.

        - even if try_ready is activated all indices are used. And the already ready ones are immediately recognized as ready.
        - they are just not submitted again.
        - note that jobids are not used!
        - function returns nothing but returns when all jobs are ready! this function therefore can take very long!
    """
    if jobids is None:
        jobids = []  # this to avoid the mutable default gotcha
    if myrun.nosub == 2 or len(mols_tocal) == 0:
        print "Job tester skipped because jobs are evaluated on login node or no jobs to be calculated"
        return
    print "len(mols_tocal):", len(mols_tocal)
    test_ready = myrun.test_ready
    path = myrun.path
    fileparameters = myrun.__dict__
    if test_ready == 1:
        test_ready1(mols_tocal, myrun)
    elif test_ready == 2:  # default
        test_ready2(mols_tocal, myrun, jobids)
    elif test_ready == 3:
        test_ready3(mols_tocal, myrun)
    else:
        raise SystemExit('no valid test_ready value')
    return


def test_ready1(indices, myrun):
    ''' test ready based on the presence of inputfile.o$$ file '''
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = []  # here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + \
            myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab'] == 1:  # property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + \
                    '_' + str(pos) + myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
                paths.append(path2)
    while True:  # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths == []:
            break
        print "time/h:", tijdje / 3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return


def test_ready2(mols_tocal, myrun, jobids=None):
    """
    this function tests if the jobs are still queing or running based on the output of the 'qsta' command
    function needs:
    mols:
    -jobs
    myrun:
    -timelimit
    -timestep
    -extrawaittime
    """
    completedjobs = []
    files = []
    fileparameters = myrun.__dict__

    # 1. get the list of entries that are in the queue
    if myrun.worker:
        #files = ['my-gaussian-worker-job']
        #files = ['CINDES_ATOOLS']
        files = jobids
    else:
        for mol in mols_tocal:
            jobnames = [job.filename for job in mol.jobs]
            files.extend(jobnames)
    print "files:", files

    # 2. and wait until all are completed or not anymore in queue
    tijdje = 0
    while True:
        count = 0
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        filescopy = files[:]
        njobs = len(filescopy)
        qsta_raw = qsta()
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
        if debug:
            print "states:", states
            print "jobs:", jobs
            print "files:", files
        for filetje in filescopy:
            for state, job in zip(states, jobs):
                if filetje == job:
                    if state in ['Q', 'R']:  # so if job still in queue and not has state=='C'
                        count += 1  # so count all the jobs still in queue
                    elif state in ['H', 'E']:
                        print "ERROR jobs on hold or Error"
                        count += 1
                    else:
                        assert state == 'C'
                        if job not in completedjobs:
                            print 'JOB completed:', job
                            completedjobs.append(job)
                    if debug:
                        print "found a job: ", state, job, filetje
        t = "{:.2f}".format(tijdje / 3600.)
        print "njobs -running: {:d} -ready: {:d} | waittime={} hrs".format(count, njobs - count, t)
        if count == 0:  # so no jobs anymore in queue
            break
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return


def test_ready3(indices, myrun):
    '''this could be something using a line as :
        touch ${PBS_JOBID}.completed
    in the jobscript ID_gauss
    '''
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = []  # here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.completed'
        paths.append(path1)
        if fileparameters['stab'] == 1:  # property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + \
                    indices[i] + '_' + str(pos) + '.completed'
                paths.append(path2)
    while True:  # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths == []:
            break
        print "time/h:", tijdje / 3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return

if __name__ == "__main__":
    import sys
    file = sys.argv[1]
    jobid = str(submit(index))
    print "jobid: ", jobid
    print "output of 'qstat | grep 020':\n", jobstatus(jobid)
