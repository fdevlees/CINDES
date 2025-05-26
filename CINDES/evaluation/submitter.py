''' this module submits the jobs and monitors them '''

import os
import random
import time
import string
import tempfile
import subprocess
import logging
import glob
import json

from CINDES.utils.writings import log_io
once = 0
debug = False
safe = False
counter = 0

def submit(job, script='ID_gauss'):
    #print('filename:'),
    #print(job.filename)
    command = './' + script
    try:#@@@david: changes so you can have variable walltime per job
        #jobid = subprocess.check_output([command, job.filename], cwd=job.path)
        try:
            walltime=(vars(job))["calc"]["walltimelimit"]
            if walltime[-1]=='m' or walltime[-1]=='M':
                walltime='00:'+walltime[:-1]
            elif walltime[-1]=='h' or walltime[-1]=='H':
                walltime=walltime[:-1]+':00'
        except:
            print("No walltime keyword found for job: default time of 24h assumed")
            walltime='24:00'
        print('processed walltime: ',walltime)
        print('ddebug: jobid=|command:'+str(command)+"|job.filename:"+str(job.filename)+"|walltime="+walltime+"|cwd=:"+job.path)#@@@david
        jobid = subprocess.check_output([command, job.filename,walltime], cwd=job.path)#@@

        if not jobid:
            raise RuntimeError('no jobid')
    except subprocess.CalledProcessError as e:
        print("submitting error:", repr(e))
    except OSError as e:
        print("submission error; submission script might not be executable?")
        scriptpath = "{}/{}".format(job.path, script)
        st = os.stat(scriptpath)
        print("permissions of {} is:".format(scriptpath), st)
        import stat
        os.chmod(scriptpath, st.st_mode | stat.S_IEXEC)
    except Exception as e:
        print("unforeseen submission error:", repr(e))
        raise
    time.sleep(1)
    return jobid

def submitworker(nprocs):
    # command to submit on 2 procs: wsub --batch myworker9.pbs --data loglist.csv -threaded 2 -master
    if nprocs==1:
        jobids = subprocess.check_output(['wsub', '-batch', 'CINDES_worker.pbs', '-data', 'loglist.csv'])
    elif nprocs>1:
        jobids = subprocess.check_output(['wsub', '-batch', 'CINDES_worker.pbs', '-data', 'loglist.csv', '-threaded', str(nprocs), '-master'])
    else:
        print("nprocs:", nprocs)
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
    print("fakename:", fakename)
    with open(inputname, 'r') as inp, open(outname, 'w') as out, open(fakename, 'w') as err:
        try:
            p = subprocess.Popen(cmd, stdin=inp, stdout=out, stderr=err, cwd=path)
            p.wait()
        except OSError as e:
            print("module probably not loaded.")
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
        output = subprocess.check_output(['qsta']).decode('utf-8')
        jobs = json.loads(output)
#        print("Jobs ontvangen:", jobs)
        return jobs
    except subprocess.CalledProcessError as e:
        print("Error with executing qsta:", e)
    except json.JSONDecodeError as e:
        print("Error with parsing JSON:", e)
    except Exception as e:
        print("Unexpected error:", e)
        return []


# 2. submission
@log_io()
def submission(mols_tocal, myrun):
    global once
    if once == 1 and myrun.no1sub == 1:
        print("submit skipped")
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
                        print("already calculated:", name)
                        continue
                else:
                    name1 = job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                    name2 = job.filepath + '.o[0-9][0-9][0-9][0-9]*'
                    nameG16log= job.filepath[:-4] + '.log' #@@@david: added so it also works with .log Gaussian output files (so no o files needed)
                    if glob.glob(name1) or glob.glob(name2) or glob.glob(nameG16log):
                        print("already calculated:", job.filepath[:-4])#@@
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
        print("there are too few jobs to use the worker efficiently so jobs will be submitted to full nodes")
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
    if myrun.nosub or len(mols_tocal) == 0:
        print("Job tester skipped because jobs are evaluated on login node or no jobs to be calculated")
        return
    print("len(mols_tocal):", len(mols_tocal))
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
    """
    Controleert of jobs klaar zijn op basis van het bestaan van outputbestanden (*.o******).
    """
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = []

    # ✅ Veilige standaardwaarden
    try:
        timestep = float(fileparameters.get('timestep', 30))
        if timestep <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for timestep - default value of 30s is being used.")
        timestep = 300

    try:
        timelimit = float(fileparameters.get('timelimit', 86400))  # 24 uur
        if timelimit <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for timestep — default value of 86400s (=24h) is being used.")
        timelimit = 86400

    try:
        extrawait = float(fileparameters.get('extrawaittime', 10))
        if extrawait < 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for extrawaittime - default value of 10s is being used.")
        extrawait = 10

    # ✅ Bouw padlijst op
    for i in range(len(indices)):
        base = f"{path}/{fileparameters['identify']}{indices[i]}{myrun.extension}.o[0-9][0-9][0-9][0-9][0-9][0-9]"
        paths.append(base)

        if fileparameters.get('stab', 0) == 1:
            for pos in fileparameters.get('positions', []):
                sub = f"{path}/{indices[i]}/{fileparameters['identify']}{indices[i]}_{pos}{myrun.extension}.o[0-9][0-9][0-9][0-9][0-9][0-9]"
                paths.append(sub)

    # ✅ Wacht tot alle outputbestanden bestaan
    while True:
        if tijdje > timelimit:
            print("⏰ Time limit has been reached.")
            break

        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print("✅ Ready:", pathje[:-25])

        if not paths:
            break

        print(f"⏳ Time (hour): {tijdje / 3600:.2f} | Still to be checked: {len(paths)}")
        time.sleep(timestep)
        tijdje += timestep

    logging.info("✅ All jobs are ready.")
    time.sleep(extrawait)
    return


def test_ready2(mols_tocal, run, jobids=None, debug=False):
    """
    Controleert of jobs nog in de wachtrij staan of aan het draaien zijn op basis van de JSON-output van 'qsta'.
    """
    completedjobs = []
    files = []


    # ✅ Veilige standaardwaarden
    try:
        timestep = float(getattr(run, 'timestep', 300))
        if timestep <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for timestep - default value of 30s is being used.")
        timestep = 300

    try:
        timelimit = float(getattr(run, 'timelimit', 86400))  # 24 uur
        if timelimit <= 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for timestep — default value of 86400s (=24h) is being used.")
        timelimit = 86400

    try:
        extrawait = float(getattr(run, 'extrawaittime', 10))
        if extrawait < 0:
            raise ValueError
    except (ValueError, TypeError):
        print("⚠️ Invalid value for extrawaittime - default value of 10s is being used.")
        extrawait = 10


    # 1. Bepaal welke jobnamen we moeten controleren
    if run.worker:
        files = jobids
    else:
        for mol in mols_tocal:
            jobnames = [job.filename for job in mol.jobs]
            files.extend(jobnames)

    print("Files to be checked:", files)

    # 2. Wacht tot alle jobs klaar zijn
    tijdje = 0

    while True:
        count = 0
        if tijdje > timelimit:
            print("Time limit had been reached.")
            break

        filescopy = files[:]
        njobs = len(filescopy)
        qsta_data = qsta()

        if not qsta_data:
            print("qsta does not work!")
            time.sleep(timestep)
            continue

        # 3. Verwerk de JSON-output
        states = []
        jobs = []
        for job in qsta_data:
            state = job.get("State")
            name = job.get("Name")
            if state and name:
                states.append(state)
                jobs.append(name)

 #       print("Statussen:", states)
 #       print("Jobnamen:", jobs)

        if debug:
            print("Status:", states)
            print("Job name:", jobs)
            print("Files to be checked:", files)

        for filetje in filescopy:
            for state, job in zip(states, jobs):
                if filetje == job:
                    if state in ['PD', 'R', 'CG']:  # Pending or Running or Completing
                        count += 1
                    elif state in ['F', 'TO', 'NF', 'SE']:  # Failed or Timeout or Node Fail or Special Exit
                        print("ERROR: Job failed or timeout:", job)
                        count += 1
                    else:
                        if job not in completedjobs:
                            print('JOB completed:', job)
                            completedjobs.append(job)
                    if debug:
                        print("Job found:", state, job)

        t = "{:.2f}".format(tijdje / 3600.)
#        print(f"⏳ Jobs totaal: {njobs} | Nog bezig: {count} | Klaar: {njobs - count} | Wachttijd: {t} uur")

        if count == 0:
#            print("🟡 Alle jobs lijken klaar  breek uit de lus.")
            break

        time.sleep(timestep)
        tijdje += timestep

#    print("✅ Alle jobs zijn klaar.")
    logging.info("✅ All jobs are ready.")
    time.sleep(extrawait)
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
            print("time is up")
            break
        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print("ready: ", pathje[:-25])
        if paths == []:
            break
        print("time/h:", tijdje / 3600, "len paths:", len(paths), end=' ')
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return

if __name__ == "__main__":
    qsta_out = qsta()
    print(json.dumps(qsta_out, indent=2))
