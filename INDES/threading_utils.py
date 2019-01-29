
import threading

def runjobs_threading(mols_tocal, myrun, i):
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

    # From here can be threaded:


    # 2. now the jobs have to be submitted (this function contains a try_ready test)
    jobids = subm.submission(mols_calc, myrun)

    # 3. test of all jobs are ready
    if jobids:
        subm.jobtester(mols_calc, myrun, jobids)
    else:
        print "no jobids so assume no jobs submitted"

    # 4. test normal termination and read jobs (only myrun variable used is actually debug)
    datareader.datareader(mols_calc, myrun)

    # until here

    return

def single_job_worker(run, job):

    pass
