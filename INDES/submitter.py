''' this module submits the jobs and monitors them '''

import os
import random
import time
import string
import tempfile
import subprocess


def submit(job, script='ID_gauss'):
    print('filename:'),
    print(job.filename)
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
        # cmd='g09'
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
    except subprocess.CalledProcessError, e:
        print "subprocess.CalledProcessError"
        print repr(e)
        p1 = False
    return p1


if __name__ == "__main__":
    import sys
    file = sys.argv[1]
    jobid = str(submit(index))
    print "jobid: ", jobid
    print "output of 'qstat | grep 020':\n", jobstatus(jobid)
