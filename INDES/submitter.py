''' this module submits the jobs and monitors them '''

import os
import random
import time
import string
import tempfile
import subprocess

def submit(path,index,identify,script='ID_gauss'):
    filename = identify + index
    print( 'filename:'),
    print( filename)
    command = './' + script
    try:
        jobid = subprocess.check_output([command,filename],cwd=path)
    except subprocess.CalledProcessError as e:
        print "submitting error:", repr(e)
    time.sleep(1)
    return jobid

def nosubmit_orca(path,index,identify):
    inputname = identify + index
    outname =  path + '/' + identify + index + '.out'
    fakename = path + '/' + identify + index + '.o123456'
    with open(outname,'w') as out, open(fakename,'w') as err:
        p = subprocess.Popen(['orca',inputname],stdout=out,stderr=err,cwd=path)
        p.wait()
    return 123456

def nosubmit(path,index,identify, extension='.com'): # not tested
    #gaussiancmd='g09'
    gaussiancmd='g16'
    inputname = path + '/' + identify + index + extension
    outname =   path + '/' + identify + index + '.log'
    fakename =  path + '/' + identify + index + '.com.o123456'
    with open(inputname,'r') as inp, open(outname,'w') as out, open(fakename,'w') as err:
        p = subprocess.Popen(gaussiancmd, stdin=inp, stdout=out, stderr=err, cwd=path)
        p.wait()
    return 123456

def jobstatus(jobid):
    """ this function is not used anymore because it occured once that jobids changed during the calculation
    or that jobs need to be resubmitted because of torque problems. In that case jobids change
    Jobs are now declared ready when there exists a <name>.com.<6digits> file. (which can be empty)"""
    p1 = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep','-w','%s' % jobid],stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close() # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0]
    return output

def qsta():
    try:
        p1 = subprocess.check_output(['qsta'])
    except subprocess.CalledProcessError,e:
        print "subprocess.CalledProcessError"
        print repr(e)
        p1 = False
    return p1

if __name__ == "__main__":
    import sys
    file = sys.argv[1]
    jobid = str( submit(index))
    print "jobid: ", jobid
    print "output of 'qstat | grep 020':\n" , jobstatus(jobid)
