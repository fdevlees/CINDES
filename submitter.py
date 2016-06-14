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
    jobid = subprocess.check_output([command,filename],cwd=path)
    return jobid

def nosubmit_orca(path,index,identify):
    inputname = identify + index
    outname =  path + '/' + identify + index + '.out'
    fakename = path + '/' + identify + index + '.o123456'
    with open(outname,'w') as out, open(fakename,'w') as err:
        p = subprocess.Popen(['orca',inputname],stdout=out,stderr=err,cwd=path)
        p.wait()
    return 123456

def nosubmit(path,index,identify): # not tested
    #print "path", path
    inputname = path + '/' + identify + index + '.com'
    outname =   path + '/' + identify + index + '.log'
    fakename =  path + '/' + identify + index + '.com.o123456'
    with open(inputname,'r') as inp, open(outname,'w') as out, open(fakename,'w') as err:
        p = subprocess.Popen('g09',stdin=inp,stdout=out,stderr=err,cwd=path)
        p.wait()
    return 123456

def jobstatus(jobid):
    p1 = subprocess.Popen(['qstat'], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['grep','-w','%s' % jobid],stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close() # Allow p1 to receive a SIGPIPE if p2 exits.
    output = p2.communicate()[0]
    return output

if __name__ == "__main__":
    import sys
    file = sys.argv[1]
    jobid = str( submit(index))
    print "jobid: ", jobid
    print "output of 'qstat | grep 020':\n" , jobstatus(jobid)
