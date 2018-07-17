#!/bin/env python
'''a BaseClass for any kind of Job'''

import time
import submitter as subm

class BaseJob(object):
    '''BaseJob of any kind of submittable job'''
    script = ''

    def __init__(self, filepath, calc=dict()):
        '''There are different combinations of path and file to be set:
        -filepath = full path to INPUT file
        -path     = path where file is located
        -filename = name of inputfile with extension without path
        -name     = name of file without extension without path
        -logfile  = name of outputfile without path
        -logpath  = full path to outputfile
        '''
        path, filename=filepath.rsplit('/',1)

        self.filepath=filepath
        self.path=path
        self.filename=filename

        if '.' in filename:
            self.name=filename.split('.')[0]
        else: self.name=self.filename

        self.logfile=self.name + '.log'
        self.logpath=self.path + '/' + self.logfile

        self.calc=calc
        self.errorpath = None
        self.errorfile = None
        self.ignorejob = False
        self.IsReady = False
        return

    def __str__(self):
        return "Job: {}".format(self.name)

    def __repr__(self):
        return self.__str__()

    def write(self):
        '''write the job input file'''
        pass

    def submit(self):
        '''submit the job input file'''
        jobid = None
        try:
            if self.calc['nosub'] ==2:
                time.sleep(1)
                subm.nosubmit(self, self.cmd)
                jobid = '0'
            else:
                jobid = subm.submit(self, self.script).strip()
        except Exception as e:
            print repr(e)
            print "retry submit after 5s..."
            time.sleep(5)
            jobid = subm.submit(self, self.script).strip()
        print '{} submitted'.format(self)
        return jobid

    def read(self):
        ''' reads the job output'''
        pass

    def normaltermination(self, debug=True, **kwargs):
        ''' checks if the job has ran correctly '''
        for attempt in range(3):
            # these attempts are made because sometimes
            # files are sometimes not yet copied back.
            try:
                ret = self.termination(self.logpath)
            except IOError as e:
                print "read error with:", self.logpath
                time.sleep(10)
                ret = 3
                continue
            else:
                break

        if ret==3:
            print "open-new-file submit-problem. trying to resubmit"
            self.submit()
        elif not ret==1:
            print "Error termination:", self.logpath
            self.errortermination(debug)
        else:
            self.IsReady=True
        return

    def ready(self, ignore=0, extratime=0, **kwargs):
        ''' if job had no normal termination, this function checks how to proceed:
            1. by a normal termination of the logfile
            2. by a normal termination of the errorlogfile (*zzz.log)
            3. by IGNORE statement
        '''
        # try a normal termination of errorpath
        try:
            ret = self.termination(self.logpath)
        except IOError:
            ret = 0

        if self.errorpath:
            try:
                ret_zzz=self.termination(self.errorpath)
            except IOError:ret_zzz=0
        else: ret_zzz=0

        if ret==1:
            self.IsReady=True
        elif ret==2:
            print "\n{0}\n\t  IGNORED: {1} IGNORED!\n{0}\n".format("    --oOo--"*10, self.logpath)
            self.ignorejob=True
            self.IsReady=True
        elif self.errorpath and ret_zzz==1:
            import shutil
            shutil.copyfile(self.errorpath, self.logpath)
            time.sleep(1)
            print "*zzz.log file with normal termination copied back to original logfile."
            pass # in the following iteration of while true the termination(path) should return 1
        elif (not self.errorpath is None) and ignore and extratime>ignore:
            self.ignorejob=True
            self.IsReady=True
            print "\n\n{0}\n    AUTOMATICALLY IGNORED after {2} seconds of waiting: {1}\n{0}\n".format("    --oOo--"*10, self.logpath, str(ignore))
            try:
                with open(self.logpath, 'a') as f:
                    f.write(' IGNORED')
            except IOError as e:
                print "not able to write IGNORE to file"
        return

