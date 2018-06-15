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
        try:
            if self.calc['nosub'] ==2:
                time.sleep(1)
                subm.nosubmit(self)
                jobid = '0'
            else:
                jobid = subm.submit(self, self.script).strip()
        except Exception as e:
            print repr(e)
            print "retry submit..."
            subm.submit(self, self.script).strip()
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
                time.sleep(10)
            else:
                break
            finally:
                if ret==3 and once==0:
                    print "open-new-file submit-problem. trying to resubmit"
                    import submitter
                    submitter.submit(self)
                elif not ret==1:
                    print "Error termination:", self.logpath
                    self.errortermination(debug)
                else:
                    self.IsReady=True
        else:
            raise e
        return

    def ready(self, ignore=0, **kwargs):
        ''' if job had no normal termination, this function checks how to proceed:
            1. by a normal termination of the logfile
            2. by a normal termination of the errorlogfile (*zzz.log)
            3. by IGNORE statement
        '''
        pass


