#!/bin/env python
'''a BaseClass for any kind of Job'''

class BaseJob(object):
    '''BaseJob of any kind of submittable job'''

    def __init__(self, filepath, calc=dict()):
        '''There are different combinations of path and file to be set:
        -filepath = full path to INPUT file
        -path     = path where file is located
        -filename = name of inputfile with extension
        -name     = name of file without extension
        -logfile  = name of outputfile
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
        print "self.path:", self.path
        print "self.filename:", self.filename
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
        pass

    def read(self):
        ''' reads the job output'''
        pass

    def check_normal_termination(self):
        ''' checks if the job has ran correctly '''



