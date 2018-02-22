#!/bin/env python
'''a BaseClass for any kind of Job'''

class BaseJob(object):
    '''BaseJob of any kind of submittable job'''

    def __init__(self, *args, **kwargs):
        return

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


