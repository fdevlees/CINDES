#!/bin/env python
''' this module contains all functions related to the NWChem program '''
import re
import time
import string
import random
import os
multiplicity = {1: 'singlet', 2: 'doublet', 3: 'triplet', 4: 'quartet'}

from .job import BaseJob


class VASPJob(BaseJob):
    extension = ''
    script = 'ID_VASP'
    cmd = './vasp'

    def write(self):
        pass
    # -----

    def termination(self, logpath):
        try:
            with open(logpath, 'r') as fid:
                text = fid.readlines()[-4:]
                raise NotImplementedError('VASP determinations not yet implemented')
        except IOError:
            return 0
    # -----

    def errortermination(self, debug=False):
        raise NotImplementedError('VASP error handling not implemented')


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(mol, calc, pos=None):
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    # ------------
    # this function uses globals: identify, path
    # ------------
    index = mol.index
    paras = calc
    jobs = paras['jobs']

    name = paras['identify'] + str(index)
    filename = name
    filepath = paras['path'] + '/' + filename + '/POSCAR'
    try:
        geom = calc['geom']
    except KeyError:
        geom = None
    Job = VASPJob(filepath, calc)
    mol.addjob(Job)

    fid = open(filepath, 'w')

    # JOB 1
    # extract info
    job1 = jobs[0]

    # write info
    fid.write('write job_info\n')
    for xyz in mol.cartesian:
        fid.write("{:2s} {:12.6f} {:12.6f} {:12.6f}\n".format(xyz[0], *xyz[1]))

    fid.close()
    return Job


