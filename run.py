#!/bin/env python
''' This file contains a parameter class called BaseRun and FrameRun which holds all the input parameters

'''

# this line must be at the beginning of the file!
from __future__ import division

# debug flag
debug = 1

# import python libraries
import time  # for getting time/date and time delays
start = time.clock()
from inspect import stack
from math import ceil
import shutil  # module to copy files
from platform import node
from pprint import pformat
import os  # for getting window width and testing existence of files
import inspect # to see if function is a class
import logging  # instead of the large amount of print statements not using it at the moment
from copy import deepcopy  # for keeping matrices while changing others
import numpy as np

# import my own modules
from evaluation import construction as zcon  # all functions needed for constructing new geometries
from evaluation import reader as r  # this reads the zmatrix in gaussian format

# import utils
from CINDES.utils.writings import dump
from CINDES.utils.utils import *
from CINDES.evaluation.job import BaseJob

print "time for imports:", time.clock() - start

class BaseRun(object):
    ''' This is the main object for all the parameters used during any process
    this object is initiated with a dictionary from the inputreader '''

    def __init__(self, **entries):
        self.__dict__.update(entries)  # here all the key/value pairs in entries are converted to attributes.

        # give some logging information:
        logging.debug("loggers:" + pformat(logging.Logger.manager.loggerDict))
        logging.debug("rootlogger level" + pformat(logging.getLogger().getEffectiveLevel()))
        logging.debug("rootlogger handlers:" + pformat(logging.getLogger().handlers))

        # set system variables
        self.script = stack()[0][1]
        self.node = node()
        self.starttime = time.time()
        self.directory = os.getcwd()
        self.pid = os.getpid()
        self.ppid = os.getppid()

        if not self.nosub==2:
            self.setup_filesystem()
            self.set_calcs()

        # for Gaussian this is default. otherwise it has to be switched on
        if self.assign_geom==True or self.program=='gaussian':
            BaseJob.assign_geom=True

        # sometimes complicated property functions have to be initialized:
        if self.property=='func' and inspect.isclass(self.function):
            print "initializing function..."
            self.function = self.function(self)


    def __str__(self):
        sb = ['BaseRun object with the following attributes:']
        empty_attributes = []
        for key, value in sorted(self.__dict__.items()):
            try:
                if not value:  # i.e. value is either None, False, zero, empty list/string
                    empty_attributes.append(key)
                    continue
            except ValueError:
                pass
            if key=='array':continue
            if key in ['predictions']:
                sb.append("{key:20}=".format(key=key))
                sb.append(dump(value))
            elif key in ['TZmat', 'TZmatrices', 'genalg', 'jobs', 'stabjobs', 'prejobs', 'extrajobs', 'calcs', 'pso']:
                sb.append("{key:20}=".format(key=key))
                sb.append(pformat(value, width=150))
            elif key == 'function' and callable(value):  # i.e. the value is a lambda function
                if value.__doc__:
                    sb.append("{key:20}={value}".format(key=key, value=value.__doc__))
                else:
                    sb.append("{key:20}= lambda function".format(key=key))
            elif key in ['adj']:
                sb.append("{key:20}=".format(key=key))

                def f(v): return ''.join([('0', '1')[int(item)] for item in v])
                sb.append('\n'.join(map(f, value)))
            else:
                sb.append("{key:20}= {value}".format(key=key, value=value))

        # now a table of items that are empty is printed in rows of 5 items
        l = 5
        inrows = [ empty_attributes[l*i:l*i+l] for i in range(int(ceil(len(empty_attributes)/l))) ]
        asstr = "\n".join([ "".join(map(" {:20}".format, item)) for item in inrows])
        sb.append("empty attributes    =\n{}".format(asstr))
        return '\n'.join(sb)

    def __repr__(self):
        return self.__str__()

    # the the calculationskeyword:
    def set_calcs(self):
        def tocalc(calc):
            # path only set after setup_filesystem!
            tohavekeys = ['program', 'nprocs', 'identify', 'path', 'nosub']
            for key in tohavekeys:
                if not key in calc:
                    calc[key] = paras[key]
                elif key=='identify':
                    identify = self.identify.rstrip('_') + calc[key]
                    if not identify[-1]=='_':
                        identify += '_'
                    calc[key]=identify
            return calc

        def check(cal, i):
            if cal['identify'] in identifiers:
                cal['identify'] = "{}{}_".format(cal['identify'], str(i))
                assert not cal['identify'] in identifiers
                i += 1
            identifiers.append(cal['identify'])
            return cal, i

        paras = self.__dict__

        # if no calcs are given
        if not (paras['jobs'] or paras['calcs']):
            self.calcs = []
            return

        # if calcs are given as jobs input:
        if paras['jobs']:
            calc = paras['jobs']
            calcs = [tocalc(calc)]
        else:
            # we modify the jobs to calcs this is a bit confusing since calcs are actually
            # still jobs
            precalcs = paras['calcs']
            calcs = []
            for calc in precalcs:
                if isinstance(calc, list):
                    C = []
                    for cal in calc:
                        C.append(tocalc(cal))
                    calcs.append(C)
                else:
                    calcs.append(tocalc(calc))
        if paras['prejobs']:
            calc = paras['prejobs']
            calcs.insert(0, tocalc(calc))
        for key in ['stabjobs', 'extrajobs']:
            if paras[key]:
                calc = paras[key]
                calcs[-1] = [calcs[-1], tocalc(calc)]

        # verify that there are no similar identifiers
        identifiers = []
        i = 1
        for calc in calcs:
            if isinstance(calc, list):
                for cal in calc:
                    cal, i = check(cal, i)
            else:
                calc, i = check(calc, i)
        self.calcs = calcs
        return

    def setup_filesystem(self):
        param = self.__dict__
        if self.nosub == 2:
            path = ''
        else:
            #param['workdir'] = os.getcwd()
            self.workdir = os.getcwd()
            path = self.workdir + '/CALC'
            param["path"] = str(path)
            if not os.path.exists(path):
                os.makedirs(path)
            if param['program'] == 'gaussian':
                self.script = 'ID_gauss'
                self.extension = '.com'
                shutil.copy(os.getcwd() + '/ID_gauss', path)
            elif param['program'] == 'orca':
                self.script = 'ID_orca'
                shutil.copy(os.getcwd() + '/ID_orca', path)
            elif param['program'] == 'nwchem':
                self.script = 'ID_NWChem'
                self.extension = ''
                shutil.copy(os.getcwd() + '/ID_NWChem', path)
            elif param['program'] == 'vasp':
                self.script = 'ID_VASP'
                self.extension = ''
                shutil.copy(os.getcwd() + '/ID_VASP', path)
            else:
                raise SystemExit('ERROR: No valid program specified')
        self.path = path
        return param, path

    def currenttime(self):
        return "Current time %s" % str(time.time() - self.starttime)

class XYZRun(BaseRun):

    def __init__(self, **entries):
        super(XYZRun, self).__init__(**entries)
        if self.nsites and not self.nosub==2:
            self.cartesian = r.get_cartesian(**entries)
        else:
            logging.info("NO ZMAT NOR XYZ")
        return


class FrameRun(BaseRun):
    ''' This inherites from BaseRun and is the main object for all BFS/SD molecular frame based
    procedures.
    '''

    def __init__(self, **entries):
        super(FrameRun, self).__init__(**entries)
        # specific for FrameRun:
        if self.nsites and not self.nosub==2:
            if isinstance(self.zmatrixfile, list):
                self.TZmatrices={}
                for zmatrixfile in self.zmatrixfile:
                    tzmat = r.geometry(zmatfile=zmatrixfile, **entries)
                    self.TZmatrices[zmatrixfile]=tzmat
                self.TZmat = self.TZmatrices[self.zmatrixfile[0]]
            else:
                self.TZmat = r.geometry(zmatfile=self.zmatrixfile, **entries)
            self.adj = self.set_adj(self.TZmat['core'], self.TZmat['active'])
            self.corresp = self.set_corresp(self.TZmat['active'], self.TZmat['passive'])
        else:
            logging.info("NO ZMAT")
        return

    def set_adj(self, core, active):
        debug = 0
        from CINDES.utils.converter import Converter
        import numpy as np
        conv = Converter()
        conv.read_zmalist(core)
        xyz = np.asarray([atom[1] for atom in conv.zmatrix_to_cartesian()])
        #print xyz
        ncore = len(xyz)
        adj = np.zeros([ncore, ncore])
        for i in range(ncore):
            for j in range(i, ncore):
                adj[i][j] = 0.1 < np.linalg.norm(xyz[i] - xyz[j]) < 2.0
                adj[j][i] = adj[i][j]
        sites = [int(methyl[0][1]) - 1 for methyl in active]
        sites_adj = adj[sites][:, sites]
        if debug:
            print "adjacency matrix of core:", adj
            print "self.sites:", self.sites
            print "active: ", active
            print "sites: ", sites
            print "sites_adj:", sites_adj
        return sites_adj

    def set_corresp(self, active, passive):
        '''makes a dictionary that gives the correspondance of sites with position in core matrix'''
        corresp = dict()
        correspI= dict()
        for site in active:
            corresp[int(site[0][1])] = int(site[1][1])
            correspI[int(site[1][1])] = int(site[0][1])
        for site in passive:
            corresp[int(site[0][1])] = int(site[1][1])
            correspI[int(site[1][1])] = int(site[0][1])
        return correspI

