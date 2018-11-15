#!/bin/env python

# python modules
import sys
import numpy as np
import json
import pprint
from copy import deepcopy
import random
import itertools
import logging

# my own modules
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.utils.molecule import Molecule
from CINDES.utils.table import set_table, get_property_table
from CINDES import INDES
from CINDES.INDES.predictions import predictor
from CINDES.predictor import learning_int as ml_i
from CINDES.INDES.GA import Fitness_Function
from CINDES.INDES.construction import indtocon, contoind

from CINDES.pyevolve import DBAdapters

###### Set global variables:
debug = False
pp = pprint.PrettyPrinter(width=200)
np.set_printoptions(linewidth=120)

class MyDBSQLiteAdapter(DBAdapters.DBSQLite):
    '''This class only overwrites the insert method of the normal DBAdapters given in pyevolve'''

    def insert(self, pso_engine):
        """ Inserts the statistics data to database

        :param ga_engine: the GA Engine

        .. versionchanged:: 0.6
           The method now receives the *ga_engine* parameter.
        """
        stats      = pso_engine.getStatistics()
        generation = pso_engine.iter


        c = self.getCursor()

        logging.debug("in insert with:" + repr(stats) + repr(c))

        pstmt = "insert into %s values (?, ?, " % ("statistics")
        for i in xrange(len(stats)):
            pstmt += "?, "
        pstmt = pstmt[:-2] + ")"
        c.execute(pstmt, (self.getIdentify(), generation) + stats)

        pstmt = "insert into %s values(?, ?, ?, ?, ?, ?)" % ("population",)
        tups = []
        for particle in pso_engine.swarm:
            tups.append((self.getIdentify(), generation, particle.index, particle.localbestP, particle.P, contoind(particle.localbestsample)))
        tups.sort(key=lambda x:x[4])

        c.executemany(pstmt, tups)

        if (generation % self.commitFreq == 0):
            self.commit()
        return


###### CALL(s) from __main__.py ###########

# 1. setup system
from CINDES.INDES import procedures

def rounder(history):
    n=2
    roundn = lambda x:round(x, n)
    for i in range(len(history)):
        history[i] = map(roundn, history[i])
    return history

def stringify_array(array):
    newarray = []
    for site in array:
        newsite = []
        for group in site:
            newgroup = "".join(group)
            newsite.append(newgroup)
        newarray.append(newsite)
    return newarray


def main(param):
    global w1, w2, c1

    # -1. random seeds:
    np.random.seed(param['seed'])
    random.seed(param['seed'])

    # 0. setup
    mprms = procedures.FrameRun(**param)
    print mprms
    table = set_table(mprms, mprms.array)
    mprms.array = stringify_array(mprms.array)

    # 1. define fitness function
    FF = Fitness_Function(mprms, table=table, array=mprms.array)

    mypso = run_pso(mprms, function=FF)
    return mypso


def run_pso(mprms, function):
    # 2.
    w1 = mprms.pso['w1']
    w2 = mprms.pso['w2']
    c1 = mprms.pso['c1']

    # 3. Initialize Algorithm:
    minimize = mprms.optimum == 'minimum'
    if mprms.pso['type']=='concrete':
        from concretePSO import ConcretePSO
        mypso = ConcretePSO(array=mprms.array,
                            npop= mprms.pso['npopulation'],
                            function=FF,
                            maxiter=mprms.pso['ngenerations'],
                            parallel=True,
                            minimize=minimize)
    else:
        options = { k:v for k,v in mprms.pso.iteritems() if k in ('w', 'c1', 'c2', 'epsilon') }
        from probabilityPSO import ProbabilityPSO
        mypso = ProbabilityPSO(array=mprms.array,
                            npop =mprms.pso['npopulation'],
                            function = function,
                            maxiter=mprms.pso['ngenerations'],
                            parallel=True,
                            minimize=minimize,
                            options = options)

    # initialize DB:
    sqlite_adapter = MyDBSQLiteAdapter(
            dbname='psostats.db',
            identify=mprms.pso['db_identify'],
            resetDB=False,
            resetIdentify=True,
            commit_freq=1)
    mypso.setDBAdapter(sqlite_adapter)

    # 4. Run Algorithm:
    mypso.evolve()

    # 5. Some Logging
    logging.info("global bestX:" + repr(np.round(mypso.globalbestX, decimals=3)))
    logging.warning("global bestP:" + repr(mypso.globalbestP))

    totalhistory=[]
    for particle in mypso.swarm:
        totalhistory.append(particle.history)
    logging.info("particle history:" + str(rounder(totalhistory)))

    localhistory=[]
    for particle in mypso.swarm:
        localhistory.append(particle.localhistory)
    logging.info("particle.localbest history:" + str(rounder(localhistory)))

    logging.info("global best history:" + str(map(lambda x:round(x,8), mypso.globalhistory)))

    # do only when run on commandline and logging level INFO/DEBUG
    if sys.stdin.isatty() and logging.getLogger().isEnabledFor(logging.INFO):
        import matplotlib.pyplot as plt
        for i, a in enumerate(totalhistory):
            plt.plot(np.array(a)+i*0.05, alpha=0.9)
        plt.show()
        for i, a in enumerate(localhistory):
            plt.plot(np.array(a)+i*0.05, alpha=0.9)
        plt.show()


    return mypso


if __name__ == "__main__":
    pass
