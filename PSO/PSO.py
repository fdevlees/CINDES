#!/bin/env python

# python modules
import sys
import numpy as np
import json
import pprint
from copy import deepcopy
import random
import itertools

# my own modules
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.utils.molecule import Molecule
from CINDES.utils.table import set_table, get_property_table
from CINDES import INDES
from CINDES.INDES.predictions import predictor
from CINDES.predictor import learning_int as ml_i
from CINDES.INDES.GA import Fitness_Function
from CINDES.INDES.construction import indtocon, contoind

###### Set global variables:
debug = False
pp = pprint.PrettyPrinter(width=200)
np.set_printoptions(linewidth=120)

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
                            function = FF,
                            maxiter=mprms.pso['ngenerations'],
                            parallel=True,
                            minimize=minimize,
                            options = options)
    print mypso
    #print mypso.swarm

    # 4. Run Algorithm:
    mypso.evolve()

    # 5. Some Logging
    print "global bestX:", mypso.globalbestX
    print "global bestP:", mypso.globalbestP

    print "particle history:"
    totalhistory=[]
    for particle in mypso.swarm:
        totalhistory.append(particle.history)
    print(rounder(totalhistory))

    print "particle.localbest history:"
    localhistory=[]
    for particle in mypso.swarm:
        localhistory.append(particle.localhistory)
    print(rounder(localhistory))

    print "global best history:"
    print(map(lambda x:round(x,8), mypso.globalhistory))

    if sys.stdin.isatty():
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
