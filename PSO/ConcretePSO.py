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
from predictions import predictor
from CINDES.predictor import learning_int as ml_i
from CINDES.INDES.GA import Fitness_Function
from CINDES.INDES.construction import indtocon, contoind

###### Set global variables:
debug = False
pp = pprint.PrettyPrinter(width=200)
np.set_printoptions(linewidth=120)

# weight
w1, w2 = None, None
c1 = None

###### The algorithm:


def rounder(history):
    n=2
    roundn = lambda x:round(x, n)
    for i in range(len(history)):
        history[i] = map(roundn, history[i])
    return history

def permutate(array):
    nrows, ncols = array.shape
    print "ncols:", ncols
    raise NotImplementedError('this is too much!')
    perms = np.array(list(itertools.permutations(range(ncols))))
    choices = np.random.randint(len(perms), size=nrows)
    i = np.arange(nrows).reshape(-1, 1)
    return array[i, perms[choices]]

def permutate(array):
    array = np.array(array)
    nrows, ncols = array.shape
    print "ncols:", ncols
    for i in range(len(array)):
        np.random.shuffle(array[i])
    return array


def priority_score(dG, dL, dA, progress):
    if progress<=0.5:
        score = (1.0 - 4./3.*progress) * float(dA)
        score -= w1 * (2/3.*progress) * float(dL)
    else:
        score = 1./3. * float(dA)
        score -= w1 * (2./3. - 2/3*progress) * float(dL)
    score -= w2 * (2./3.* progress) * float(dG)
    #score = (1.0-progress)*dA - w * progress * dG
    # test only converge to global minimum
    #print "score:", score,
    return score + c1*random.random()


def distance(X1, X2):
    # if two indices are pased: calculate the number of not shared groups
    score=0
    for i in range(len(X1)):
        if not X1[i]==X2[i]:
            score+=1
    return score


def average_distance(DB, X):
    total_distance = 0.0
    for x in DB:
        total_distance += distance(x, X)
    return total_distance / len(DB)


class Particle(object):

    def __init__(self, X, index=None):
        self.X = X
        self.P = None
        self.fitness = None
        self.index = index
        self.history = []
        self.nDim = len(X)

        self.localbestX = X
        self.localbestP = None
        self.localhistory = []
        return

    def __repr__(self):
        return "<P{}: {}>".format(self.index, repr(self.X))

    def __getitem__(self, i):
        return self.X[i]

    def evaluate(self, function):
        self.P = function(self.X)
        self.history.append(self.P)
        return

    def setproperty(self, P):
        self.P = P
        self.history.append(self.P)
        return

    def update_local_best(self):
        self.localbestP = deepcopy(self.P)
        self.localbestX = deepcopy(self.X)

    def updateX(self, globalbestX, array, progress, DB):
        # here do a new mutant
        # for initial test just mutate one site to the globalbestX
        mutants = self.gen_single_mutants(array)
        #print "self.X:", self.X

        # for every mutant evaluate priority
        best_P  = None
        best_X  = None
        for mutant in mutants:
            dG = distance(globalbestX, mutant)
            dL = distance(self.localbestX, mutant)
            dA = average_distance(DB, mutant)
            score = priority_score(dG, dL, dA, progress)
            if score >= best_P or best_P is None:
                if score == best_P:
                    if np.random.random()<0.5:continue
                best_P = score
                best_X = deepcopy(mutant)
            else:
                pass
                #print
        #print "best_P", best_P
        #print "best_X", best_X
        self.X = best_X
        return


    def gen_single_mutants(self, array):
        mutants = []
        for i in range(self.nDim):
            for group in array[i]:
                if not self.X[i]==group:
                    mutant = deepcopy(self.X)
                    mutant[i]=group
                    mutants.append(mutant)
        return mutants


class ConcretePSO(object):
    def __init__(self, array, npop, function,
            parallel=False,
            maxiter=10,
            minimize=True,
            options=None):

        # initiate variables
        self.iter = 0
        self.n = npop
        self.array = array
        self.nDim = len(array)
        self.function = function
        self.parallel = parallel
        self.maxiter = maxiter
        self.minimize = minimize
        self.globalbestX = None
        self.globalbestP = None
        self.globalhistory = []

        self.DB = []

        # initiate swarm with positions and velocities
        self.init_swarm()

        return

    def __repr__(self):
        return "<PSO object at {}>".format(id(self))

    def __getitem__(self, index):
        return self.swarm[index]

    def init_concrete_swarm(self):
        shuffles=[]
        ndim = len(self.array)
        maxdim = max(map(len, self.array))
        div = list(divmod(self.n, maxdim)) # returns (n, rest)
        if div[1]: div[0]+=1 # if there is a rest an extra shuffle is needed
        print "nshuffle", div[0]
        for i in range(div[0]):
            # generate a randomly ordered array. we do not need permutate only shuffle!
            newarray = map(list, zip(*permutate(np.array(deepcopy(self.array)))))
            print "newarray:", newarray
            shuffles.append(newarray)
        print "shuffles:", shuffles
        swarm = np.concatenate(shuffles)[:self.n]
        print "swarm:", swarm
        return map(list,swarm)

    def init_swarm(self):
        concrete_swarm = self.init_concrete_swarm()
        swarm = []
        for idx, X in enumerate(concrete_swarm):
            particle = Particle(X, index=idx)
            swarm.append(particle)
        self.swarm = swarm
        return

    def update_global_best(self, particle):
        # set globalbest Property
        self.globalbestP = deepcopy(particle.P)
        self.globalbestX = deepcopy(particle.X)
        return

    def set_fitness(self):
        ''' This is sigma truncation scaling! But do I need it?'''
        c = 2.0
        Ps = [ particle.P for particle in self.swarm ]
        pop_rawAve = np.average(Ps)
        pop_rawDev = np.std(Ps)
        for i in xrange(len(self.swarm)):
            f = self.swarm[i].P - pop_rawAve
            f+= c * pop_rawDev
            if f < 0:
                f = 0.0
            self.swarm[i].fitness = f

    @log_io()
    def evaluate(self):
        ''' this is my evaluate function '''
        # 1. get configurations to calculate:
        populationlist = []
        for particle in self.swarm:
            populationlist.append(particle.X)
        print "populationlist", populationlist, "len:", len(populationlist)

        # 1.1. make confs hashable to make it a set and make it list again
        new_confs = tuple( '_'.join(item) for item in populationlist )
        print "new_confs:", new_confs
        unique_confs = [ indtocon(item) for item in set(new_confs)]
        print "n unique_confs:", len(unique_confs)
        print "unique_confs:", unique_confs

        # 2. call CINDES via FF to calculate the configurations
        mols = self.function.evaluate_multi(unique_confs, gen=self.iter)

        # 3. set the calculations to the correct indivual score
        y_dict = {mol.index: mol.Pvalue for mol in mols}
        print "y_dict:", y_dict
        for particle in self.swarm:
            index = "_".join(particle.X)
            particle.setproperty(y_dict[index])
        return

    def evolve(self):
        self.iter=1
        while True:
            print "----------- Generation {} -----------".format(self.iter)

            # evaluate and set global best
            if self.parallel:
                self.evaluate()
            else:
                for particle in self.swarm:
                    particle.evaluate(self.function)

            for particle in self.swarm:
                if not particle.X in self.DB:
                    self.DB.append(particle.X)

                if self.minimize:
                    if particle.P <= particle.localbestP or particle.localbestP is None:
                        particle.update_local_best()
                    if particle.P <= self.globalbestP or self.globalbestP is None:
                        self.update_global_best(particle)
                else:
                    if particle.P >= particle.localbestP or particle.localbestP is None:
                        particle.update_local_best()
                    if particle.P >= self.globalbestP or self.globalbestP is None:
                        self.update_global_best(particle)

            # since the score is only calculated based on similarity no scaling needs to be done at the moment
            #self.set_fitness()

            # log values
            # to screen
            print "globalbest:", self.globalbestX, self.globalbestP
            # and to history
            self.globalhistory.append(self.globalbestP)
            for particle in self.swarm:
                particle.localhistory.append(particle.localbestP)

            # mutate X
            progress = float(self.iter) / float(self.maxiter)
            for particle in self.swarm:
                particle.updateX(self.globalbestX, self.array, progress, self.DB)

            # decide if last iteration    
            self.iter+=1
            if self.iter>self.maxiter:
                break
        return


