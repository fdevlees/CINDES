#!/bin/env python
debug = False
verbose = False
# python modules
import numpy as np
random = np.random.random
#import pickle
import json
from itertools import izip
import pprint

# needed by evolve
from time import time
import logging
from sys import platform as sys_platform
from sys import stdout as sys_stdout
from CINDES.pyevolve.GPopulation import GPopulation
from CINDES.pyevolve.GPopulation import Util
import random as rrandom

# my own modules
#from writings import log_io, sprint, print_title
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.utils.molecule import Molecule
from CINDES.utils.table import set_table, get_property_table
from CINDES.utils.utils import skipper
from CINDES import INDES
from CINDES.INDES.run import FrameRun
from CINDES.INDES.GA import Fitness_Function
from predictions import predictor
from calculator import evaluate_mols
from CINDES.INDES import procedures


def main(param, array=None):
    if array is None:
        array = param['array']
    GArun = FrameRun(**param)
    print GArun
    table = set_table(GArun, array)
    function = Fitness_Function(GArun, table=table, array=array)
    final_genome = run_deap(array, GArun, function=function)
    print "final genome"
    return final_genome


def run_deap(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    import numpy as np
    np.random.seed(1)

    from deap import base
    from deap import creator
    from deap import tools
    from deap import algorithms

    minimize = options.optimum == 'minimize'

    gen = 0

    #creator.create("FitnessMulti", base.Fitness, weights=(-1.0,-1.0))
    creator.create("FitnessMulti", base.Fitness, weights=options.genalg['weights'])
    # this creates an inidivual class. instances will automatically be 'individual'
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    logbook = tools.Logbook()
    halloffame = tools.HallOfFame(maxsize=1)
    paretofront= tools.ParetoFront()

    toolbox = base.Toolbox()

    # a function that gives a group. This function is called by initRepeat
    groups = options.array[0]
    groups_joined = [ "".join(item) for item in groups ]
    get_group = lambda : np.random.choice(groups_joined)
    print "WARNING: assuming all sites have the same group!"
    toolbox.register("get_group", get_group)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
            toolbox.get_group, n=options.nsites)

    # make a population function. now calling: 'toolbox.population(n=3) will return 3 fresh indivuals
    toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=options.genalg['npopulation'])

    def evaluate_multi(population, gen=None):
        print "population:", population
        P = [ map(list,p) for p in population[:] ]
        mols = function.evaluate_multi(P, gen=gen)

        print "props:", [ mol.props for mol in mols ]

        vals = [ mol.Pvalue for mol in mols ]
        return vals

    def mutateF(individual, MUP):
        if np.random.random() < MUP:
            i = np.random.randint(0, len(individual))
            individual[i] = np.random.choice(groups_joined)
        return (individual,)

    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])

    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])
    toolbox.register("select", tools.selNSGA2)

    population = toolbox.population()
    fits = evaluate_multi(population)
    for fit, ind in zip(fits, population):
        try:
            ind.fitness.values = fit
        except TypeError as e:
            print "fit:", fit

    for gen in range(1, options.genalg['ngenerations']):
        print_title("Generation No.: " + str(gen), outline='l', signator="-")

        offspring = algorithms.varOr(population, toolbox,
            lambda_=options.genalg['npopulation'],
            cxpb=options.genalg['CXP'],
            mutpb=options.genalg['MUP'])

        fits = evaluate_multi(offspring, gen=gen)

        print "fits:", fits
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit
        halloffame.update(population)
        paretofront.update(population)
        print "hall of fame:", halloffame
        record = stats.compile(population)
        print "record:", record
        logbook.record(gen=gen, **record)

        population = toolbox.select(offspring + population,
            k=options.genalg['npopulation'])

    print "pareto front", paretofront
    print population
    print logbook

    return population


if __name__ == "__main__":
    pass
