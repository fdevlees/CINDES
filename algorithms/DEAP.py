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
from CINDES.algorithms.pyevolve.GPopulation import GPopulation
from CINDES.algorithms.pyevolve.GPopulation import Util
import random as rrandom

# my own modules
#from writings import log_io, sprint, print_title
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.utils.molecule import Molecule
from CINDES.utils.table import set_table, get_property_table
from CINDES.utils.utils import skipper
from CINDES.run import FrameRun
from CINDES.algorithms.GA import Fitness_Function
from CINDES.evaluation.predictions import predictor
from CINDES.evaluation.calculator import evaluate_mols
from CINDES.evaluation.construction import indtocon
import procedures


def main(run):
    print run
    table = set_table(run, run.array)
    function = Fitness_Function(run, table=table, array=run.array)
    final_genome = run_deap(run.array, run, function=function)
    print "final genome"
    return final_genome


def run_deap(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    import numpy as np
    np.random.seed(options.seed)

    from deap import base
    from deap import creator
    from deap import tools
    from deap import algorithms

    minimize = options.optimum == 'minimize'


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

    toolbox = base.Toolbox()

    # a function that gives a group. This function is called by initRepeat
    groups = options.array[0]
    get_group = lambda : np.random.choice(groups)
    print "WARNING: assuming all sites have the same group!"
    toolbox.register("get_group", get_group)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
            toolbox.get_group, n=options.nsites)

    # make a population function. now calling: 'toolbox.population(n=3) will return 3 fresh indivuals
    toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=options.genalg['npopulation'])

    def evaluate_multi(population, gen=None):
        P = [ list(p) for p in population[:] ]
        mols = function.evaluate_multi(P, gen=gen)

        print "props:", [ mol.props for mol in mols ]

        vals = [ mol.Pvalue for mol in mols ]
        return vals

    def mutateF(individual, MUP):
        if np.random.random() < MUP:
            i = np.random.randint(0, len(individual))
            individual[i] = np.random.choice(groups)
        return (individual,)

    # register mutation
    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])

    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])
    toolbox.register("select", tools.selNSGA2)

    # make initial population
    population = toolbox.population()
    if options.genalg['restart']:
        startgen = restart_population(population)
        print "start generation at:", startgen
        print "with population:\n",
        for i in population:
            print i
    else:
        startgen = 1

    # evaluate initial population
    fits = evaluate_multi(population)
    for fit, ind in zip(fits, population):
        ind.fitness.values = fit

    for gen in range(startgen, options.genalg['ngenerations']):
        print_title("Generation No.: " + str(gen), outline='l', signator="-")

        # 1. new generation
        offspring = algorithms.varOr(population, toolbox,
            lambda_=options.genalg['npopulation'],
            cxpb=options.genalg['CXP'],
            mutpb=options.genalg['MUP'])

        # 2. evaluate fitness
        fits = evaluate_multi(offspring, gen=gen)
        print "fits:", fits
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit

        # 3. loggings
        halloffame.update(population)
        print "hall of fame:", halloffame

        record = stats.compile(population)
        print "record:", record
        logbook.record(gen=gen, **record)

        # 4. select new population
        population = toolbox.select(offspring + population,
            k=options.genalg['npopulation'])

    print population
    print logbook

    return population

def restart_population(population):

    with open('cyclesinfo') as f:
        data = [ line.strip().split() for line in f.readlines() ]

    # get generation of last compound. always index -3
    maxgen = int(data[-1][-3])

    confs = []

    i=0
    while True:
        # loop backwards getting indices -1, -2, -3 etc while i is going from 0, 1, 2, etc
        newline = data[-(i+1)]

        # take last lines until maxgen changes
        if not int(newline[-3])==maxgen:
            break

        # get index from line and get as conf
        conf = indtocon(newline[0].strip("'"))

        confs.append(conf)

        i+=1

    # check if same size
    if not len(population) == len(population):
        raise NotImplementedError('length of old and new population do not match!')

    # set old confs as new confs
    for individual, conf in zip(population, confs):
        individual[:] = conf

    return maxgen




if __name__ == "__main__":
    pass
