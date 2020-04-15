#!/bin/env python
debug = False
verbose = False
test = True
# python modules
import numpy as np
random = np.random.random
#import pickle
import json

import pprint

# needed by evolve
from time import time
import logging
from sys import platform as sys_platform
from sys import stdout as sys_stdout
from CINDES.algorithms.pyevolve.GPopulation import GPopulation
from CINDES.algorithms.pyevolve.GPopulation import Util
import random as rrandom
from copy import deepcopy

# my own modules
#from writings import log_io, sprint, print_title
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.molecule import Molecule
from CINDES.utils.table import set_table, get_property_table
from CINDES.utils.utils import skipper
from CINDES.run import FrameRun
from CINDES.algorithms.GA import Fitness_Function
from CINDES.evaluation.predictions import predictor
from CINDES.evaluation.calculator import evaluate_mols
from CINDES.evaluation.construction import indtocon, contoind
from . import procedures

from .deap import tools
from .deap import base
from .deap import creator
from .deap import algorithms
from .DEAP_utils import *


def run_GA(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    raise NotImplementedError('bla')
    history = []
    import numpy as np
    np.random.seed(options.seed)

    from .deap import base
    from .deap import creator
    from .deap import tools
    from .deap import algorithms


    minimize = options.optimum == 'minimize'

    creator.create("FitnessSingle", base.Fitness, weights=options.genalg['weights'])

    # this creates an inidivual class. instances will automatically be 'individual'
    creator.create("Individual", list, fitness=creator.FitnessSingle)

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
    groups_joined = [ "".join(item) for item in groups ]
    get_group = lambda : np.random.choice(groups_joined)
    print("WARNING: assuming all sites have the same group!")
    toolbox.register("get_group", get_group)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
            toolbox.get_group, n=options.nsites)

    # make a population function. now calling: 'toolbox.population(n=3) will return 3 fresh indivuals
    toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=options.genalg['npopulation'])

    def evaluate_multi(population, gen=None):
        print("population:", population)
        P = [ list(map(list,p)) for p in population[:] ]
        mols = function.evaluate_multi(P, gen=gen)

        print("props:", [ mol.props for mol in mols ])

        vals = [ mol.Pvalue for mol in mols ]
        return vals

    def mutateF(individual, MUP):
        newindividual = deepcopy(individual)
        if np.random.random() < MUP:
            i = np.random.randint(0, len(individual))
            newindividual[i] = np.random.choice(groups_joined)
        return (newindividual,)

    # register mutation
    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])
    toolbox.register("select", tools.selTournament, tournsize=3)


    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])

    # make initial population
    population = toolbox.population()
    if options.genalg['restart']:
        startgen = restart_population(population)
        print("start generation at:", startgen)
        print("with population:\n", end=' ')
        for i in population:
            print(i)
    else:
        startgen = 1



    # evaluate initial population
    fits = evaluate_multi(population)
    for fit, ind in zip(fits, population):
        try:
            ind.fitness.values = [fit]
        except TypeError as e:
            print("fit:", fit)

    for gen in range(startgen, options.genalg['ngenerations']):
        print_title("Generation No.: " + str(gen), outline='l', signator="-")

        # 0.
        offspring = toolbox.select(population + halloffame[:], 
                len(population))

        # 1. new generation
        offspring = algorithms.varAnd(offspring, toolbox,
            cxpb=options.genalg['CXP'],
            mutpb=options.genalg['MUP'])

        # 2. evaluate fitness
        fits = evaluate_multi(offspring, gen=gen)
        print("fits:", fits)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = [fit]

        # 3. loggings
        halloffame.update(offspring)
        print("hall of fame:", halloffame)

        population[:] = offspring

        record = stats.compile(population)
        print("record:", record)
        logbook.record(gen=gen, **record)
        best_ind = tools.selBest(population, 1)[0]
        history.append({
            'gen':gen,
            'ncalcs':function.ncalcs,
            'index':contoind(best_ind),
            'p':best_ind.fitness.values[0]
            })

        # 4. select new population
        population = toolbox.select(offspring + halloffame[:],
            k=options.genalg['npopulation'])

    print("hall of fame:", halloffame)
    print("population:", population)
    print("logbook:", logbook)
    best_ind = tools.selBest(population, 1)[0]
    best_index = contoind(best_ind)
    print("Best individual is %s, %s" % (best_index, best_ind.fitness.values[0]))
    print("ncalcs:", function.ncalcs)

    result = {'p':best_ind.fitness.values[0], 'index':best_index, 'history':history, 'ncalcs':function.ncalcs, 'gen':options.genalg['ngenerations']}

    return result

def run_test_GA(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    history = []
    import numpy as np
    #np.random.seed(options.seed)

    from .deap import base
    from .deap import creator
    from .deap import tools
    from .deap import algorithms


    minimize = options.optimum == 'minimize'

    creator.create("FitnessSingle", base.Fitness, weights=options.genalg['weights'])

    # this creates an inidivual class. instances will automatically be 'individual'
    creator.create("Individual", list, fitness=creator.FitnessSingle)

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
    groups_joined = [ "".join(item) for item in groups ]
    get_group = lambda : np.random.choice(groups_joined)
    print("WARNING: assuming all sites have the same group!")
    toolbox.register("get_group", get_group)
    toolbox.register("individual", tools.initRepeat, creator.Individual,
            toolbox.get_group, n=options.nsites)

    # make a population function. now calling: 'toolbox.population(n=3) will return 3 fresh indivuals
    toolbox.register("population", tools.initRepeat, list, toolbox.individual, n=options.genalg['npopulation'])

    def mutateF(individual, MUP):
        newindividual = deepcopy(individual)
        if np.random.random() < MUP:
            i = np.random.randint(0, len(individual))
            newindividual[i] = np.random.choice(groups_joined)
        return (newindividual,)

    # register mutation
    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])
    toolbox.register("select", tools.selTournament, tournsize=3)

    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])

    # make initial population
    population = toolbox.population()

    toolbox.register("evaluate", options.function)
    population, logbook, history = eaSimple(
            population,
            toolbox,
            cxpb=options.genalg['CXP'],
            mutpb=options.genalg['MUP'],
            ngen=options.genalg['ngenerations'],
            stats=stats,
            halloffame=halloffame,
            verbose=False,
            nelitism=options.genalg['nelitism'])

    print("hall of fame:", halloffame)
    if verbose:
        print("population:", population)
        print("logbook:", logbook)
    best_ind = tools.selBest(population, 1)[0]
    best_index = contoind(best_ind)
    print("Best individual is %s, %s" % (best_index, best_ind.fitness.values[0]))
    print("ncalcs:", function.ncalcs)

    result = {'p':best_ind.fitness.values[0], 'index':best_index, 'history':history, 'ncalcs':history[-1]['ncalcs'], 'gen':options.genalg['ngenerations']}

    return result

