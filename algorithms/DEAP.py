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
from .DEAP_GA import *
from .DEAP_utils import *

def main(run):
    print(run)
    table = set_table(run, run.array)
    function = Fitness_Function(run, table=table, array=run.array)
    print("run.genalg['algorithm']:", run.genalg['algorithm'], 'test:', test)
    if run.genalg['algorithm']=='GA':
        if test:
            final_genome = run_test_GA(run.array, run, function=function)
        else:
            final_genome = run_GA(run.array, run, function=function)
    else:
        final_genome = run_nsgaii(run.array, run, function=function)
    print("final genome")
    return final_genome


def run_nsgaii(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    import numpy as np
    np.random.seed(options.seed)

    history = []

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
    paretofront= tools.ParetoFront()

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
        for i in range(len(individual)):
            if np.random.random() < MUP:
                newindividual[i] = np.random.choice(groups_joined)
        return newindividual,

    # register mutation
    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])

    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])
    toolbox.register("select", tools.selNSGA2)

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
            ind.fitness.values = fit
        except TypeError as e:
            print("fit:", fit)

    for gen in range(startgen, options.genalg['ngenerations']):
        print_title("Generation No.: " + str(gen), outline='l', signator="-")

        # 1. new generation
        offspring = algorithms.varOr(population, toolbox,
            lambda_=options.genalg['npopulation'],
            cxpb=options.genalg['CXP'],
            mutpb=options.genalg['MUP'])

        # 2. evaluate fitness
        fits = evaluate_multi(offspring, gen=gen)
        print("fits:", fits)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit

        # 3. loggings
        halloffame.update(population)
        paretofront.update(population)
        print("hall of fame:", halloffame)

        record = stats.compile(population)
        print("record:", record)
        logbook.record(gen=gen, **record)
        best_ind = selBestSum(population, 1)[0]
        print("best_ind:", contoind(best_ind), "fitness values:", best_ind.fitness)
        history.append({
            'gen':gen,
            'ncalcs':function.ncalcs,
            'index':contoind(best_ind),
            'p':best_ind.fitness.values[0]
            })

        # 4. select new population
        population = toolbox.select(offspring + population,
            k=options.genalg['npopulation'])

    print("pareto front", paretofront)
    print(population)
    print(logbook)
    best_index = contoind(best_ind)

    result = {'p':best_ind.fitness.values[0], 'index':best_index, 'history':history, 'ncalcs':history[-1]['ncalcs'], 'gen':options.genalg['ngenerations']}

    return result


def run_nsgaii(array, options, level=None, function=None):
    np.random.seed(options.seed)

    history = []
    best = None

    minimize = options.optimum == 'minimize'
    creator.create("FitnessMulti", base.Fitness, weights=options.genalg['weights'])

    # this creates an inidivual class. instances will automatically be 'individual'
    creator.create("Individual", list, fitness=creator.FitnessMulti)
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
        for i in range(len(individual)):
            if np.random.random() < MUP:
                newindividual[i] = np.random.choice(groups_joined)
        return newindividual,

    # register mutation
    toolbox.register("mutate", mutateF, MUP=options.genalg['MUP'])

    # uniform crossover from tools
    toolbox.register("mate", tools.cxUniform, indpb=options.genalg['CXP'])
    toolbox.register("select", tools.selNSGA2)
    #toolbox.register("select", tools.selSPEA2)

    MU = 100

    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean, axis=0)
    stats.register("std", np.std, axis=0)
    stats.register("min", np.min, axis=0)
    stats.register("max", np.max, axis=0)
    
    logbook = tools.Logbook()
    logbook.header = "gen", "evals", "std", "min", "avg", "max"
    
    pop = toolbox.population(n=options.genalg['npopulation'])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in pop if not ind.fitness.valid]

    #for ind, fit in zip(invalid_ind, fitnesses):
    #    ind.fitness.values = fit
    #fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    fits = evaluate_multi(invalid_ind)
    for fit, ind in zip(fits, invalid_ind):
        if best is None or sum(best[1]) > sum(fit):
            best = (ind, fit)
        try:
            ind.fitness.values = fit
        except TypeError as e:
            print("fit:", fit)

    # This is just to assign the crowding distance to the individuals
    # no actual selection is done
    pop = toolbox.select(pop, len(pop))
    best_ind = selBestSum(pop, 1)[0]
    print("pop type:", type(pop[0]))
    print("best type:", type(best_ind))
    
    record = stats.compile(pop)
    logbook.record(gen=0, evals=len(invalid_ind), **record)
    print((logbook.stream))

    # Begin the generational process
    for gen in range(1, options.genalg['ngenerations']):
        # Vary the population
        if True:
            #offspring1 = tools.selTournamentDCD(pop, len(pop)/2)
            offspring = selTournamentSum(pop, len(pop), 2)
            #offspring = offspring1 + offspring2
        else:
            offspring = tools.selTournamentDCD(pop, len(pop))
        print("type offspring:", type(offspring), "type offspring[0]:", type(offspring[0]))
        for ind in offspring:
            print("ind.fitness", ind.fitness)
        print("best_ind", best_ind, "type:", type(best_ind))
        print(best_ind.fitness)
        
        offspring = [toolbox.clone(ind) for ind in offspring]
        
        for ind1, ind2 in zip(offspring[::2], offspring[1::2]):
            if np.random.random() <= options.genalg['CXP']:
                toolbox.mate(ind1, ind2)
            
            toolbox.mutate(ind1)
            toolbox.mutate(ind2)
            #del ind1.fitness.values, ind2.fitness.values
            for ind in (ind1, ind2):
                if hasattr(ind, 'fitness'):
                    del ind.fitness.values

        # ELITISM IN GA
        #if not best[0] in offspring:
        #    offspring.append(best[0])
        
        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        fits = evaluate_multi(invalid_ind, gen=gen)
        for fit, ind in zip(fits, invalid_ind):
            if sum(best[1]) < sum(fit):
                best = (ind, fit)
            ind.fitness.values = fit

        # Select the next generation population
        #pop = toolbox.select(pop + offspring, MU)
        pop = toolbox.select(pop + offspring, options.genalg['npopulation'])
        best_ind = selBestSum(offspring, 1)[0]
        record = stats.compile(pop)
        logbook.record(gen=gen, evals=len(invalid_ind), **record)
        print((logbook.stream))
        print("best_ind:", contoind(best_ind), "fitness values:", best_ind.fitness)
        history.append({
            'gen':gen,
            'ncalcs':function.ncalcs,
            'index':contoind(best_ind),
            'p':sum(best_ind.fitness.values)
            })

    best_index = contoind(best_ind)
    result = {'p':sum(best_ind.fitness.values), 'index':best_index, 'history':history, 'ncalcs':history[-1]['ncalcs'], 'gen':options.genalg['ngenerations']}
    print("best:", best)

    return result


