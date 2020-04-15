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


def eaSimple(population, 
             toolbox, cxpb, mutpb, ngen,
             stats=None,
             halloffame=None,
             nelitism=None,
             verbose=__debug__):
    """This algorithm reproduce the simplest evolutionary algorithm as
    presented in chapter 7 of [Back2000]_.
    :param population: A list of individuals.
    :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                    operators.
    :param cxpb: The probability of mating two individuals.
    :param mutpb: The probability of mutating an individual.
    :param ngen: The number of generation.
    :param stats: A :class:`~deap.tools.Statistics` object that is updated
                  inplace, optional.
    :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                       contain the best individuals, optional.
    :param verbose: Whether or not to log the statistics.
    :returns: The final population
    :returns: A class:`~deap.tools.Logbook` with the statistics of the
              evolution
    The algorithm takes in a population and evolves it in place using the
    :meth:`varAnd` method. It returns the optimized population and a
    :class:`~deap.tools.Logbook` with the statistics of the evolution. The
    logbook will contain the generation number, the number of evalutions for
    each generation and the statistics if a :class:`~deap.tools.Statistics` is
    given as argument. The *cxpb* and *mutpb* arguments are passed to the
    :func:`varAnd` function. The pseudocode goes as follow ::
        evaluate(population)
        for g in range(ngen):
            population = select(population, len(population))
            offspring = varAnd(population, toolbox, cxpb, mutpb)
            evaluate(offspring)
            population = offspring
    As stated in the pseudocode above, the algorithm goes as follow. First, it
    evaluates the individuals with an invalid fitness. Second, it enters the
    generational loop where the selection procedure is applied to entirely
    replace the parental population. The 1:1 replacement ratio of this
    algorithm **requires** the selection procedure to be stochastic and to
    select multiple times the same individual, for example,
    :func:`~deap.tools.selTournament` and :func:`~deap.tools.selRoulette`.
    Third, it applies the :func:`varAnd` function to produce the next
    generation population. Fourth, it evaluates the new individuals and
    compute the statistics on this population. Finally, when *ngen*
    generations are done, the algorithm returns a tuple with the final
    population and a :class:`~deap.tools.Logbook` of the evolution.
    .. note::
        Using a non-stochastic selection method will result in no selection as
        the operator selects *n* individuals from a pool of *n*.
    This function expects the :meth:`toolbox.mate`, :meth:`toolbox.mutate`,
    :meth:`toolbox.select` and :meth:`toolbox.evaluate` aliases to be
    registered in the toolbox.
    .. [Back2000] Back, Fogel and Michalewicz, "Evolutionary Computation 1 :
       Basic Algorithms and Operators", 2000.
    """
    history = []
    ncalcs = 0
    logbook = tools.Logbook()
    logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in population if not ind.fitness.valid]
    fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
    ncalcs += len(invalid_ind)
    for ind, fit in zip(invalid_ind, fitnesses):
        ind.fitness.values = fit

    if halloffame is not None:
        halloffame.update(population)

    record = stats.compile(population) if stats else {}
    logbook.record(gen=0, nevals=len(invalid_ind), **record)
    if verbose:
        print(logbook.stream)

    # Begin the generational process
    for gen in range(1, ngen + 1):
        # Select the next generation individuals
        offspring = toolbox.select(population, len(population))

        # Vary the pool of individuals
        offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
        ncalcs += len(invalid_ind)
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        # Update the hall of fame with the generated individuals
        if halloffame is not None:
            halloffame.update(offspring)

        # Replace the current population by the offspring
        if nelitism:
            population[:] = offspring + halloffame[:]
        else:
            population[:] = offspring

        # Append the current generation statistics to the logbook
        record = stats.compile(population) if stats else {}
        logbook.record(gen=gen, nevals=len(invalid_ind), **record)
        if verbose:
            print(logbook.stream)
        best_ind = tools.selBest(population, 1)[0]
        history.append({
            'gen':gen,
            'ncalcs':ncalcs,
            'index':contoind(best_ind),
            'p':best_ind.fitness.values[0]
            })

    return population, logbook, history

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

def selBestSum(individuals, k, fit_attr="fitness"):
    """Select the *k* best individuals among the input *individuals*. The
    list returned contains references to the input *individuals*.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :param fit_attr: The attribute of individuals to use as selection criterion
    :returns: A list containing the k best individuals.
    """
    #return sorted(individuals, key=attrgetter(fit_attr), reverse=True)[:k]
    print(individuals[0].fitness)
    print(repr(individuals[0].fitness))
    print(individuals[0].fitness.values)
    return sorted(individuals, key=lambda x:sum(getattr(x, fit_attr).values), reverse=True)[:k]
    #try:
    #    return sorted(individuals, key=lambda x:sum(getattr(x, fit_attr)), reverse=True)[:k]
    #except TypeError:
    #    return sorted(individuals, key=lambda x:sum(getattr(x, fit_attr).values), reverse=True)[:k]

def selRandom(individuals, k):
    """Select *k* individuals at random from the input *individuals* with
    replacement. The list returned contains references to the input
    *individuals*.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :returns: A list of selected individuals.
    This function uses the :func:`~random.choice` function from the
    python base :mod:`random` module.
    """
    return [rrandom.choice(individuals) for i in range(k)]


def selTournamentSum(individuals, k, tournsize, fit_attr="fitness"):
    """Select the best individual among *tournsize* randomly chosen
    individuals, *k* times. The list returned contains
    references to the input *individuals*.
    :param individuals: A list of individuals to select from.
    :param k: The number of individuals to select.
    :param tournsize: The number of individuals participating in each tournament.
    :param fit_attr: The attribute of individuals to use as selection criterion
    :returns: A list of selected individuals.
    This function uses the :func:`~random.choice` function from the python base
    :mod:`random` module.
    """
    chosen = []
    for i in range(k):
        aspirants = selRandom(individuals, tournsize)
        #chosen.append(max(aspirants, key=attrgetter(fit_attr)))
        chosen.append(max(aspirants, key=lambda x:sum(getattr(x, fit_attr).values)))
    return chosen


if __name__ == "__main__":
    pass
