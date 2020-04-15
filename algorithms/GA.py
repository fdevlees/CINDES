#!/bin/env python
debug = False
verbose = False
# python modules
import numpy as np
random = np.random.random
#import pickle
import json

# needed by evolve
from time import time
import logging
from sys import platform as sys_platform
from sys import stdout as sys_stdout
import random as rrandom

# my own modules
#from writings import log_io, sprint, print_title
from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.molecule import Molecule
from CINDES.utils.table import set_table
from CINDES.utils.utils import skipper
from CINDES.evaluation.calculator import evaluate_mols
from CINDES.evaluation.construction import contoind, indtocon

from CINDES.algorithms import loggings
from CINDES.algorithms.pyevolve.GPopulation import GPopulation
from CINDES.algorithms.pyevolve.GPopulation import Util
from CINDES.algorithms.pyevolve import G1DList, GSimpleGA, GAllele, Mutators, Initializators, Selectors, Consts, DBAdapters, Crossovers
from CINDES.algorithms.pyevolve import Scaling
import CINDES.algorithms.pyevolve as pyevolve


class Fitness_Function():
    '''A class for the fitness functions. The class contains the attr's needed for evaluation. Here this will be
    the database that does not change. maybe even the kernel. see which part stays here and what part has to be
    done in learning.py.
    Aim is that i can use:
        evaluator = FitnessFunction()
    and subsequently in each iteration
        evaluator.evaluate(population)
    '''

    def __init__(self, run, array=None, table=dict()):
        '''for evaluation i need at least to have the database and the core / active / passive (all in zmatrix)
        i probably should also already get a self.kernel here such that the evaluatefunction only should call predict
        '''
        self.run = run
        self.table = table
        self.array = array
        self.ncalcs = 0
        return

    def evaluate_multi(self, confs, gen=0):
        ''' this function is used by my_GSimpleGA class.my_evaluate '''
        if self.run.write: print("confs:", confs)

        individuals = [Molecule(conf=conf) for conf in confs]  # list of molecules

        if self.run.write: print("indivuals:", individuals)

        mols, nnewcalcs, made_pred = evaluate_mols(self.run, individuals, self.table, gen, nsite=0)
        self.ncalcs += nnewcalcs

        self.table = loggings.loggings(mols,
                                             self.table,
                                             gen,
                                             1, 1,
                                             made_pred=made_pred,
                                             tablename=self.run.tablename,
                                             write=self.run.write)

        return mols


class My_GSimpleGA(GSimpleGA.GSimpleGA):

    def __init__(self, genome, function=None, precalculation=True, restart=False, **kwargs):
        GSimpleGA.GSimpleGA.__init__(self, genome, **kwargs)

        self.restart = restart

        self.precalculation = precalculation
        if self.precalculation:
            self.FF = function
        return


    @log_io()
    def my_evaluate(self, step=False, population=None):
        ''' this is my evaluate function '''
        # 1. get configurations to calculate:
        populationlist = []
        if step:
            pop = population.internalPop
        else:
            pop = self.internalPop.internalPop
        for id in pop:
            populationlist.append(id.genomeList)
        #print "populationlist", populationlist, "len:", len(populationlist)

        # 1.1. make confs hashable to make it a set and make it list again
        new_confs = tuple(tuple(map(tuple, item)) for item in populationlist)
        unique_confs = [list(map(list, item)) for item in set(new_confs)]
        logging.info("n unique_confs: {:d}".format(len(unique_confs)))
        if self.track:
            try:
                self.ncalcs += len(unique_confs)
            except AttributeError:
                self.ncalcs = len(unique_confs)
        #print "unique_confs:", unique_confs

        # 2. call CINDES via FF to calculate the configurations
        logging.info("BestIndividual (score): {} ({})".format(self.bestIndividual().genomeList, self.bestIndividual().score))
        mols = self.FF.evaluate_multi(unique_confs, gen=self.currentGeneration)

        # 3. set the calculations to the correct indivual score
        y_dict = {mol.index: mol.Pvalue for mol in mols}
        for ind in population:
            index = contoind(ind.genomeList)
            if verbose: print("individual:", ind.genomeList, "y:", y_dict[index], index)
            ind.score = y_dict[index]
            ind.index = index
        return

    def step(self):
        """ Just do one step in evolution, one generation """
        genomeMom = None
        genomeDad = None

        newPop = GPopulation(self.internalPop)
        logging.debug("Population was cloned.")

        size_iterate = len(self.internalPop)

        # Odd population size
        if size_iterate % 2 != 0:
            size_iterate -= 1

        crossover_empty = self.select(popID=self.currentGeneration).crossover.isEmpty()

        for i in range(0, size_iterate, 2):
            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)

            if not crossover_empty and self.pCrossover >= 1.0:
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2):
                    (sister, brother) = it
            else:
                if not crossover_empty and Util.randomFlipCoin(self.pCrossover):
                    for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=2):
                        (sister, brother) = it
                else:
                    sister = genomeMom.clone()
                    brother = genomeDad.clone()

            sister.mutate(pmut=self.pMutation, ga_engine=self)
            brother.mutate(pmut=self.pMutation, ga_engine=self)

            newPop.internalPop.append(sister)
            newPop.internalPop.append(brother)

        if len(self.internalPop) % 2 != 0:
            genomeMom = self.select(popID=self.currentGeneration)
            genomeDad = self.select(popID=self.currentGeneration)

            if Util.randomFlipCoin(self.pCrossover):
                for it in genomeMom.crossover.applyFunctions(mom=genomeMom, dad=genomeDad, count=1):
                    (sister, brother) = it
            else:
                sister = rrandom.choice([genomeMom, genomeDad])
                sister = sister.clone()
                sister.mutate(pmut=self.pMutation, ga_engine=self)

            newPop.internalPop.append(sister)

        # EVALUATE
        logging.info("Evaluating the new created population.")
        if self.precalculation:
            self.my_evaluate(step=True, population=newPop)
        else:
            newPop.evaluate()

        # Niching methods- Petrowski's clearing
        self.clear()

        if self.elitism:
            logging.debug("Doing elitism.")
            if self.getMinimax() == Consts.minimaxType["maximize"]:
                for i in range(self.nElitismReplacement):
                    if self.internalPop.bestRaw(i).score > newPop.bestRaw(i).score:
                        newPop[len(newPop) - 1 - i] = self.internalPop.bestRaw(i)
            elif self.getMinimax() == Consts.minimaxType["minimize"]:
                for i in range(self.nElitismReplacement):
                    if self.internalPop.bestRaw(i).score < newPop.bestRaw(i).score:
                        newPop[len(newPop) - 1 - i] = self.internalPop.bestRaw(i)

        self.internalPop = newPop
        self.internalPop.sort()

        logging.debug("The generation %d was finished.", self.currentGeneration)

        self.currentGeneration += 1
        print_title('Generation: %d' % self.currentGeneration, outline='l', signator='=')

        return (self.currentGeneration == self.nGenerations)

    def evolve(self, freq_stats=10):
        """ Do all the generations until the termination criteria, accepts
        the freq_stats (default is 0) to dump statistics at n-generation

        Example:
           >>> ga_engine.evolve(freq_stats=10)
           (...)

        :param freq_stats: if greater than 0, the statistics will be
                           printed every freq_stats generation.
        :rtype: returns the best individual of the evolution

        .. versionadded:: 0.6
           the return of the best individual

        """

        stopFlagCallback = False
        stopFlagTerminationCriteria = False

        self.time_init = time()

        logging.debug("Starting the DB Adapter and the Migration Adapter if any")
        if self.dbAdapter:
            self.dbAdapter.open(self)
        if self.migrationAdapter:
            self.migrationAdapter.start()

        if self.getGPMode():
            gp_function_prefix = self.getParam("gp_function_prefix")
            if gp_function_prefix is not None:
                self.__gp_catch_functions(gp_function_prefix)

        self.initialize()

        #print "Jos in evolve"
        #print "self.internalPop:", self.internalPop
        #print "self.internalPop.internalPop", self.internalPop.internalPop[0]
        print("self.internalPop.internalPop.genomeList", self.internalPop.internalPop[1].genomeList)
        print("n start Individuals:", len(self.internalPop))

        if self.restart:
            restartPop = self.getLastPopulation()
            for startInd, restartInd in zip(self.internalPop, restartPop):
                print("i:", startInd.genomeList)
                startInd.genomeList = restartInd
            print("n restart Individuals:", len(restartPop))
            print("self.internalPop:", self.internalPop)
            #raise NotImplementedError('bla')

        if self.precalculation:
            self.my_evaluate(population=self.internalPop)
        else:
            self.internalPop.evaluate()
        self.internalPop.sort()
        logging.debug("Starting loop over evolutionary algorithm.")

        try:
            while True:  # GenAlg loop

                if self.migrationAdapter:
                    logging.debug("Migration adapter: exchange")
                    self.migrationAdapter.exchange()
                    self.internalPop.clearFlags()
                    self.internalPop.sort()

                if not self.stepCallback.isEmpty():
                    for it in self.stepCallback.applyFunctions(self):
                        stopFlagCallback = it

                if not self.terminationCriteria.isEmpty():
                    for it in self.terminationCriteria.applyFunctions(self):
                        stopFlagTerminationCriteria = it

                if freq_stats:
                    if (self.currentGeneration % freq_stats == 0) or (self.getCurrentGeneration() == 0):
                        #logging.info("freq_stats:",
                        self.printStats()

                if self.dbAdapter:
                    if self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0:
                        self.dumpStatsDB()

                if stopFlagTerminationCriteria:
                    logging.debug("Evolution stopped by the Termination Criteria !")
                    if freq_stats:
                        print("\n\tEvolution stopped by Termination Criteria function !\n")
                    break

                if stopFlagCallback:
                    logging.debug("Evolution stopped by Step Callback function !")
                    if freq_stats:
                        print("\n\tEvolution stopped by Step Callback function !\n")
                    break

                if self.interactiveMode:
                    if sys_platform[:3] == "win":
                        if msvcrt.kbhit():
                            if ord(msvcrt.getch()) == Consts.CDefESCKey:
                                print("Loading modules for Interactive Mode...", end=' ')
                                logging.debug("Windows Interactive Mode key detected ! generation=%d",
                                              self.getCurrentGeneration())
                                from CINDES.algorithms.pyevolve import Interaction
                                print(" done !")
                                interact_banner = "## Pyevolve v.%s - Interactive Mode ##\nPress CTRL-Z to quit interactive mode." % (
                                    pyevolve.__version__,)
                                session_locals = {"ga_engine": self,
                                                  "population": self.getPopulation(),
                                                  "pyevolve": pyevolve,
                                                  "it": Interaction}
                                print()
                                code.interact(interact_banner, local=session_locals)

                    if (self.getInteractiveGeneration() >= 0) and (
                            self.getInteractiveGeneration() == self.getCurrentGeneration()):
                        print("Loading modules for Interactive Mode...", end=' ')
                        logging.debug("Manual Interactive Mode key detected ! generation=%d",
                                      self.getCurrentGeneration())
                        from CINDES.algorithms.pyevolve import Interaction
                        print(" done !")
                        interact_banner = "## Pyevolve v.%s - Interactive Mode ##" % (pyevolve.__version__,)
                        session_locals = {"ga_engine": self,
                                          "population": self.getPopulation(),
                                          "pyevolve": pyevolve,
                                          "it": Interaction}
                        code.interact(interact_banner, local=session_locals)

                b_max_iter = self.step()
                if b_max_iter:
                    break  # exit if the number of generations is equal to the max. number of gens.

        except KeyboardInterrupt:
            logging.debug("CTRL-C detected, finishing evolution.")
            if freq_stats:
                print("\n\tA break was detected, you have interrupted the evolution !\n")

        if freq_stats != 0:
            self.printStats()
            self.printTimeElapsed()

        if self.dbAdapter:
            logging.debug("Closing the DB Adapter")
            if not (self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0):
                self.dumpStatsDB()
            self.dbAdapter.commitAndClose()

        if self.migrationAdapter:
            logging.debug("Closing the Migration Adapter")
            if freq_stats:
                print("Stopping the migration adapter... ", end=' ')
            self.migrationAdapter.stop()
            if freq_stats:
                print("done !")

        return self.bestIndividual()

    def getLastPopulation(self):
        generation, rawPop = self.dbAdapter.getLastPopulation()
        self.currentGeneration = generation
        pop = [ indtocon(ind) for ind in rawPop ]
        print("new Pop", pop)
        return pop

###### CALL(s) from __main__.py ###########



def main(GArun):
    table = set_table(GArun, GArun.array)
    function = Fitness_Function(GArun, table=table, array=GArun.array)
    final_genome = run_pyevolve(GArun.array, GArun, function=function)
    best = final_genome.bestIndividual()
    print("final_genome:", final_genome)
    print("best:", best)
    print("best.~ score fitness genomelist :", best.score, best.fitness, best.genomeList)
    return final_genome


def run_pyevolve(array, options, level=None, function=None):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''

    # 0.
    logging.info("options:"+ repr(options))

    # 1. Enable the logging system:
    if level is None:
        pyevolve.logEnable()
    else:
        pyevolve.logEnable(level=level)

    # 2. Set Genome instance using as allelles the sites with the different functionalisations.
    setOfAlleles = GAllele.GAlleles()
    nalleles = options.nsites
    for i in range(nalleles):
        a = GAllele.GAlleleList(array[i])
        setOfAlleles.add(a)
    genome = G1DList.G1DList(nalleles)
    genome.setParams(allele=setOfAlleles)

    # 3. Set evaluator function (objective function) or set precalculation is True! this circumvents serial evaluation
    precalculation = True
    if not precalculation:
        genome.evaluator.set(skipper)
        # is made that moment.

    # 4. Set mutator function
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    # one could also add another mutator for exapmple to swap sites
    # genome.mutator.add(Mutators.G1DListMutatorSwap)

    # 5. Set initalizator function
    genome.initializator.set(Initializators.G1DListInitializatorAllele)

    # 6. set Crossover function: types: G1DListCrossoverUniform, G1DListCrossoverSinglePoint, G1DListCrossoverTwoPoint
    if not options.genalg['CXP'] == 0.0:
        genome.crossover.set(Crossovers.G1DListCrossoverUniform)
    logging.info("genome:\n" + repr(genome))

    # 7.1 set Genetic Algorithm Instance using a defined random.seed()
    if not 'seed' in options.genalg:
        options.genalg['seed'] = np.random.randint(1, 9999)
    logging.info("seed to generate randomness: {:d}".format(options.genalg['seed']))

    # 7.2 set trackinghistory flag
    if not 'trackhistory' in options.genalg:
        options.genalg['trackhistory']=False

    # 7.3 set Genetic Algorithm Instance
    restart = options.genalg['restart']
    ga = My_GSimpleGA(
            function=function,
            genome=genome,
            precalculation=precalculation,
            seed=options.genalg['seed'],
            trackhistory=options.genalg['trackhistory'],
            restart=restart )
    pop = ga.getPopulation()

    # 8. set Selector
    #if options.genalg['selector'] == 'RouletteWheel':  # Default = GRouletteWheel
    if any(item in options.genalg['selector'] for item in ['Roulette', 'roulette','wheel']):
        ga.selector.set(Selectors.GRouletteWheel)
    elif any(item in options.genalg['selector'] for item in ['Rank', 'rank']):
        ga.selector.set(Selectors.GRankSelector)
    elif any(item in options.genalg['selector'] for item in ['Uni', 'uni']):
        ga.selector.set(Selectors.GUniformSelector)
    elif any(item in options.genalg['selector'] for item in ['Tour', 'tour']):
        ga.selector.set(Selectors.GTournamentSelector)
        pop.setParams(tournamentPool=options.genalg['tournamentPoolsize'])
    else:
        raise SystemExit('no valid selector is chosen')

    # 9. set NGEN (number of generations)
    ga.setGenerations(options.genalg['ngenerations'])

    # 10. set min / max (optimize to a maximum or to a minimum)
    if options.optimum in ['min', 'minimize', 'minimum']:
        ga.setMinimax(Consts.minimaxType["minimize"])

    # 11. set MUP (mutation probability)
    ga.setMutationRate(options.genalg['MUP'])  # i added this from another example

    # 12. set CXP (crossover probability)
    if not options.genalg['CXP'] == 0.0:
        ga.setCrossoverRate(options.genalg['CXP'])

    # 13. set termination at convergence?:
    ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)

    # 14. set population size
    ga.setPopulationSize(options.genalg['npopulation'])

    # 15. set elitism
    if options.genalg['nelitism']:
        ga.setElitism(bool(options.genalg['nelitism']))
        logging.info("n elitism: {:d}".format(options.genalg['nelitism']))
        ga.nElitismReplacement = options.genalg['nelitism']

    # 16. set scaling: to allow for negative scores we have to use
    # SigmaTruncScaling. otherwise also LinearScaling or PowerLawScaling could
    # be used
    pop.scaleMethod.set(Scaling.SigmaTruncScaling)

    # 17. for plotting / logging
    if options.write or True:
        if restart:
            resetIdentify=False
        else:
            resetIdentify=True
        sqlite_adapter = DBAdapters.DBSQLite(
            identify=options.genalg['db_identify'], resetDB=False, resetIdentify=resetIdentify, commit_freq=1)
        ga.setDBAdapter(sqlite_adapter)

    logging.info("GenAlg:"+ repr(ga))

    # Do the evolution, with stats dump
    ga.evolve(freq_stats=options.genalg['freq_stats'])
    return ga


if __name__ == "__main__":
    pass
