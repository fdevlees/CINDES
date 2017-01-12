#!/bin/env python
debug=0
# python modules
import numpy as np
random = np.random.random
import pickle
from itertools import izip
import pprint

#needed by evolve
from time  import time
import logging
from sys   import platform as sys_platform
from sys   import stdout as sys_stdout
from CINDES4.pyevolve.GPopulation  import GPopulation
from CINDES4.pyevolve.GPopulation  import Util
import random as rrandom

# my own modules
#from writings import log_io, sprint, print_title
from CINDES4.utils.writings import log_io, sprint, print_title
from CINDES4 import INDES
from CINDES4.predictor import learning
from CINDES4.predictor import learning_int as ml_i

from CINDES4.pyevolve import G1DList , GSimpleGA, GAllele, Mutators, Initializators, Selectors, Consts, DBAdapters
from CINDES4.pyevolve import Scaling
import CINDES4.pyevolve as pyevolve

def skipper(conf):
    indje = INDES.procedures.zcon.contoind(conf)
    replaced = indje.replace('_','')
    if False:
        output= len(replaced)
    else:
        import string
        output = 0
        for i in replaced:
            try:
                output += string.uppercase.index(i)
            except ValueError:
                output += string.lowercase.index(i)
    return output

def my_mutator(conf):
    pos_to_mutate = np.random.randint(0, len(individual)-1)
    conf[pos_to_mutate] = np.random.choice( array[ pos_to_mutate ] )
    return conf

def make_individual(array, *args, **kwargs):
    arlen = len(array)
    conf = []
    for i in range(arlen):
        conf.append(np.random.choice(array[i]))
    return conf

def make_population(count,array):
    population = []
    for _ in xrange(count):
        individual = make_individual(array)
        population.append(individual)
    return population

def get_database():
    import pickle
    try:
        #with open('table_unbiased','rb') as f:
        with open('table_new3.dat','rb') as f:
            table = pickle.load(f)
    except IOError as e:
        print "NO TABLE ONLY VALID IF SKIPPER IS USED:,", e
        table = []
    return table

def get_input():
    options, subs_array = INDES.inputreader.read_input('INPUTBC')
    if debug: print "len(subs_array)", len(subs_array)
    return options, subs_array

@log_io()
def get_geometry(options):
    zmatrix = INDES.reader.geometry(zmatrixfile='ZMAT', **options)
    return zmatrix

class Fitness_Function():
    '''A class for the fitness functions. The class contains the attr's needed for evaluation. Here this will be
    the database that does not change. maybe even the kernel. see which part stays here and what part has to be 
    done in learning.py.
    Aim is that i can use:
        evaluator = FitnessFunction()
    and subsequently in each iteration
        evaluator.evaluate(population)
    '''
    def __init__(self, run, array=None, table=[]):
        '''for evaluation i need at least to have the database and the core / active / passive (all in zmatrix)
        i probably should also already get a self.kernel here such that the evaluatefunction only should call predict
        '''
        #self.zmatrix    = get_geometry(options)
        self.run = run
        self.table = table
        options = run.__dict__
        if options['ml']==1:  # depending on a not yet implemented option... 
            self.table      = get_database()
            from CINDES4.utils.converter import Converter
            self.converter = Converter()
            #kwargs['converter'] = self.converter
            self.initiate_machine_learning() #sets self.my_ML
        elif options['ml']==2:
            self.array      = array
            assert array!=None, "Give Array!"
            self.initiate_ml_int(**options)
        return

    def predict_via_submit_mono(self,conf):
        index = INDES.procedures.zcon.contoind(conf)
        print "index:", index
        confs = [conf]
        indices_tocal = [index]
        data_nocal = []
        myrun = self.run
        newy = INDES.procedures.submittingprocedure(confs,indices_tocal,data_nocal,myrun,**myrun.TZmat)
        print "newy:", newy[0]
        return newy[0][2]

    def indexmaker4(self,indices, confs):
        '''checks for confs already calculated'''
        indicesfull = indices[:]
        data = []
        if not self.table == []:
            for item in self.table:
                for index,confje in izip(indices[:],confs[:]):
                    if item[0] == index:
                        # remove that from the configurations
                        indices.remove(index)
                        confs.remove(confje)
                        # add that item from table to data
                        if item[1]==1:
                            raise SystemExit('elements in tablebin shouldnt be one')
                            data.append(item)
                        else:
                            new_item = item[:]
                            new_item.insert(1,1)
                            data.append(new_item)
            if not data == []:
                logging.info('filled data with ones already calced:' + pprint.pformat(data))
        return indices,data,confs #indicesfull are all the indices. 

    def log_table( self, data):
        if debug:
            print "in log_table: data:", data
            print "table:", self.table
        # here move the new data to table except duplicates
        for item in data:
            if item[1]==1:
                if not item[0] in [tja[0] for tja in self.table]:
                    tableitem = [item[0]] + item[2:]
                    self.table.append(tableitem)
                    if debug: print "tableitem:", tableitem
            else:
                assert item[1]==0, "item[1] has to be 1 or 0 but is %s" % str(item[1])
        with open('tablebin','wb') as tfid: # write the table to a file
            pickle.dump(self.table,tfid)
            print "dumped tablebin"
        return

    def predict_via_submit_multi(self,confs):
        indices = []
        for i in range(len(confs)):
            index = INDES.procedures.zcon.contoind(confs[i])
            indices.append(index)
        indices, data_nocal, confs = self.indexmaker4(indices, confs)
        myrun = self.run
        newy = INDES.procedures.submittingprocedure(confs,indices,data_nocal,myrun,**myrun.TZmat)
        if debug: print "newy:", newy
        self.log_table( newy)
        return newy

    @log_io()
    def evaluate_skip_multi(self,confs):
        indices = []
        for conf in confs:
            indices.append(INDES.procedures.zcon.contoind(conf))
        fitnesses = [ [index, skipper(index) ] for index in indices ]
        for conf,fitness in zip(confs,fitnesses):
            fitness[0] = conf
        sprint(10,fitnesses)
        return fitnesses

    def predict_via_precalculation(self, confs, new_y=[]):
        #print "confs:", confs
        index = INDES.procedures.zcon.contoind(confs)
        for item in new_y:
            if item[0] == index:
                y = item[2]
                if debug: print "y:", y
                return y
        else:
            pass
            #print "confs:", confs
            #raise SystemExit('stop for loop completed')

class My_GSimpleGA(GSimpleGA.GSimpleGA):

   def __init__(self,genome,run, precalculation=True, table=[]):
       GSimpleGA.GSimpleGA.__init__(self,genome)
       self.FF = Fitness_Function(run, table=table)
       self.precalculation = precalculation
       return

   def in_evolve(self, step=False, population=None):
      '''called in self.evolve and self.step to get the population and evaluate them'''
      populationlist = []
      if step:
          pop = population.internalPop
      else:
          pop = self.internalPop.internalPop
      for id in pop:
          populationlist.append( id.genomeList)

      print "populationlist", populationlist
      if not populationlist:
          raise SystemExit('stop')
      # calculate the population 
      new_y = self.FF.predict_via_submit_multi(populationlist)
      print "in in_evolve"
      return new_y

   def my_evaluate(self, step=False, population=None):
      populationlist = []
      if step:
          pop = population.internalPop
      else:
          pop = self.internalPop.internalPop
      for id in pop:
          populationlist.append( id.genomeList)
      print "populationlist", populationlist
      new_y = self.FF.predict_via_submit_multi(populationlist)
      y_dict = dict( [item[0], item[1:]] for item in new_y )
      for ind in population:
          index = INDES.procedures.zcon.contoind(ind.genomeList)
          print "individual:", ind.genomeList, "y:", y_dict[index], index
          ind.score = y_dict[index][1]
      #for ind, y in zip(population, new_y):
      #    print "ind:", ind.genomeList, "y:", y
      #    ind.score = y[2]
      #pop_sort = sorted(population)
      #for ind in population:
      #    ind.score = pop_sort.index(ind) + 100
      return

   def step(self):
      """ Just do one step in evolution, one generation """
      genomeMom = None
      genomeDad = None

      newPop = GPopulation(self.internalPop)
      logging.debug("Population was cloned.")

      size_iterate = len(self.internalPop)

      # Odd population size
      if size_iterate % 2 != 0: size_iterate -= 1

      crossover_empty = self.select(popID=self.currentGeneration).crossover.isEmpty()

      for i in xrange(0, size_iterate, 2):
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

      ############################################################## EVALUATE
      logging.info("Evaluating the new created population.")
      if self.precalculation:
          new_y = self.in_evolve(population=newPop, step=True)
          if debug: print "in step; new_y:", new_y
          newPop.evaluate(new_y=new_y)
      elif True:
          self.my_evaluate(step=True, population = newPop)
      else:
          newPop.evaluate()

      #Niching methods- Petrowski's clearing
      self.clear()

      if self.elitism:
         logging.debug("Doing elitism.")
         if self.getMinimax() == Consts.minimaxType["maximize"]:
            for i in xrange(self.nElitismReplacement):
               if self.internalPop.bestRaw(i).score > newPop.bestRaw(i).score:
                  newPop[len(newPop)-1-i] = self.internalPop.bestRaw(i)
         elif self.getMinimax() == Consts.minimaxType["minimize"]:
            for i in xrange(self.nElitismReplacement):
               if self.internalPop.bestRaw(i).score < newPop.bestRaw(i).score:
                  newPop[len(newPop)-1-i] = self.internalPop.bestRaw(i)

      self.internalPop = newPop
      self.internalPop.sort()

      logging.debug("The generation %d was finished.", self.currentGeneration)

      self.currentGeneration += 1

      return (self.currentGeneration == self.nGenerations)

   def evolve(self, freq_stats=0):
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
      if self.dbAdapter: self.dbAdapter.open(self)
      if self.migrationAdapter: self.migrationAdapter.start()


      if self.getGPMode():
         gp_function_prefix = self.getParam("gp_function_prefix")
         if gp_function_prefix is not None:
            self.__gp_catch_functions(gp_function_prefix)

      self.initialize()
      print "Jos in evolve"
      #print "self.internalPop:", self.internalPop
      #print "self.internalPop.internalPop[0]", self.internalPop.internalPop
      #print "self.internalPop.internalPop.genomeList", self.internalPop.internalPop[0].genomeList

      if self.precalculation:
          new_y = self.in_evolve()
          if debug: print "new_y:", new_y
          self.internalPop.evaluate(new_y=new_y)          ######### EVALUATE statement
      elif True:
          self.my_evaluate(population = self.internalPop)
      else:
          self.internalPop.evaluate()
      self.internalPop.sort()
      logging.debug("Starting loop over evolutionary algorithm.")

      try:
         while True:                                             ####### GenAlg loop

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
                  self.printStats()

            if self.dbAdapter:
               if self.currentGeneration % self.dbAdapter.getStatsGenFreq() == 0:
                  self.dumpStatsDB()

            if stopFlagTerminationCriteria:
               logging.debug("Evolution stopped by the Termination Criteria !")
               if freq_stats:
                  print "\n\tEvolution stopped by Termination Criteria function !\n"
               break

            if stopFlagCallback:
               logging.debug("Evolution stopped by Step Callback function !")
               if freq_stats:
                  print "\n\tEvolution stopped by Step Callback function !\n"
               break

            if self.interactiveMode:
               if sys_platform[:3] == "win":
                  if msvcrt.kbhit():
                     if ord(msvcrt.getch()) == Consts.CDefESCKey:
                        print "Loading modules for Interactive Mode...",
                        logging.debug("Windows Interactive Mode key detected ! generation=%d", self.getCurrentGeneration())
                        from pyevolve import Interaction
                        print " done !"
                        interact_banner = "## Pyevolve v.%s - Interactive Mode ##\nPress CTRL-Z to quit interactive mode." % (pyevolve.__version__,)
                        session_locals = { "ga_engine"  : self,
                                           "population" : self.getPopulation(),
                                           "pyevolve"   : pyevolve,
                                           "it"         : Interaction}
                        print
                        code.interact(interact_banner, local=session_locals)

               if (self.getInteractiveGeneration() >= 0) and (self.getInteractiveGeneration() == self.getCurrentGeneration()):
                        print "Loading modules for Interactive Mode...",
                        logging.debug("Manual Interactive Mode key detected ! generation=%d", self.getCurrentGeneration())
                        from pyevolve import Interaction
                        print " done !"
                        interact_banner = "## Pyevolve v.%s - Interactive Mode ##" % (pyevolve.__version__,)
                        session_locals = { "ga_engine"  : self,
                                           "population" : self.getPopulation(),
                                           "pyevolve"   : pyevolve,
                                           "it"         : Interaction}
                        print
                        code.interact(interact_banner, local=session_locals)

            b_max_iter = self.step()
            if b_max_iter: break #exit if the number of generations is equal to the max. number of gens.

      except KeyboardInterrupt:
         logging.debug("CTRL-C detected, finishing evolution.")
         if freq_stats: print "\n\tA break was detected, you have interrupted the evolution !\n"

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
         if freq_stats: print "Stopping the migration adapter... ",
         self.migrationAdapter.stop()
         if freq_stats: print "done !"

      return self.bestIndividual()

###### CALL(s) from __main__.py ###########

# 1. setup system
from CINDES4.INDES import procedures

def main(param, array):
    GArun = procedures.Run(**param)
    table = procedures.set_table(GArun)
    print "run object:\n", GArun
    get_genome(array, table, GArun)
    return

def get_genome(array,table, options):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    print "options:", options
    precalculation = False

    # Enable the logging system:
    pyevolve.logEnable()

    # Genome instance
    setOfAlleles = GAllele.GAlleles()
    for i in xrange(options.nsites):
       a = GAllele.GAlleleList(array[i])
       setOfAlleles.add(a)
    #for i in xrange(11, 20):
    #   # You can even add an object to the list
    #   a = GAllele.GAlleleList(['a','b', 'xxx', 666, 0])
    #   setOfAlleles.add(a)
    genome = G1DList.G1DList(options.nsites)
    genome.setParams(allele=setOfAlleles)

    # The evaluator function (objective function)
    if precalculation:
        FF = Fitness_Function(array=array,run = options)
        genome.evaluator.set(FF.predict_via_precalculation)
        #genome.evaluator.set(FF.int_predict_mono)
        #def predict_via_submit_mono(self,conf):
        #genome.evaluator.set(FF.predict_via_submit_mono)
    else:
        genome.evaluator.set(skipper)
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    genome.initializator.set(Initializators.G1DListInitializatorAllele)
    print "genome:\n", genome

    # Genetic Algorithm Instance
    ga = My_GSimpleGA(run = options, genome=genome, precalculation=precalculation, table=table)
    #ga.setMultiProcessing(True) #gives error thread.error: can't start new thread
    ga.selector.set(Selectors.GRouletteWheel)
    ga.setGenerations(100)
    #ga.setMinimax(Consts.minimaxType["minimize"])
    ga.setMutationRate(0.05) #i added this from another example
    # termination at convergence?:
    ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
    print "GenAlg:", ga

    # for negative fitness results:
    ga.setPopulationSize(5)
    pop = ga.getPopulation()
    pop.scaleMethod.set(Scaling.SigmaTruncScaling)

    # for plotting?
    sqlite_adapter = DBAdapters.DBSQLite(identify="ex4")
    ga.setDBAdapter(sqlite_adapter)

    # Do the evolution, with stats dump
    # frequency of 10 generations
    ga.evolve(freq_stats=50)
    return ga

if __name__ == "__main__":
    pass
