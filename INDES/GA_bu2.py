#!/bin/env python
debug=1
# python modules
import numpy as np
random = np.random.random

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

#class Fitness_Function(learning.ML):
class Fitness_Function():
    '''A class for the fitness functions. The class contains the attr's needed for evaluation. Here this will be
    the database that does not change. maybe even the kernel. see which part stays here and what part has to be 
    done in learning.py.
    Aim is that i can use:
        evaluator = FitnessFunction()
    and subsequently in each iteration
        evaluator.evaluate(population)
    '''
    def __init__(self, run, array=None):
        '''for evaluation i need at least to have the database and the core / active / passive (all in zmatrix)
        i probably should also already get a self.kernel here such that the evaluatefunction only should call predict
        '''
        #self.zmatrix    = get_geometry(options)
        self.run = run
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

    def predict_via_submit_multi(self,confs):
        indices = []
        for i in range(len(confs)):
            index = INDES.procedures.zcon.contoind(confs[i])
            indices.append(index)
        data_nocal = []
        myrun = self.run
        newy = INDES.procedures.submittingprocedure(confs,indices,data_nocal,myrun,**myrun.TZmat)
        if debug: print "newy:", newy
        return newy

    def initiate_ml_int(self, sigma=1e2, labda=1e-7, **options):
        self.sigma = sigma
        self.labda = labda
        self.X = [ ml_i.indtoint(item[0], self.array) for item in self.table ]
        Y = np.asarray([ item[2] for item in self.table ])
        K = ml_i.get_kernel_int(self.X,sigma = sigma )
        self.alpha = ml_i.solver_int(K,Y, labda = labda)
        print "machine learning based on integer list initiated"
        return

    def int_predict_mono(self,conf):
        newx  = INDES.zcon.contoint(conf, self.array)
        newy  = ml_i.testnew_int(self.alpha, self.X, [newx], sigma=self.sigma)[0]
        return newy

    def initiate_machine_learning(self,printlevel=1,**kwargs):
        ''' or this function will be called by CINDES'''
        from CINDES4.utils.converter import Converter
        converter = Converter()
        kwargs['converter'] = converter

        self.my_ML = learning.MachineLearning('name',type='norm3')

        #now alpha and kernel are constructed
        self.my_ML.ML2(self.table,converter=self.converter,**self.zmatrix)

        if printlevel==1:
            print "alpha coefficients are calculated"
            sprint(5,self.my_ML.alpha)
            print "kernel[0:1]:"
            sprint(2,self.my_ML.kernel)
        return

    def predict_ml(self, population):
        indices = []
        for conf in population:
            indices.append(INDES.zcon.contoind(conf) )

        # now the indices has to be converted to xyz coordinates to. 
        new_y = self.my_ML.predict2(indices,self.converter,**self.zmatrix)
        fitnesses = zip(population, new_y)
        return fitnesses

    def predict_ml_mono(self, conf):
        index = INDES.procedures.zcon.contoind(conf)
        print "index:", index
        newy  = self.my_ML.predict2([index],self.converter,**self.zmatrix)[0]
        print "newy:", newy
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
        print "confs:", confs
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


def evaluate_skip_mono(conf):
    index = INDES.zcon.contoind(conf)
    print "index:", index
    fitness = INDES.skipper([index],iprint=False)
    print "fitness:", fitness
    fitness[0][0] = conf
    fitness = fitness[0]
    print "new_fitness:", fitness
    return fitness

def evolve(pop, array, function=skipper, target=0, retain=0.2, random_select=0.05, mutate=0.01):
    #print "pop:",
    #sprint(5,pop)
    #print "fitnessf( pop[1] ):", skipper(pop[2])
    #graded = []
    #for individual in pop:
    #    print "individual:", individual
    #    fitness = skipper(individual)
    #    print "fitness individual:", fitness
    #    graded.append( [ individual, fitness ] )
    if False:
        graded = [ [ individual , function(individual) ] for individual in pop ]
    else:
        graded = function(pop)
    graded = [ x for x in sorted(graded,key = lambda x:x[1])]
    graded, values = zip(*graded)
    graded = list(graded)
    print "values:", values
    average = np.mean(values)
    std     = np.std(values)
    print " average :" , average
    print "     std :" , std

    # keep 20% of best performing #lowest (sorted = lowest first)
    retain_length = int(len(graded)*retain) #define what amount is retain fraction
    parents = graded[:retain_length] #take that amount

    # randomly add other individuals to change = 5%
    # promote genetic diversity
    for individual in graded[retain_length:]:
        if random_select > random():
            parents.append(individual)
    # mutate some individuals
    for individual in parents:
        if mutate > random(): # = 1% here
            pos_to_mutate = np.random.randint(0, len(individual)-1)
            # this mutation is not ideal, because it
            # restricts the range of possible values,
            # but the function is unaware of the min/max
            # values used to create the individuals,
            #individual[pos_to_mutate] = np.random.randint( min(individual), max(individual))
            individual[pos_to_mutate] = np.random.choice( array[ pos_to_mutate ] )
            #print "mutation performed"

    # crossover parents to create children
    parents_length = len(parents)
    desired_length = len(pop) - parents_length
    children = []
    while len(children) < desired_length:
        male = np.random.randint(0, parents_length-1)
        female = np.random.randint(0, parents_length-1)
        if male != female:
            male = parents[male]
            female = parents[female]
            half = len(male) / 2
            child = male[:half] + female[half:]
            children.append(child)
    parents.extend(children)
    return parents, average, std

def procedure1():
    if True:
        p_count         = 50
        n_generations   = 20
        options, array  = get_input()
        population      = make_population(p_count, array)
        if debug:
            print "first 10 of population:"
            sprint(10,population)
        FF              = Fitness_Function(options)
        #FF.initiate_machine_learning() #this is done automatically
        #fitnesses = FF.predict_ml(population)
        #print "fitnesses:"
        #print(2,fitnesses)
        averages = []
        stds     = []
        for i in xrange(n_generations):
            print_title("ITERATION: "+str(i),outline='l',signator='^')
            new_population, average, std = evolve(population, array, function=FF.predict_ml)
            #new_population, average, std = evolve(population, array, function=FF.evaluate_skip_multi, mutate=0.05)
            print "average:", average
            averages.append(average)
            stds.append(std)
            population = new_population
            if debug:
                print "population:"
                sprint(10,population)

        print averages
        #sprint(10, population)
        import matplotlib.pyplot as plt
        import seaborn
        plt.errorbar(range(len(averages)),averages,stds)
        plt.show()

def procedure2():
    if True:
        p_count = 20
        n_generations = 1000
        options, array = get_input()
        population     = make_population(p_count,array)
        #zmat           = get_geometry(options)
        print "population:", 
        sprint(10,population)
        print "*"*10
        FF              = Fitness_Function(options)
        averages = []
        stds=[]
        for _ in xrange(n_generations):
            #new_population, average = evolve(population,array,mutate=0.02)
            new_population, average, std = evolve(population, array, function=FF.evaluate_skip_multi, mutate=0.05)
            print "average:", average
            averages.append(average)
            stds.append(std)
            population = new_population

        print averages
        sprint(10, population)
        import matplotlib.pyplot as plt
        import seaborn
        plt.errorbar(range(len(averages)),averages,stds)
        plt.show()


class My_GSimpleGA(GSimpleGA.GSimpleGA):

   def __init__(self,genome,run):
       GSimpleGA.GSimpleGA.__init__(self,genome)
       self.FF = Fitness_Function(run)

   def in_evolve(self, step=False, population=None):
      '''called in self.evolve and self.step to get the population and evaluate them'''
      populationlist = []
      if step:
          pop = population.internalPop
          #pop = population
          #print "pop:", pop
          #raise SystemExit('stop')
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
      new_y = self.in_evolve(population=newPop, step=True)
      

      if debug: print "in step; new_y:", new_y
      newPop.evaluate(new_y=new_y)

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
      print "self.internalPop:", self.internalPop
      print "self.internalPop.internalPop[0]", self.internalPop.internalPop
      print "self.internalPop.internalPop.genomeList", self.internalPop.internalPop[0].genomeList

      new_y = self.in_evolve()
      if debug: print "new_y:", new_y
      self.internalPop.evaluate(new_y=new_y)          ######### EVALUATE statement
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

def test_pyevolve():
    # This function is the evaluation function, we want
    # to give high score to more zero'ed chromosomes
    def eval_func(chromosome):
        score = 0.0
        # iterate over the chromosome elements (items)
        for value in chromosome:
            if value==0:
                score += 1.0
        return score
    from pyevolve import G1DList , GSimpleGA
    genome = G1DList.G1DList(20)
    print "genome:", genome
    genome.evaluator.set(eval_func)
    ga = GSimpleGA.GSimpleGA(genome)
    print "ga:", ga
    ga.evolve(freq_stats=10)
    print ga.bestIndividual()

def test_pyevolve2():
    from pyevolve import G1DList , GSimpleGA
    genome = G1DList.G1DList(size=10)
    genome.evaluator.set(skipper)
    genome.mutator.set(my_mutator)
    genome.initializator.set(make_individual)
    print "genome:", genome
    ga = GSimpleGA.GSimpleGA(genome)
    print "ga:", ga
    ga.evolve(freq_stats=10)
    print ga.bestIndividual()

def test_pyevolve3(*args,**kwargs):
    '''going to try allele'''
    from pyevolve import G1DList , GSimpleGA, GAllele, Mutators, Initializators, Selectors, Consts, DBAdapters
    from pyevolve import Scaling
    import pyevolve
    # Enable the logging system:
    pyevolve.logEnable()

    # get input from INPUTBC inputfile
    options, array = get_input()

    # Genome instance
    setOfAlleles = GAllele.GAlleles()
    for i in xrange(10):
       a = GAllele.GAlleleList(array[i])
       setOfAlleles.add(a)
    #for i in xrange(11, 20):
    #   # You can even add an object to the list
    #   a = GAllele.GAlleleList(['a','b', 'xxx', 666, 0])
    #   setOfAlleles.add(a)
    genome = G1DList.G1DList(10)
    genome.setParams(allele=setOfAlleles)

    # The evaluator function (objective function)
    if False:
        FF = Fitness_Function(array,**options)
        #genome.evaluator.set(FF.int_predict_mono)
        genome.evaluator.set(FF.predict_ml_mono)
    else:
        genome.evaluator.set(skipper)
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    genome.initializator.set(Initializators.G1DListInitializatorAllele)
    print "genome:\n", genome

    # Genetic Algorithm Instance
    ga = GSimpleGA.GSimpleGA(genome)
    ga.setMultiProcessing(True)
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
    sqlite_adapter = DBAdapters.DBSQLite(identify="ex3")
    ga.setDBAdapter(sqlite_adapter)

    # Do the evolution, with stats dump
    # frequency of 10 generations
    ga.evolve(freq_stats=50)

    # Best individual
    best =  ga.bestIndividual()
    print "\n Best individual score: %.2f" % best.score
    print best
    print "index:", INDES.procedures.zcon.contoind(best)

###### CALL(s) from __main__.py ###########

# 1. setup system
from CINDES4.INDES import procedures

def main(param, array):
    GArun = procedures.Run(**param)
    table = procedures.set_table(GArun)
    print "run object:\n", GArun
    get_genome(array, GArun)
    pass

def get_genome(array,options):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    print "options:", options

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
    if True:
        FF = Fitness_Function(array=array,run = options)
        #genome.evaluator.set(FF.int_predict_mono)
        #def predict_via_submit_mono(self,conf):
        #genome.evaluator.set(FF.predict_via_submit_mono)
        genome.evaluator.set(FF.predict_via_precalculation)
    else:
        genome.evaluator.set(skipper)
    genome.mutator.set(Mutators.G1DListMutatorAllele)
    genome.initializator.set(Initializators.G1DListInitializatorAllele)
    print "genome:\n", genome

    # Genetic Algorithm Instance
    ga = My_GSimpleGA(run = options, genome=genome)
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
    #procedure2()
    test_pyevolve3()
    #options , array = get_input()
    #print array
