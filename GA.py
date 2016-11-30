#!/bin/env python
debug=1
# python modules
import numpy as np
random = np.random.random

# my own modules
from writings import log_io, sprint, print_title
import INDES
import learning
import learning_int as ml_i

class Genetic_Algorithm(object):
    def __init__(self, evaluator):
        self.f_cross = 0.0
        self.evaluator = evaluator # = fitness function

def skipper(conf):
    indje = INDES.zcon.contoind(conf)
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
    options, subs_array = INDES.read_input('INPUTBC')
    if debug: print "len(subs_array)", len(subs_array)
    return options, subs_array

@log_io()
def get_geometry(options):
    zmatrix = INDES.geometry('ZMAT',options,ilogging=False)
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
    def __init__(self, array=None, **options):
        '''for evaluation i need at least to have the database and the core / active / passive (all in zmatrix)
        i probably should also already get a self.kernel here such that the evaluatefunction only should call predict
        '''
        self.zmatrix    = get_geometry(options)
        self.table      = get_database()
        self.array      = array
        if options['ml']==1:  # depending on a not yet implemented option... 
            from converter import Converter
            self.converter = Converter()
            #kwargs['converter'] = self.converter
            self.initiate_machine_learning() #sets self.my_ML
        elif options['ml']==2:
            assert array!=None, "Give Array!"
            self.initiate_ml_int(**options)
        return

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
        from converter import Converter
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
        index = INDES.zcon.contoind(conf)
        print "index:", index
        newy  = self.my_ML.predict2([index],self.converter,**self.zmatrix)[0]
        print "newy:", newy
        return newy

    @log_io()
    def evaluate_skip_multi(self,confs):
        indices = []
        for conf in confs:
            indices.append(INDES.zcon.contoind(conf))
        fitnesses = [ [index, skipper(index) ] for index in indices ]
        for conf,fitness in zip(confs,fitnesses):
            fitness[0] = conf
        sprint(10,fitnesses)
        return fitnesses

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
    if True:
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
    ga.selector.set(Selectors.GRouletteWheel)
    ga.setGenerations(100)
    #ga.setMinimax(Consts.minimaxType["minimize"])
    ga.setMutationRate(0.05) #i added this from another example
    # termination at convergence?:
    ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
    print "GenAlg:", ga

    # for negative fitness results:
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
    print "index:", INDES.zcon.contoind(best)

if __name__ == "__main__":
    #procedure2()
    test_pyevolve3()
    #options , array = get_input()
    #print array
