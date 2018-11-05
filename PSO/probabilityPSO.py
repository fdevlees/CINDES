import numpy as np
import random
import pprint
import sys
from copy import deepcopy

from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.INDES.construction import indtocon, contoind

seed = 10
np.random.seed(seed)
random.seed(seed)
pp = pprint.PrettyPrinter(width=200)
np.set_printoptions(linewidth=120)

debug = False

def normalize(matrix):
    l1 = np.linalg.norm(matrix, ord=1, axis=1)
    if 0.0 in l1:
        raise ZeroDivisionError('bla')

    matrix = matrix / np.expand_dims(l1, axis=1)
    return matrix

class Particle(object):

    def __init__(self, X, V=None, index=None):
        self.X = X
        self.V = V
        self.P = None
        self.index = index
        self.history = []
        self.nDim, self.nGroups = X.shape

        self.localbestX = deepcopy(X)
        self.localbestP = None
        self.localbestsample = None
        self.localhistory = []

    def __repr__(self):
        return "<P{}:\n{}>".format(self.index, repr(self.X))

    def __getitem__(self, i):
        return self.X[i]

    def setsample(self, array):
        sample = []
        for sitegroup, siteprobs in zip(array, self.X):
            try:
                group = np.random.choice(sitegroup, p=siteprobs)
            except ValueError:
                print siteprobs
            sample.append(group)
        self.sample = sample
        return

    def evaluate(self, array, function):
        self.setsample(array)
        self.P = function(self.sample)
        self.history.append(self.P)
        return

    def setproperty(self, P):
        self.P = P
        self.history.append(self.P)
        return

    def reset(self):
        self.sample = None
        self.P = None

    def update_local_best(self, array, epsilon):

        self.localbestP = self.P
        self.localbestsample = self.sample
        #pp.pprint(self.localbestX)

        # for every site
        for i in range(self.nDim):
            n=0
            # get the index of the site that was performing best
            sampleindex = list(array[i]).index(self.sample[i])
            if debug: print "site:", i, "epsilon:", epsilon, "sampleindex", sampleindex
            totalcorrection = 0.0

            # change every group
            for j in range(self.nGroups):
                # if the group was not the good group
                if not j==sampleindex:
                    # decrease velocity of local best
                    new_prob = epsilon * self.localbestX[i][j]
                    totalcorrection += (1.0 - epsilon) * self.localbestX[i][j]
                    self.localbestX[i][j] = new_prob
                else:
                    n+=1

            self.localbestX[i][sampleindex] += totalcorrection
            if not n==1:
                 print "Assertion Error!"
                 print "n:", n
                 print "self.sample", self.sample
                 raise AssertionError()

        if debug:
            print "new self.local_best:", self.name
            pp.pprint(self.localbestX)
        self.localbestX = normalize(self.localbestX)
        return

    def updateV(self, globalbestX, options):
        if debug:
            print "old V:", self.name
            pp.pprint(self.V)
        #cognitive term:
        c1 = options['c1']
        #social term:
        c2 = options['c2']
        #inertia weight
        w  = options['w']

        if debug:
            print "X:"
            pp.pprint(self.X)
            print "localbestX:"
            pp.pprint(self.localbestX)

        for i in range(self.nDim):
            # here or within next loop?
            r1=np.random.random()
            r2=np.random.random()
            for j in range(self.nGroups):
                vel_cognitive = c1*r1*(self.localbestX[i][j]-self.X[i][j])
                vel_social    = c2*r2*(globalbestX[i][j]    -self.X[i][j])
                vel_inertia   =  w*self.V[i][j]
                if abs(vel_cognitive)>1.0 or abs(vel_social)>1.0 or abs(vel_inertia)>1.0 and debug:
                    print "inertia, cognitive, social, w, r1, r2, c1, c2"
                    print "vel corrections!:", vel_inertia, vel_cognitive, vel_social, w, r1, r2, c1, c2
                newv  = vel_inertia + vel_cognitive + vel_social
                #if newv>1.0:
                #    newv=1.0
                #elif newv<-1.0:
                #    newv=-1.0
                self.V[i][j]=newv

        if debug:
            print "updated V:"
            pp.pprint(self.V)

        # do we have to normalize V here?"
        #self.V = normalize(self.V)
        return

    def updateX(self):
        for i in range(self.nDim):
            for j in range(self.nGroups):
                newx=self.X[i][j]+self.V[i][j]

                if newx<0.0: newx=0.0
                if newx>1.0: newx=1.0

                self.X[i][j] = newx
        try:
            self.X = normalize(self.X)
        except ZeroDivisionError:
            print "V:"
            print self.V
            print "X:"
            print self.X

            raise
        return

class ProbabilityPSO(object):
    def __init__(self, array, npop, function,
            parallel=False,
            maxiter=10,
            options=None,
            minimize=True):

        # initiate variables
        self.iter = 0
        self.n = npop
        self.array = np.array(array)
        self.nDim, self.nGroups = np.array(array).shape
        self.function = function
        self.parallel = parallel
        self.maxiter = maxiter
        self.minimize = minimize
        self.globalbestX = None
        self.globalbestP = None
        self.globalbestsample = None
        self.globalhistory = []
        self.options = {'epsilon':0.75, 'w':0.8, 'c1':1.4, 'c2':1.4}
        if options:
            self.options.update(options)

        self.dbAdapter = None

        # initiate swarm with positions and velocities
        self.init_swarm()

        return

    def __repr__(self):
        return "<Probability-PSO object at {}>".format(id(self))

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
            newarray=deepcopy(self.array)
            for j in range(ndim):
                np.random.shuffle(newarray[j])
            shuffles.append(newarray)
        return np.concatenate(shuffles, axis=1).T[:self.n]

    def init_swarm(self):
        concrete_swarm = self.init_concrete_swarm()
        swarm = []
        for idx, Cparticle in enumerate(concrete_swarm):
            X = np.ones(self.array.shape)

            # change the ones to twos if position matches group
            for i, group in enumerate(Cparticle):
                index = list(self.array[i]).index(group)
                X[i][index]=3.0

            X = normalize(X)

            #get velocity either as random or as zero matrix
            if True:
                V = np.random.randn(*self.array.shape)
                V = normalize(V)
            else:
                V = np.zeros(self.array.shape)

            particle = Particle(X, V, index=idx)

            swarm.append(particle)
        self.swarm = swarm
        return

    def update_global_best(self, particle):
        # set globalbest Property
        self.globalbestP = particle.P
        self.globalbestsample = particle.sample
        # the first time a global best is set it gets the coordinates from the best particle
        # else only the probaility of the global best is changed 
        if self.globalbestX is None:
            self.globalbestX = particle.X
        else:
            for i in range(particle.nDim):
                sampleindex = list(self.array[i]).index(particle.sample[i])
                totalcorrection = 0.0
                for j in range(particle.nGroups):
                    if not j==sampleindex:
                        # decrease probability by multiplication with 0<epsilon<1
                        new_prob = self.options['epsilon'] * self.globalbestX[i][j]
                        totalcorrection += (1.0 - self.options['epsilon']) * self.globalbestX[i][j]
                        self.globalbestX[i][j] = new_prob
                self.globalbestX[i][sampleindex] += totalcorrection
            self.globalbestX = normalize(self.globalbestX)
        return

    @log_io()
    def evaluate(self):
        ''' this is my evaluate function '''
        # 1. get configurations to calculate:
        populationlist = []
        for particle in self.swarm:
            particle.setsample(self.array)
            populationlist.append(particle.sample)
        print "len populationlist:", len(populationlist),

        # 1.1. make confs hashable to make it a set and make it list again
        new_confs = tuple( '_'.join(item) for item in populationlist )
        unique_confs = [ indtocon(item) for item in set(new_confs)]
        print "n unique_confs:", len(unique_confs)
        print "new_confs:"
        for mol in set(new_confs):print mol

        # 2. call CINDES via FF to calculate the configurations
        mols = self.function.evaluate_multi(unique_confs, gen=self.iter)

        # 3. set the calculations to the correct indivual score
        y_dict = {mol.index: mol.Pvalue for mol in mols}
        print "y_dict:", y_dict
        for particle in self.swarm:
            index = "_".join(particle.sample)
            particle.setproperty(y_dict[index])
        return len(unique_confs)==1

    def evolve(self):
        self.iter=1
        nconvergence=0
        while True:
            print "----------- Generation {} -----------".format(self.iter)

            # evaluate and set global best
            if self.parallel:
                if self.evaluate():
                    nconvergence += 1
            else:
                for particle in self.swarm:
                    particle.evaluate(self.array, self.function)

            for particle in self.swarm:
                if debug:
                    print "after evaluate for:", particle.name, "property:", particle.P, "from "\
                          "sample:",particle.sample

                if self.minimize:
                    if particle.P <= particle.localbestP or particle.localbestP is None:
                        particle.update_local_best(self.array, self.options['epsilon'])
                    if particle.P <= self.globalbestP or self.globalbestP is None:
                        self.update_global_best(particle)
                else:
                    if particle.P >= particle.localbestP or particle.localbestP is None:
                        particle.update_local_best(self.array, self.options['epsilon'])
                    if particle.P >= self.globalbestP or self.globalbestP is None:
                        self.update_global_best(particle)

            # log values
            self.globalhistory.append(self.globalbestP)
            for particle in self.swarm:
                particle.localhistory.append(particle.localbestP)
            self.log()

            # update velocities and position
            for particle in self.swarm:
                particle.updateV(self.globalbestX, self.options)
                particle.updateX()
                particle.reset()

            print "global best Pvalue:", self.globalbestP
            print "global best sample:", self.globalbestsample

            # decide if last iteration    
            self.iter+=1
            if self.iter>self.maxiter:
                print "Run terminated. Max number of generations"
                break
            if nconvergence>1:
                print "Run Terminated due to onvergence criteria"
                break

        self.dbAdapter.commitAndClose()
        return

    def getStatistics(self):
        scores = [ particle.P for particle in self.swarm ]
        rawMin = np.min(scores)
        rawVar = np.var(scores)
        rawDev = np.std(scores)
        rawAve = np.average(scores)
        rawMax = np.max(scores)
        localscores = [ particle.localbestP for particle in self.swarm ]
        fitAve = np.average(localscores)
        fitMin = np.min(localscores)
        fitMax = np.max(localscores)
        return (rawMin, fitAve, fitMin, rawVar, rawDev, rawAve, fitMax, rawMax)

    def log(self):
        self.dbAdapter.insert(self)

    def setDBAdapter(self, dbAdapter):
        self.dbAdapter = dbAdapter
        self.dbAdapter.open(self)


