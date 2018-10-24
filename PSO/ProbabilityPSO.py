import numpy as np
import random
import pprint
import sys
from copy import deepcopy

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

    def __init__(self, X, V=None, name=None):
        self.X = X
        self.V = V
        self.P = None
        self.name = name
        self.history = []
        self.nDim, self.nGroups = X.shape

        self.localbestX = deepcopy(X)
        self.localbestP = None
        self.localhistory = []

    def __repr__(self):
        return "<P{}:\n{}>".format(self.name, repr(self.X))

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

    def reset(self):
        self.sample = None
        self.P = None

    def update_local_best(self, array, epsilon):

        self.localbestP = self.P
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
                if abs(vel_cognitive)>1.0 or abs(vel_social)>1.0 or abs(vel_inertia)>1.0 or debug:
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

class CategoricalPSO(object):
    def __init__(self, array, n, function,
            parallel=False,
            maxiter=10,
            options=None):

        # initiate variables
        self.n = n
        self.array = array
        self.nDim, self.nGroups = np.array(array).shape
        self.function = function
        self.parallel = parallel
        self.maxiter = maxiter
        self.globalbestX = None
        self.globalbestP = None
        self.globalhistory = []

        # initiate swarm with positions and velocities
        self.init_swarm()

        # for now just set hyperparameters here:
        self.options = {'c1':1.49618, 'c2':1.49618, 'epsilon':0.75, 'w':0.729}
        if options:
            self.options.update(options)

        return

    def __repr__(self):
        return "<PSO object at {}>".format(id(self))

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

            particle = Particle(X, V, name=idx)

            swarm.append(particle)
        self.swarm = swarm
        return

    def update_global_best(self, particle):
        # set globalbest Property
        self.globalbestP = particle.P
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

    def evolve(self):
        iter=1
        while True:
            print "----------- Generation {} -----------".format(iter)

            # evaluate and set global best
            if self.parallel:
                pass
            else:
                for particle in self.swarm:
                    particle.evaluate(self.array, self.function)
                    if debug:
                        print "after evaluate for:", particle.name, "property:", particle.P, "from "\
                              "sample:",particle.sample

                    if particle.P <= particle.localbestP or particle.localbestP is None:
                        particle.update_local_best(self.array, self.options['epsilon'])

                    if particle.P <= self.globalbestP or self.globalbestP is None:
                        self.update_global_best(particle)

            # log values
            print "global bestP:", self.globalbestP
            self.globalhistory.append(self.globalbestP)
            for particle in self.swarm:
                particle.localhistory.append(particle.localbestP)

            # update velocities and position
            for particle in self.swarm:
                particle.updateV(self.globalbestX, self.options)
                particle.updateX()
                particle.reset()

            # decide if last iteration    
            iter+=1
            if iter>self.maxiter:
                break
        return

class Function4(object):
    def __init__(self):
        n = len(sitearray)
        self.k = 2
        self.NK = np.random.randn(nsites, *self.k*(n,))
        print "NK.shape:", self.NK.shape
        return

    def __call__(self, sample):
        p = 0.0
        for i in range(len(sample)):
            index=[i]
            for j in range(self.k):
                f = i+j
                if f > len(sample)-1:
                    f-= len(sample)
                g = list(array[f]).index(sample[f])
                index.append(g)
            try:
                p += self.NK[tuple(index)]
            except (IndexError, ValueError) as e:
                print "index:", index
        return p
function4 = Function4()

options = {'epsilon':0.65, 'c1':1.4, 'c2':1.4, 'w':0.8}
mypso = CategoricalPSO(array, 20, function4, maxiter=40, options=options)
print mypso
print mypso.swarm

mypso.evolve()

print "global bestX:", mypso.globalbestX
bestsample = []
for i in range(len(mypso.globalbestX)):
    igroup = list(mypso.globalbestX[i]).index(max(mypso.globalbestX[i]))
    group = array[i][igroup]
    bestsample.append(group)
print "best sample:", bestsample

print "X of best particle:", mypso[1]

print "global bestP:", mypso.globalbestP

totalhistory=[]
for particle in mypso.swarm:
    totalhistory.append(particle.history)

localhistory=[]
for particle in mypso.swarm:
    localhistory.append(particle.localhistory)

if debug:
    print "particle history:"
    pp.pprint(totalhistory)

    print "particle.localbest history:"
    pp.pprint(localhistory)

    print "global best history:"
    pp.pprint(mypso.globalhistory)

if sys.stdout.isatty():
    import matplotlib.pyplot as plt
    for i, a in enumerate(totalhistory):
        plt.plot(np.array(a)+i*0.02, alpha=0.9)
    plt.plot(mypso.globalhistory)
    plt.show()
    for i, a in enumerate(localhistory):
        plt.plot(np.array(a)+i*0.02, alpha=0.9)
    plt.show()
