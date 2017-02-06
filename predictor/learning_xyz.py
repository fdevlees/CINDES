#!/bin/env python
'''module for machine learning for CINDES2.py program'''
debug=1
BoB = False

#import pyximport; pyximport.install()
#import cython_sum

import logging
import sys
import pickle
import math as m
import gc
import numpy as np
import matplotlib.pyplot as plt

#for a progres bar:
import time
import progressbar
bar = progressbar.ProgressBar()

from sklearn.model_selection import train_test_split
from sklearn import preprocessing


from copy import deepcopy
from pprint import pprint

from CINDES4.utils.writings import log_io, sprint, print_title
from CINDES4.utils.converter import Converter
from CINDES4.INDES import construction as zcon
#from converter import Converter
#import construction as zcon

class MachineLearning(object):
    ''' Class for making training set / kernel / coulomb / predictions etc. '''
    def __init__(self,
                 converter=Converter(), 
                 descriptor='normal',
                 kerneltype='laplacian',
                 sigma = 1e8, 
                 labda = 1e-5, 
                 cutoff = [],
                 doscaling = False,
                 validation_fraction = 0.0,
                 inputfile = 'table.xyz',
                 **kwargs):
        self.type = descriptor
        self.kerneltype = kerneltype
        self.converter  = converter
        self.sigma = sigma
        self.labda = labda
        self.cutoff= cutoff
        self.inputfile = inputfile
        self.validation_fraction = validation_fraction
        self.doscaling = doscaling
        print_title("As a kernel: "+kerneltype+" is used",outline='l',signator='k',newlines=True)
        self.get_input(self.converter) # sets self.X / self.y
        if self.doscaling:
            self.scaler = preprocessing.StandardScaler().fit(self.X)
            X_scaled = scaler.transform(self.X)
            self.X = X_scaled
        if validation_fraction:
            self.set_validation_set()
        return

    def __repr__(self):
        ret = "\nMachine Learning class using the following parameters:\n"
        ret += "  inputfile          : " + str(self.inputfile) + '\n'
        ret += "  kernel             : " + str(self.kerneltype) + '\n'
        ret += "  descriptor         : " + str(self.type)       + '\n'
        ret += "  sigma              : " + str(self.sigma)      + '\n'
        ret += "  labda              : " + str(self.labda)      + '\n'
        ret += "  use scaling        : " + str(bool(self.doscaling)) + '\n'
        ret += "  use validation set : " + str(bool(self.validation_fraction))
        ret += ': ' + str(self.validation_fraction) + '\n'
        ret += "  use crossvalidation: " + str(False) + '\n'
        ret += "and the following data characteristics:\n"
        assert len(self.X)==len(self.y)
        ret += "  ndata     : " + str(len(self.X))     + '\n'
        ret += "-----"
        return ret

    def __str__(self):
        return self.__repr__()

    def get_input(self,converter, data=1):
        y = []
        xyzs = []
        with open(self.inputfile) as fid:
            while True: #for all xyzs
                headerline = fid.readline()
                if not headerline: break
                natoms = int(headerline.split()[1])
                if data:
                    line = fid.readline().split()
                    y.append( float(line[1]) )   # 0: gap, 1: homo 2: lumo 3: Etotal
                xyz = []
                for i in xrange(natoms): #for each atom
                    line = fid.readline().split()
                    #print 'line:', line
                    assert not line=='\n'
                    xyztje = np.zeros([3])
                    #print "xyztje", xyztje
                    for j in xrange(3): #for x,y,z
                        xyztje[j] = line[j+1]
                    xyz.append([line[0],xyztje,converter.masses[line[0]]])
                xyzs.append(xyz)
                fid.readline() #empty line

        if self.cutoff and data:
            valcutoffmin, valcutoffplus = args.cutoff
            print "cutoff applied of ", str(valcutoffmin), "and", str(valcutoffplus), " eV"
            nbefore = len(y)
            xyzs, y =  zip ( *[ item for item in zip(xyzs,y) if not ( item[1]<valcutoffmin  or item[1]>valcutoffplus )] )
            nafter = len(y)
            print str( nbefore - nafter ) ,"elements were removed from list"
            #data = [ [ item[0], float(item[1])] for item in datar if float(item[1])<cutoff ] 

        if data:  #total set. with data
            self.xyzs = xyzs
            self.y = np.asarray(y)
        else:
            self.xyzs_t = xyzs

        # get coulombs of training set
        if self.type=='BoB':
            self.X = np.asarray( tuple( self.BoB(item) for item in self.xyzs) )
        else:
            self.X = np.asarray( tuple( self.coulomb(item) for item in self.xyzs) )
        return

    def coulomb(self, xyz):
        #print "xyz:", xyz
        l = len(xyz)
        C = np.zeros([l,l])
        for i in xrange(l):
            for j in range(l): #changed this from i+1 to j
                if i==j:
                    C[i][i]= 0.5*xyz[i][2]**(2.4)
                else:
                    t = xyz[i][2] * xyz[j][2]
                    n = np.sqrt(
                            np.sum(
                                np.square(
                                    xyz[i][1] - xyz[j][1] ) ) )
                    C[i][j] = t/n
                    # let not make it symmetric. because we don't use these elements                  
        print "&",
        if self.type=='norm1': #return a sorted Coulomb matrix based on norm
            return symsort(C)
        elif self.type=='norm2':
            return np.sort(np.linalg.norm(C,axis=0))[::-1]
        elif self.type=='norm3':
            #make a sorted Coulomb matrix
            C = symsort(C)
            if True:
                maxn = 56   # max no of atoms for the adamantane and diamantane derivatives. 
                C = padzeros(C,maxn=maxn)
                return C[ np.tril_indices(maxn) ]
            else:
                # return a lower triangular matrix of the coulomb matrix. 
                return C[ np.tril_indices(l) ]
        elif self.type=='norm4':
            # here no triangularization
            C = symsort(C)
            C = padzeros(C, maxn=56)
            #plotmat(C,log=False)
            return C.flatten()
        else:
            assert self.type=='normal'
            return C

    def BoB(self, xyz):
        d = False
        from collections import OrderedDict
        l = len(xyz)
        Bag = []
        # max no of atoms per type
        typef = OrderedDict((( 'H' , 36 ),
                             ( 'C' , 20 ),
                             ( 'O' , 20 ),
                             ( 'N' , 10 ),
                             ( 'F' , 30 ),
                             ( 'S' , 10 ),
                             ( 'Cl', 10 ) ) )
        types = typef.keys()
        # make a list of typef with max no of combination of atom1 with atom2 
        def trianglen(typef,key1,key2):
            if key1==key2: return int( .5 * typef[key1] * ( typef[key1] - 1 ) )
            else:          return typef[key1]*typef[key2]
        ncombis = OrderedDict( ( (''.join(sorted((key1,key2),key = lambda x: typef.keys().index(x)) ), trianglen(typef,key1, key2) ) for key1 in
            types for key2 in types ) )
        ntypes= len(types)
        monos = OrderedDict(( (key,[]) for key in types ))
        # the next line makes dicts of every possible atom combination with combined keys. combined in order as in typef!
        duos  = OrderedDict( ( (''.join(sorted((key1,key2),key = lambda x: types.index(x)) ), [] ) for key1 in types for key2 in types ) )
        if d:
            print "monos, duos:", monos, duos
            print "typef, ncombis", typef, ncombis
        for i in xrange(l): # for every atom
            for j in range(i,l): #so for every combination with that atom not yet visited
                if i==j:
                    nuclear = 0.5*xyz[i][2]**(2.4)
                    monos[xyz[i][0]].append(nuclear)
                else:
                    t = xyz[i][2] * xyz[j][2]
                    n = np.sqrt(
                            np.sum(
                                np.square(
                                    xyz[i][1] - xyz[j][1] ) ) )
                    force = t/n
                    duo_indices = ( xyz[i][0], xyz[j][0] )
                    duo_key = ''.join(sorted(duo_indices ,key=lambda x: types.index(x) ))
                    duos[duo_key].append(force)
        # now sorted every item in the dictionaries and pad with zeros
        if d:
            print "monos, duos:", monos, duos
        for dictio, ntypes in ( (monos, typef  ),
                                (duos , ncombis)):
            for key,value in dictio.iteritems():
                N = ntypes[key]
                dictio[key] = sorted(dictio[key])[::-1] + [0.0] * ( N - len(value) )
        if d:
            print "monos, duos:", monos, duos
        # merge everything together orderly
        monos_flat = [ x for v in monos.itervalues() for x in v ]
        duos_flat  = [ x for v in duos.itervalues()  for x in v ]
        Bag = monos_flat + duos_flat
        print "B",
        if d:
            print "Bag:", Bag
            print "len(Bag):", len(Bag)
            raise SystemExit('stop')
        return Bag

    def set_estimator(self, estimator=None):
        if estimator:
            self.clf = estimator
        else:
            self.clf = KernelRidge( kernel = self.kerneltype, alpha = self.labda )
        return

    def set_validation_set(self):
        X = self.X[:]
        y = self.y[:]
        self.X, self.X_val, self.y, self.y_val = train_test_split( X, y, test_size= self.validation_fraction, random_state=3)
        return

    def GridSearch(self,param_grid=None, set_new_clf = True, **kwargs):
        from sklearn.model_selection import GridSearchCV
        from sklearn.kernel_ridge import KernelRidge
        if not param_grid:
            param_grid = {'alpha' : [1e4, 1e2, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5,1e-6 ],
                          'gamma' : np.logspace(-6, 6, 10) }
        krr = GridSearchCV( self.clf ,
                            cv    = 5,
                            param_grid = param_grid,
                            iid   = False, # no assumption of universal distribution of data 
                            refit = True,  # refit the best estimator with the entire dataset
                            return_train_score = True )
        print "krr:", krr
        krr_out = krr.fit(self.X, self.y)
        print "krr_out:", krr_out
        #print "krr_out results:", krr_out.cv_results_
        print "best parameters:", krr.best_params_
        print "best score:", krr.best_score_

        if set_new_clf:
            self.clf = krr_out
        return krr.best_params_

    def predict(self,X=None):
        if X==None:
            y_pred = self.clf.predict(self.X_val)
        else:
            y_pred = self.clf.predict( X )
        return y_pred


def KDE(data1,data2):
    scores = np.array([ i - j for i,j in zip(data1,data2) ])
    if True:
        import matplotlib.pyplot as plt
        from sklearn.neighbors import KernelDensity
        fig, ax = plt.subplots()
        if True:
            from sklearn.grid_search import GridSearchCV
            grid = GridSearchCV( KernelDensity(),
                                 { 'bandwidth': np.logspace(-3, 1,20) } , cv=8)
            grid.fit(scores[:,np.newaxis])
            print "gridsearch cv for best bandwidth for KDE:", grid.best_params_
            xlim = grid.best_params_['bandwidth']* 20
            kde = grid.best_estimator_
        else:
            kde = KernelDensity(kernel='gaussian', bandwidth=args.bandwidth).fit(scores[:,np.newaxis])
            xlim = 20*args.bandwidth
        X_plot = np.linspace(-xlim, xlim, 1000)[:, np.newaxis]
        log_dens = kde.score_samples(X_plot)
        print "logdens:"
        sprint(10,log_dens)
        ax.plot(X_plot, np.exp(log_dens), 'r-')
        #plt.show()
    if True:
        import matplotlib.pyplot as plt
        import matplotlib.mlab as mlab
        from scipy.stats import norm
        (f_mu, f_std) = norm.fit(scores)
        nbars = 10
        n, bins, patches = plt.hist(scores, nbars, normed=1, facecolor='green', alpha=0.75)
        mu = np.mean(scores)
        std = np.std(scores)
        print "mu:", mu
        print "std:", std
        y = mlab.normpdf( bins,mu,std)
        l = ax.plot(bins, y, 'r--', linewidth=1)
        f_y = mlab.normpdf( bins,f_mu,f_std)
        f_l = ax.plot(bins, f_y, 'b--', linewidth=1)
        #plt.show()
    return log_dens

def get_sigma( Cs ):
    ''' trying to get sigma defined from the maximum distance
        see paper Lilienfled and Ramakrishnan
        sigma_opt = Dmax / Log(2) '''
    l1 = len(Cs)
    maxd = None
    for i in range(l1):
        for j in range(i):
            if i==j: continue
            d = distance( Cs[i], Cs[j] )
            if d > maxd:
                maxd = d
    print "maxd:", maxd
    return maxd / np.log(2)

def distance(C1, C2, kerneltype='gaussian'):
    '''calculate euclidian distance between two coulomb matrices'''
    assert C1.shape==C2.shape, 'coulomb matrices not the same size'
    # make a difference array D
    D = C1.copy()
    D -= C2
    if kerneltype == "gaussian":
        # sum of quadratic differences
        #total = np.sqrt( np.sum( D**2 ) )
        total = np.sum( D**2 )
        # now the funny thing is euclidean norm is np.sqrt(sum(elements)) but this norm has to be squared. and so sqrt(x)**2 = x 
    elif kerneltype == "laplacian":
        total = np.sum( np.abs(D) )
    return total

def symsort(mat):
    indexlist = np.argsort(np.linalg.norm(mat,axis=1))[::-1]
    return mat[indexlist][:,indexlist]

def padzeros(M, maxn=76 ):
    ''' pads the matrix M until size is size*size '''
    npad = maxn - M.shape[0]
    padM = np.pad(M, (0, npad), 'constant', constant_values=0.0)
    return padM

def MAE(data1,data2):
    assert len(data1)==len(data2)
    npoints = len(data1)
    mae = 0
    for i in xrange(npoints):
        mae += abs( data1[i] - data2[i] )
    return mae / float(npoints)

def RMSE(data1,data2):
    assert len(data1)==len(data2)
    npoints = len(data1)
    rmse = 0
    for i in xrange(npoints):
        rmse += ( data1[i] - data2[i] )**2
    return np.sqrt( rmse / float(npoints) )

def pearsonr(data1,data2):
    import scipy.stats as ss
    R = ss.pearsonr(data1,data2)
    return R[0]

def plot_error(data1,data2, *args, **kwargs):
    import matplotlib.pyplot as plt
    plt.plot(data1,data2, *args, **kwargs)
    return

def plotmat(mat, log=1):
    if log:
        mat = np.log(mat)
    import matplotlib.pyplot as plt
    plt.matshow(mat)
    #ax = plt.gca()
    plt.colorbar()
    plt.show()
    return

################ UTILS #########################

def zmatoxyz(a,mat):
    #if debug: print "a:", a
    zmat = a.read_zmalist(mat)
    return a.zmatrix_to_cartesian()

def contozma(conf,core,active,passive,**kwargs):
    c = deepcopy(core)
    a = deepcopy(active)
    p = deepcopy(passive)
    mat = zcon.constructor2(conf,c,a,p)
    return mat

############################ FOR GENERATE TABLE.XYZ PROCEDURE ########################

def generate_xyz(indices,converter,outputfile='table.xyz',y=0,*args,**kwargs):
    # get a list of configurations:
    confs = [ zcon.indtocon(item) for item in indices ]
    printindices = 1
    try:
       if type(y[0]) == str: 
           ty = 1
       elif type(y[0]) in (tuple,list):
           ty = 2
    except TypeError:
       ty=0
       pass
    # textfile open
    with open(outputfile,'w') as fid:

        for i in range(len(confs)):
            print i, "indices[i]", indices[i]
            logging.debug("i=" + str(i))
            mat = contozma(confs[i],**kwargs)
            xyz = zmatoxyz(converter,mat)
            if printindices==1:
                fid.write('{:04d} {:4d} {:s}\n'.format(i+1,len(xyz), indices[i]))
            else:
                fid.write('{:03d} {:4d}\n'.format(i+1,len(xyz)))
            if ty==1:
                fid.write('{:12.8f}\n'.format(y[i]))
            elif ty==2:
                for item in y[i]:
                    fid.write( ' {:12.8f} '.format(item) )
                fid.write('\n')
            for item in xyz:
                fid.write('{:3s} {:10.4f} {:10.4f} {:10.4f}\n'.format(item[0],item[1][0],item[1][1],item[1][2]))
            fid.write('\n')
            # now for each mat transform to xyz
    return

def generate1(converter,table=[],**kwargs):
    if table==[]:
        with open('tablebin','rb') as f:
            table = pickle.load(f)

    # get the property vector Y
    #Y = [ item[1] for item in table ]
    Y = [ item[1:] for item in table ]

    # for each index in tablebin get zmat
    # get a list of indices
    indices = [ item[0] for item in table ]

    #Y, indices = zip( * [ item for item in zip(Y,indices) if not 'turned' in H    

    generate_xyz(indices=indices,y=Y,converter=converter,**kwargs)
    return

################## GRID PARAMETER SEARCH ####################

def search_grid( my_ML, sigmas, labdas ):
    total_scores = []
    for sigma in sigmas:
        for labda in labdas:
            score, std = my_ML.CV(sigma=sigma, labda=labda)
            total_scores.append( ( sigma, labda, score, std ) )
    return total_scores

