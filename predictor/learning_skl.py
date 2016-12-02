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

from copy import deepcopy
from pprint import pprint

from CINDES4.utils.writings import log_io, sprint, print_title
from CINDES4.utils.converter import Converter
from CINDES4.INDES import construction as zcon
#from converter import Converter
#import construction as zcon

class MachineLearning(object):
    ''' Class for making training set / kernel / coulomb / predictions etc. '''
    def __init__(self,converter, type='normal',kerneltype='laplacian',table=[],sigma = 1e8, labda = 1e-5, **kwargs):
        self.type = type
        self.kerneltype = kerneltype
        self.converter  = converter
        self.sigma = sigma
        self.labda = labda
        print_title("As a kernel: "+kerneltype+" is used",outline='l',signator='k',newlines=True)
        if table==[]:
            self.get_input(self.converter)
        else:
            self.setXY(table, **kwargs)
        return

    def get_input(self,converter, data=1, inputfile='table.xyz'):
        y = []
        xyzs = []
        with open(inputfile) as fid:
            while True: #for all xyzs
                headerline = fid.readline()
                if not headerline: break
                #print "headerline:", headerline
                natoms = int(headerline.split()[1])
                #print "natoms:", natoms
                if data:
                    #y.append( float(fid.readline()) )
                    line = fid.readline().split()
                    #print line
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
        if args.cutoff and data:
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
        if BoB:
            self.coulombs = np.asarray( tuple( self.BoB(item) for item in self.xyzs) )
        else:
            self.coulombs = np.asarray( tuple( self.coulomb(item) for item in self.xyzs) )
        return

    def setXY(self,table, core,active,passive):
        ''' sets self.coulombs which is self.X and self.y '''
        tableindex = 2
        print "sigma: ", str(self.sigma), ' ', 'labda: ', str(self.labda)

        #### MAKE Y
        # get input from inputfile table
        self.y =  np.fromiter((item[tableindex] for item in table ),np.float)

        #### MAKE X
        ## X.1: convert indices to confs and convert confs to ZMATRICES
        mats =  tuple( contozma(zcon.indtocon(item[0]),core,active,passive) for item in table )
        ## X.2: convert ZMATRICES to XYZ coordinates
        try:
            if debug: print "self.converter:", self.converter
            self.xyzs = [ zmatoxyz(self.converter,item) for item in mats ]
            if False:
                from CINDES4.utils.molecule import Molecule
                #from molecule import Molecule
                for i in range(len(self.xyzs)):
                    mol = Molecule()
                    mol.set_xyz(self.xyzs[i])
                    mol.set_OBMol()
                    self.xyzs[i] = mol.optimize()
                    #print self.xyzs[i]
                    print "!",

        except KeyError:
            print "Error: with:", item
            i = mats.index(item)
            print "index:", table[i]
            raise
        ## X.3: convert XYZs to COULOMB MATRICES
        # calculate all coulomb matrices
        #self.coulombs = np.asarray( tuple( self.coulomb(item) for item in self.xyzs) )
        if BoB:
            self.coulombs = np.asarray( tuple( self.BoB(item) for item in self.xyzs) )
        else:
            self.coulombs = np.asarray( tuple( self.coulomb(item) for item in self.xyzs) )

        # analyze coulombs:
        #if debug:
        #    print "ncoulombs:", len(self.coulombs)
        #    print "coulombs[3]:"
        #    sprint(100,self.coulombs[3])
        #    raise SystemExit('stop')

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


    @log_io()
    #@profile
    def get_kernel(self, sigma=1e4):
        self.sigma = sigma
        l = len(self.coulombs)
        K = np.zeros([l,l])
        #for i in bar(range(l)):
        for i in range(l):
            #for j in range(l):
            #    if not i==j:
            #        K[i][j] = np.exp( self.distance(self.coulombs[i], self.coulombs[j]) / ( 2 * sigma**2) )
            #    else:
            #        K[i][j] = 1.0
            for j in range(i+1): #for i in range(0) gives []
                if not i==j:
                    K[i][j] = np.exp( self.distance(self.coulombs[i], self.coulombs[j]) / ( 2 * sigma**2) )
                    K[j][i] = K[i][j]
                else:
                    K[i][j] = 1.0
            print "#",
        self.kernel = K
        return

    @log_io()
    def solver(self,kernel,ny, labda=None):
        if not labda:
            labda = self.labda
        l = len(kernel)
        I = np.identity(l)
        Ka = kernel + labda * I
        #print "ny:", ny
        alg = 7
        if alg == 1:
            alpha = np.dot( np.linalg.inv(Ka) , ny.T)
        elif alg == 2:
            from scipy import linalg
            alpha = np.dot( linalg.inv(Ka) , ny.T )
        elif alg == 3:
            L = np.linalg.cholesky(Ka)
        elif alg == 4:
            print "Ka shape:", np.shape(Ka)
            print "b  shape:", np.shape(ny)
            alpha = np.linalg.solve(Ka,ny) #gives error: matrix not positive-definite # that is not all eigvals are positive
        elif alg == 5:
            alpha = np.linalg.lstsq(Ka,ny)[0]
        elif alg == 6:
            from scipy import linalg
            cho = linalg.cho_factor(Ka)
            alpha = linalg.cho_solve(cho,ny)
            #print 'residual', linalg.norm(np.dot(Ka, alpha) - ny)/ linalg.norm(Ka)
            # gives: numpy.linalg.linalg.LinAlgError: 2-th leading minor not positive definite
        elif alg == 7:
            from scipy import linalg
            lu = linalg.lu_factor(Ka)
            alpha = linalg.lu_solve(lu, ny)
        print "alphashape:", alpha.shape
        self.alpha = alpha
        return alpha

    def testnew(self,newX):
        outtest = []
        for coulombje in self.coulombs:
        #for j in bar(  xrange( len(self.coulombs_t) )  ):
            ans = 0
            for i in xrange(len(self.coulombs)):
                #ans += self.alpha[i] * np.exp( self.distance( self.coulombs_t[j], self.coulombs[i] ) / ( 2 * self.sigma **2) )
                ans += self.alpha[i] * np.exp( - self.distance( coulombje, self.coulombs[i] ) / ( 2 * self.sigma **2) )
            outtest.append(ans)
            print "$",
        return outtest

    def testnew_no_skl(self,kernel_new, alpha=None):
        if alpha==None:
            alpha = self.alpha
        # #
        if debug:
            print "shape   kernel_new:", kernel_new.shape
            print "(self?)alpha.shape:", alpha.shape
        newY = np.dot( kernel_new , alpha )
        if debug:
            print "newY.shape:", newY.shape
        return newY

    def predict2(self, indices,**kwargs):
        # 1. convert new indices to confs to ZMAT
        mat_inds = tuple( contozma(zcon.indtocon(item),**kwargs) for item in indices)
        # 2. convert new ZMATs to XYZs
        self.xyzs_t = [ zmatoxyz(self.converter, item) for item in mat_inds ]
        # 3. convert XYZs to COULOMBs
        if BoB:
            self.coulombs_t = tuple( self.BoB(item) for item in self.xyzs_t )
        else:
            self.coulombs_t = tuple( self.coulomb(item) for item in self.xyzs_t )
        return self.coulombs_t


    def learn(self):
        # Set SKLEARN 
        from sklearn.kernel_ridge import KernelRidge
        if True:   # USE SKLEARN - IMPLEMENTED KERNEL
            sigma = get_sigma( self.coulombs )
            #gamma = 1. / ( 2 * sigma**2 )
            gamma = 1. / ( sigma )
            print "sigma:", sigma
            clf = KernelRidge( alpha = self.labda, kernel='rbf' , gamma = gamma)
        elif True:      # USE own KERNEL as callable. NOT WORKING!
            if True:
                #determine sigma
                # 1. get all distances. 
                sigma = get_sigma( self.coulombs )
                if self.kerneltype == 'gaussian': sigma = np.sqrt(sigma)
                print "sigma used:", sigma
            else:
                sigma = self.sigma
            if self.kerneltype == 'gaussian':
                factor = 1. / ( 2 * sigma**2 )
            elif self.kerneltype == 'laplacian':
                factor = 1. / sigma
            clf = KernelRidge( alpha = self.labda, kernel=my_kernel3 , kernel_params = { 'factor': factor , 'kerneltype': self.kerneltype } )
        print "clf:", clf
        clf_out = clf.fit( self.coulombs, self.y )
        print "clf_out:", clf_out
        self.clf = clf
        return clf

    def predict(self,new_indices,**kwargs):
        newX = self.predict2( new_indices, **kwargs)
        newy = self.clf.predict( newX )
        return newy

    def learn1(self):
        '''learn with using precomputed matrices'''
        def validate_score(X,y,sigma,labda):
            from sklearn.model_selection import KFold
            kf = KFold(n_splits=3)
            scores = []
            for train, test in kf.split(self.coulombs): # so for 
                # split
                X_train, X_test, y_train, y_test = X[train], X[test], y[train], y[test]
                if debug:
                    print "X_train.shape:", X_train.shape
                    print "X_test.shape:", X_test.shape
                    print "Y_train.shape:", y_train.shape
                    print "Y_test.shape:", y_test.shape
                # make kernel
                clf = KernelRidge( alpha = labda, kernel= 'precomputed' )
                kernel = my_kernel( X_train, X_train , sigma=sigma )
                if debug:
                    print "kernel.shape:", kernel.shape
                # train
                fit_out = clf.fit( kernel, y_train)
                # test new
                test_kernel = my_kernel( X_test, X_train, sigma=sigma)
                score = clf.score( test_kernel, y_test)
                scores.append(score)
            print "average score", np.average(scores)
            return np.average(scores)
        # Set SKLEARN 
        from sklearn.kernel_ridge import KernelRidge
        clf = KernelRidge( alpha = self.labda, kernel= 'precomputed' )
        self.clf = clf
        if True:
            # calculate kernel
            if True:
                sigma = get_sigma( self.coulombs )
                if self.kerneltype=='gaussian': sigma = np.sqrt(sigma)
                self.sigma = sigma
            print "use sigma:", sigma
            kernel = my_kernel2( self.coulombs, self.coulombs, sigma=self.sigma , kerneltype = self.kerneltype)
            if debug:
                print "first 2 of kernel:"
                sprint(2, kernel)
            # calculate alpha coefficients
            fit_out = clf.fit( kernel , self.y )
            print "fit_out:", fit_out
            # evaluate 
            print "train_score:", clf.score( kernel, self.y )
            print "params:", clf.get_params()
            print "coef0", clf.coef0
        else:
            # cross validation
            X = np.asarray(self.coulombs)
            y = self.y
            sigmas = [ 10**i for i in np.arange(-5, 10,2) ]
            labdas = [ 10**i for i in range(-10, 5,2) ]
            sup_scores = np.zeros([10,10])
            for i in range(len(sigmas)):
                for j in range(len(labdas)):
                    sup_scores[i,j] = validate_score(X,y,sigmas[i],labdas[j])
                    print "labda %s sigma: %s score: %s" %(labdas[j], sigmas[i], sup_scores[i,j] )
            print "sub_scores:"
            print sup_scores
            import matplotlib.pyplot as plt
            plt.matshow(sup_scores)
            plt.show()

    def predict1(self, new_indices,**kwargs):
        # 1 set_coulomb_t:
        newX = self.predict2( new_indices, **kwargs)
        pred_kernel = my_kernel2( newX, self.coulombs, sigma=self.sigma, kerneltype = self.kerneltype)
        #pred_kernel = my_kernel2( self.coulombs, newX , sigma=self.sigma)
        test_out = self.clf.predict( pred_kernel )
        print "test_out:", test_out
        #raise SystemExit('stop')
        return test_out

    def learn_no_skl(self):
        ''' do ML without the sklearn package '''

        # 1. calculate ideal sigma:
        if False:
            sigma = get_sigma( self.coulombs)
            self.sigma = sigma

        # 2. make kernel
        kernel = my_kernel2( self.coulombs, self.coulombs, sigma=self.sigma, kerneltype = self.kerneltype)

        # #
        if debug:
            print "first 2 of kernel:"
            sprint(2, kernel)

        # 3. learn: make alphas. sets at self.alpha
        alpha = self.solver(kernel, self.y)
        return

    def predict_no_skl(self, new_indices, alpha=None, **kwargs):
        newX = self.predict2( new_indices, **kwargs)
        pred_kernel = my_kernel2( newX, self.coulombs, sigma= self.sigma, kerneltype = self.kerneltype)

        test_out = self.testnew_no_skl(pred_kernel)

        return test_out

    def CV(self, labda, sigma, n_splits=4, skl=False):
        def CV_calcs():
            kernel_train= my_kernel2( X_train, X_train, sigma=sigma, kerneltype = self.kerneltype )
            if skl:
                clf = KernelRidge( alpha = labda, kernel= 'precomputed' )
                fit_out = clf.fit( kernel_train , self.y )
            else:
                alpha       = self.solver( kernel_train, y_train , labda = labda)
            kernel_test = my_kernel2( X_test , X_train, sigma=sigma, kerneltype = self.kerneltype )
            if skl:
                y_test_pred = clf.predict( kernel_test)
            else:
                y_test_pred = self.testnew_no_skl( kernel_test , alpha = alpha )
            return y_test_pred

        from sklearn.model_selection import KFold
        kf = KFold(n_splits=n_splits)
        scores = []
        for train,test in kf.split(self.coulombs, self.y):
            X_train = self.coulombs[train]
            y_train = self.y[train]
            X_test  = self.coulombs[test]
            y_test  = self.y[test]
            if False:
                y_test_pred = CV_calcs()
                if debug:
                    print "test y:"
                    sprint(10, y_test)
                    print "test y pred:"
                    sprint(10, y_test_pred)
            else:
                from sklearn.kernel_ridge import KernelRidge
                gamma = 1. / (2* sigma**2 )
                print "sigma:", sigma
                clf = KernelRidge( alpha = labda, kernel='rbf' , gamma = gamma)
                clf_out = clf.fit( X_train, y_train )
                print "clf_out:", clf_out
                y_test_pred = clf.predict( X_test )
            score = RMSE( y_test_pred, y_test )
            scores.append(score)
            print "sigma, labda, score:", sigma, labda, score
        avg_score = np.average(scores)
        std_score = np.std(scores)
        return avg_score, std_score

    def GridSearch(self, indices=[],**kwargs):
        from sklearn.model_selection import GridSearchCV
        from sklearn.kernel_ridge import KernelRidge
        krr = GridSearchCV( KernelRidge( kernel= 'rbf', gamma=0.1 ),
                            cv    = 5,
                            param_grid = {'alpha' : [1e4, 1e2, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5,1e-6 ],
                                          'gamma' : np.logspace(-6, 6, 10) } ) 
        krr_out = krr.fit(self.coulombs, self.y)
        print "krr_out:", krr_out
        print "krr:", krr
        print "best parameters:", krr.best_params_
        #print "cv_results:", krr.cv_results_

        newX = self.predict2( indices, **kwargs)
        newY = krr.predict(newX)

        return newY


    def fraction_learn(self,fraction=0.8, skl=False):
        labda = self.labda
        sigma = self.sigma
        from sklearn.cross_validation import train_test_split
        X_train, X_test, y_train, y_test = train_test_split(self.coulombs,self.y, train_size = fraction)
        print "using a training set of", len(X_train), "molecules"
        print "using a test set of", len(X_test), "molecules"
        if skl:
            from sklearn.kernel_ridge import KernelRidge
            sigma = get_sigma( self.coulombs )
            gamma = 1. / ( sigma )
            print "sigma:", sigma
            if self.kerneltype == 'gaussian':
                kernel = 'rbf'
            elif self.kerneltype == 'laplacian':
                #from sklearn.metrics.pairwise import laplacian_kernel
                #kernel = laplacian_kernel
                kernel = 'laplacian'
            else: kernel = self.kerneltype
            clf = KernelRidge( alpha = self.labda, kernel=kernel , gamma = gamma)
            #clf = KernelRidge( alpha = self.labda, kernel=my_kernel3 , kernel_params = { 'factor': factor , 'kerneltype': self.kerneltype } )
            print "clf:", clf
            clf_out = clf.fit( X_train, y_train )
            y_train_pred = clf.predict( X_train )
            print "clf_out:", clf_out
            y_test_pred = clf.predict( X_test )
        else:
            kernel_train= my_kernel2( X_train, X_train, sigma=sigma, kerneltype = self.kerneltype )
            alpha       = self.solver( kernel_train, y_train , labda = labda)
            kernel_test = my_kernel2( X_test , X_train, sigma=sigma, kerneltype = self.kerneltype )
            kernel_test = my_kernel2( X_test , X_train, sigma=sigma, kerneltype = self.kerneltype )
            y_test_pred = self.testnew_no_skl( kernel_test , alpha = alpha )
        score = RMSE( y_test_pred, y_test )
        print "pearsonr:", pearsonr( y_test_pred, y_test)
        if True:
            plot_error( y_test_pred, y_test, 'bo', alpha= 0.5  )
            plot_error( y_train_pred, y_train, 'ro')
            plt.show()
            KDE( y_test_pred, y_test )

        return score

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
        n, bins, patches = plt.hist(scores, 50, normed=1, facecolor='green', alpha=0.75)
        mu = np.mean(scores)
        std = np.std(scores)
        print "mu:", mu
        print "std:", std
        y = mlab.normpdf( bins,mu,std)
        l = ax.plot(bins, y, 'r--', linewidth=1)
        f_y = mlab.normpdf( bins,f_mu,f_std)
        f_l = ax.plot(bins, f_y, 'b--', linewidth=1)
        plt.show()
    return


def my_kernel(X1, X2, sigma=1e2):
    def distance(C1, C2, type='norm3', kerneltype='gaussian'):
        '''calculate euclidian distance between two coulomb matrices'''
        def dist(arg1,arg2):
            l = min( ( len(arg1) , len(arg2) ) ) #minimum length of both
            d = np.sqrt( sum( [ (arg1[i] - arg2[i])**2 for i in range(l) ] ) )
            return d

        l1 = len(C1)
        l2 = len(C2)
        if l1 < l2:
            D = - C2.copy()
            D[:l1] += C1
        else:
            D = C1.copy()
            D[:l2] -= C2
        if kerneltype == "gaussian":
            total = np.sum( D**2 )
        elif kerneltype == "laplacian":
            total = np.sum( np.abs(D) )
        return np.sqrt(total)
    l1 = len(X1)
    l2 = len(X2)
    K = np.zeros([l1,l2])
    print "Kshape:", K.shape
    #for i in bar(range(l)):
    for i in range(l1):
        for j in range(l2): #for i in range(0) gives []
            K[i][j] = np.exp( - distance(X1[i], X2[j]) / ( 2 * sigma**2) )
        print "#",
    return K

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

def my_kernel2(X1, X2, sigma=1e2, kerneltype='gaussian'):
    l1 = len(X1) # n coulomb matrices in X1
    l2 = len(X2) # n coulomb matrices in X2
    K = np.zeros([l1,l2])
    print "Kshape:", K.shape
    # calculate denominator:
    if kerneltype == 'gaussian':
        factor = 1. / ( 2. * (sigma**2) )
    elif kerneltype == 'laplacian':
        factor = 1. / sigma
    for i in range(l1):
        for j in range(l2): #for i in range(0) gives []
            d = distance(X1[i], X2[j])
            K[i][j] = np.exp( - d * factor )
        #if debug:  print "#", d, K[i][j], factor
        #else:      print '#',
        print "#",
    return K

n_k3=0.0
def my_kernel3(C1,C2, factor, kerneltype='gaussian'):
    ''' Coulomb matrices as input instead of Xs '''
    d = distance(C1,C2)
    kij = np.exp( - distance(C1, C2) * factor )
    if debug:
        global n_k3
        n_k3+=1.0
        if n_k3/10. == int(n_k3/10):
            print "#",
            #print kij,
            #print "c1:", C1.shape
            #print "C1:", sprint(100,C1)
            #print "c2:", C2.shape
            #print "d:", d
            #time.sleep(1)
    return kij

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

def bf(U,y):
    l = len(U)
    #FORWARD
    alpha = np.zeros([l])
    for i in range(l):
        v = y[i]
        if not i==0:
            for j in range(i-1):
                v = v - U[j][i]*alpha[j]
        alpha[i] = v / U[i][i]
    #BACKWARD
    for i in range(l,0,-1):
        v = alpha[i]
        for j in range(l,1,-1):
            v = v - u[i][j]*alpha[j]
        alpha[i] = v / u[i][i] 
    return alpha

def cholesky(A):
    return np.linalg.cholesky(A)

def ML(zmas,y):
    '''machine learning on molecules represented by the zmas and training data in y vector '''
    xyzs = tuple( toxyz(item) for item in zmas )
    coulombs = tuple( coulomb(item) for item in xyzs)
    K = kernel(coulombs)
    alpha = solver(K,y)
    return alpha,coulombs

def plotmat(mat, log=1):
    if log:
        mat = np.log(mat)
    import matplotlib.pyplot as plt
    plt.matshow(mat)
    #ax = plt.gca()
    plt.colorbar()
    plt.show()
    return

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




######################################################################################
def machinelearning3(*args,**kwargs):
    import concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        result = executor.submit(machinelearning2, *args, **kwargs).result()
    return result

@log_io()
def machinelearning2(indices=[],table=[],sigma=1e5, labda= 1e-5, printlevel=1,**kwargs):
    ''' or this function will be called by CINDES'''
    #from converter import Converter
    from CINDES4.utils.converter import Converter
    converter = Converter()
    kwargs['converter'] = converter

    my_ML = MachineLearning(type='norm4',table=table,**kwargs)

    if False: #precomputed gram matrix
        my_ML.learn1()
        new_y = my_ML.predict1(indices,**kwargs)
    elif True: # "rbf"
        if False: # CV and gridsearch:
            new_y = my_ML.GridSearch(indices,**kwargs)
        elif False: # CV and gridsearch 2
            sigmas = np.logspace(0,15,15)
            labdas = np.logspace(-10,0,10)
            totalscores = search_grid( my_ML, sigmas=sigmas, labdas=labdas)
            print "totalscores:", totalscores
            raise SystemExit('stop')
        my_ML.learn()
        new_y = my_ML.predict(indices,**kwargs)
    else:    # kernel as implemented or as callable
        my_ML.learn() # learns on training set self.X , self.y gets alpha coefficients
        new_y = my_ML.predict(indices,**kwargs)
    print "new_y:", new_y
    return new_y

@log_io()
def normal_machinelearning(indices=[], table=[], sigma=1e6, labda= 1e-6, printlevel=1, **kwargs):
    ''' function to do ML without sklearn '''
    from CINDES4.utils.converter import Converter
    #from converter import Converter
    converter = Converter()
    kwargs['converter'] = converter

    # set self.coulombs and self.y and self.kerneltype, self.sigma, self.labda, self.type('norm4')
    # norm4 means coulombs are row/column sorted based on their norm > padded with zeros > flattened
    # @__init__ immediately setXY is ran. 
    my_ML = MachineLearning(type='norm4',table=table,sigma= sigma, labda= labda, **kwargs)

    if False:
        score, std = my_ML.CV(sigma=sigma, labda = labda)
        print "score:", score
    elif False:
        ''' Grid Search '''
        sigmas = np.logspace(  8,20, 7)
        labdas = np.logspace(-12, 0, 7)
        totalscores = search_grid( my_ML, sigmas=sigmas, labdas = labdas )
        print "totalscores:", totalscores
        raise SystemExit('stop')

    # learns self.alpha is now constructed
    my_ML.learn_no_skl()

    # predict
    newy = my_ML.predict_no_skl(indices,**kwargs)

    # #
    print "newy:", newy

    #
    return newy

############### COPY of for call as __main__ ######

@log_io()
def Amachinelearning2(indices=[],table=[],sigma=1e4, labda= 1., printlevel=1,fraction=0.5, kernel='gaussian', descriptor='norm4', **kwargs):
    ''' or this function will be called by CINDES'''
    #from converter import Converter
    from CINDES4.utils.converter import Converter
    converter = Converter()
    kwargs['converter'] = converter
    if descriptor in ['bob','BoB']:
        global BoB
        BoB = True
    else: BoB = False
    my_ML = MachineLearning(type=descriptor,table=table,kerneltype= kernel, **kwargs)
    if args.neural and args.fraction:
        print "X shape:", my_ML.coulombs.shape
        print "Y shape:", my_ML.y.shape
        print "X[0] shape:", my_ML.coulombs[0].shape
        import neural
        new_y = neural.main(X=my_ML.coulombs, y = my_ML.y , fraction = fraction)
        raise SystemExit('stop')


    gridsearch = False
    crossval   = False
    if gridsearch:
        pass
    else:
        if crossval:
            score,std = my_ML.CV( sigma=sigma, labda=labda, skl=False)
            print "score,std:", score, std
        else:
            fraction = fraction
            score = my_ML.fraction_learn(fraction=fraction, skl=True)
            print "score:", score
    return


if __name__=='__main__':
    class Unbuffered(object):
        def __init__(self,stream):
            self.stream = stream
        def write(self,data):
            self.stream.write(data)
            self.stream.flush()
        def __getattr__(self,attr):
            return getattr(self.stream, attr)
    sys.stdout = Unbuffered(sys.stdout)
    import argparse
    parser = argparse.ArgumentParser(description="reads cycles data stored in cyclesinfo")
    parser.add_argument("-i","--interactive",action="store_true",help="to be implemented")
    parser.add_argument("-p","--plot",action="store_true",help="make a property vs property plot of the data")
    parser.add_argument("-a","--anatrain",action="store_true",help="analyze and make a property vs property plot of the training data")
    parser.add_argument("-N","--neural",action="store_true",help="going to use Neural Networks")
    parser.add_argument("-S","--use_sklearn",action="store_true",help="analyze and make a property vs property plot of the training data")
    parser.add_argument("-T","--timer",action="store_true",help="perform some time analyses")
    parser.add_argument("-k","--kernel",action="store",type = str,default='gaussian',help="which kernel to use: (laplacian, gaussian)")
    parser.add_argument("-d","--descriptor",action="store",type = str,default='norm4',help="which descriptor to use: (BoB, Coulomb(norm1/norm3/norm4))")
    parser.add_argument("-r","--random",action="store_true",help="use a randomly selected test set and training set")
    parser.add_argument("-s","--sigma",action="store",nargs='?',type=float,default=1.e2,const=1e7,help="do a sigma default 1e2 KRR")
    parser.add_argument("-B","--bandwidth",action="store",nargs='?',type=float,default=0.5,const=1e7,help="do a sigma default 1e2 KRR")
    parser.add_argument("-c","--cutoff",nargs=2, type = float, help="cutoff values min max")
    parser.add_argument("-l","--labda",action="store",nargs='?',type=float,default=1.e-5,const=1e-5,help="do a labda default 1e-5 KRR")
    parser.add_argument("-f","--fraction",action="store",nargs='?',type=float,default=1,const=1,help="between 0-1 use this fraction as training set")
    args=parser.parse_args()
    # get a test c,a,p

    #if args.timer:
    if True:
        print "use sklearn:", args.use_sklearn
        from timer import Timer
        with Timer() as t:
            Amachinelearning2( sigma = args.sigma,
                               labda = args.labda,
                               fraction = args.fraction,
                               descriptor = args.descriptor,
                               kernel = args.kernel )
        print "=> elapsed learning5: %s s" % t.secs
    #else:
    #        Amachinelearning2( sigma = args.sigma,
    #                           labda = args.labda,
    #                           fraction = args.fraction,
    #                           descriptor = args.descriptor,
    #                           kernel = args.kernel )
    print "DONE LEARNING.PY"
    # load table.xyz


