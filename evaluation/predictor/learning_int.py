#!/bin/env python
'''module for machine learning for CINDES2.py program'''
debug=1
#import pyximport; pyximport.install()
#import cython_sum

#from writings import log_io, sprint, print_title
import logging
import sys
import pickle
import math as m
import gc
import numpy as np

#for a progres bar:
import time
#import progressbar
#bar = progressbar.ProgressBar()

from copy import deepcopy
from pprint import pprint
#from converter import Converter
#import construction as zcon

from CINDES.utils.writings import log_io, sprint, print_title
from CINDES.utils.converter import Converter
from CINDES.evaluation import construction as zcon


class MachineLearning(object):
    ''' Class for making training set / kernel / coulomb / predictions etc. '''
    def __init__(self,type='normal',kerneltype='laplacian'):
        self.name = name
        self.type = type
        self.kerneltype = kerneltype
        print_title("As a kernel: "+kerneltype+" is used",outline='l',signator='k',newlines=True)


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
                    y.append( float(line[0]) )   # 0: gap, 1: homo 2: lumo 3: Etotal
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
            self.y = y
        else:
            self.xyzs_t = xyzs
        return

    #@profile
    def coulomb(self, xyz):
        l = len(xyz)
        #print "xyz:", xyz
        C = np.zeros([l,l])
        alg = 2
        if alg == 1:
             for i in range(l):
                 #for j in range(l):
                 #    if i==j:
                 #        C[i][i]= 0.5*xyz[i][2]**2 
                 #    else:
                 #        t = xyz[i][2] * xyz[j][2]
                 #        xd = xyz[i][1][0] - xyz[j][1][0]
                 #        yd = xyz[i][1][1] - xyz[j][1][1]
                 #        zd = xyz[i][1][2] - xyz[j][1][2]
                 #        n = np.sqrt( xd**2 + yd**2 + zd**2 )
                 #        C[i][j] = t/n
                 for j in range(i+1):
                     if i==j:
                         #C[i][i]= 0.5*xyz[i][2]**2 
                         C[i][i]= 0.5*xyz[i][2]**(2.4)
                         
                     else:
                         t = xyz[i][2] * xyz[j][2]
                         xd = xyz[i][1][0] - xyz[j][1][0]
                         yd = xyz[i][1][1] - xyz[j][1][1]
                         zd = xyz[i][1][2] - xyz[j][1][2]
                         n = np.sqrt( xd**2 + yd**2 + zd**2 )
                         C[i][j] = t/n
                         C[j][i] = C[i][j]
        elif alg == 2:
            for i in xrange(l):
                for j in range(i+1):
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
            # return a lower triangular matrix of the coulomb matrix. 
            return C[ np.tril_indices(l) ]
        else:
            assert self.type=='normal'
            return C

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

        
    #@profile
    def distance(self, C1, C2):
        '''calculate euclidian distance between two coulomb matrices'''
        #take eigenvalues of symmetric hermitian matrices. Sort them in descending order
        def dist(arg1,arg2):
            l = min( ( len(arg1) , len(arg2) ) ) #minimum length of both
            d = np.sqrt( sum( [ (arg1[i] - arg2[i])**2 for i in range(l) ] ) )
            # test contribution of rest values. 

            return d
        if self.type=='eig':
            arg1 = np.linalg.eigvalsh(C1)[::-1]
            arg2 = np.linalg.eigvalsh(C2)[::-1]
            return dist(arg1,arg2)
        elif self.type=='norm1':
            # assume already sorted based on norm
            arg1 = np.linalg.norm(C1,axis=0)
            arg2 = np.linalg.norm(C2,axis=0)
            return dist(arg1,arg2)
        elif self.type in ['norm2']:
            #assume that the coulombs are already sorted 1d vectors of norms
            # or that they are just a linear vector of coulomb matrix entries
            return dist(C1,C2)
        elif self.type == 'norm3':
            l1 = len(C1)
            l2 = len(C2)
            alg = 2
            if alg == 1:
                #distance between two coulomb vectors. 
                d = 0
                for i in xrange( max((l1,l2)) ):
                    try:
                        d += ( C1[i] - C2[i] )**2
                    except IndexError:
                        if l1>l2:
                           d += C1[i] ** 2
                        if l1<l2:
                           d += C2[i] ** 2
                dis = np.sqrt(d)
                #print "distance: ", dis,
                return dis
            if alg == 2:
                if l1 < l2:
                    D = - C2.copy()
                    D[:l1] += C1
                else:
                    D = C1.copy()
                    D[:l2] -= C2
                if self.kerneltype == "gaussian":
                    total = np.sum( D**2 )
                elif self.kerneltype == "laplacian":
                    total = np.sum( np.abs(D) )
                return np.sqrt(total)
        else:
            raise SystemExit('No valid distance calculation specified')
    
    @log_io()
    def solver(self,labda=1):
        l = len(self.kernel)
        I = np.identity(l)
        Ka = self.kernel + labda * I
        #print "labda*I", labda*I
        #U = cholesky(Ka)
        #alpha = bf(U,y)
        if hasattr(self,'training_y'):
            print "has a training set. so use self.training_y instead of self.y"
            ny = np.array(self.training_y)
        else:
            ny = np.array(self.y)
        #print "ny:", ny
        alg = 7
        if alg == 1:
            alpha = np.dot( np.linalg.inv(Ka) , (ny.T) )
        elif alg == 2:
            from scipy import linalg
            alpha = np.dot( linalg.inv(Ka) , (ny.T) )
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
        return

    #@profile
    def testnew(self,coulombs):
        outtest = []
        for coulombje in coulombs:
        #for j in bar(  xrange( len(self.coulombs_t) )  ):
            ans = 0
            for i in xrange(len(self.kernel)):
                #ans += self.alpha[i] * np.exp( self.distance( self.coulombs_t[j], self.coulombs[i] ) / ( 2 * self.sigma **2) )
                ans += self.alpha[i] * np.exp( self.distance( coulombje, self.coulombs[i] ) / ( 2 * self.sigma **2) )
            outtest.append(ans)
            #print "$",
        return outtest

    def predict(self, inputfile='chemspace.xyz'):
        # get input from inputfile. stored in self.xyzs_t
        self.get_input(data=0, inputfile=inputfile)
        # get coulombs
        self.coulombs_t = tuple( self.coulomb(item) for item in self.xyzs_t )
        # calculate 
        self.y_t = self.testnew(self.coulombs_t)
        return self.y_t       

    def predict2(self,indices,converter,**kwargs):
        mat_inds = tuple( contozma(zcon.indtocon(item),**kwargs) for item in indices)
        self.xyzs_t = [ zmatoxyz(converter, item) for item in mat_inds ]
        # get coulombs
        self.coulombs_t = tuple( self.coulomb(item) for item in self.xyzs_t )
        # calculate 
        self.y_t = self.testnew(self.coulombs_t)
        return self.y_t

    def predict_conf(self,ind):
        #convert ind to zma 

        #convert zma to xyz

        #convert xyz to coulomb

        #testnew
        pass

    def ML(self):
        from converter import Converter
        converter = Converter()
        # get input from inputfile table.xyz
        self.get_input(converter)
        # calculate all coulomb matrices
        self.coulombs = tuple( self.coulomb(item) for item in self.xyzs)
        # calculate kernel
        self.get_kernel()
        # calculate alpha coefficients
        self.solver()
        return 

    def ML0(self,fraction,sigma,labda):
        print "sigma: ", str(sigma), ' ', 'labda: ', str(labda)
        percentage = fraction
        size = len(self.xyzs)

        #select data sets
        n = int(fraction*size)
        if args.random:
            from sklearn.model_selection import train_test_split
            self.training_xyzs, self.test_xyzs, self.training_y, self.test_y = train_test_split(self.xyzs,self.y, train_size = fraction)
        else:
            #training set:
            self.training_xyzs = self.xyzs[: n ]
            self.training_y = self.y[ : n ]
            # test set:
            self.test_xyzs = self.xyzs[n:]
            self.test_y = self.y[n:]
        #print "i will use ", len(self.test_xyzs), " test molecules"
        print "using a training set of", len(self.training_xyzs), "molecules"
        print "using a test set of", len(self.test_xyzs), "molecules"

        # get coulombs of training set
        self.coulombs = tuple( self.coulomb(item) for item in self.training_xyzs)
        # make the kernel from self.coulombs. 
        if args.timer:
            with Timer() as t:
                self.get_kernel(sigma=sigma)
            print "=> elapsed make kernel data: %s s" % t.secs
        else:
            self.get_kernel(sigma=sigma)
        print "kernel is made. first few entries of self.kernel[0] look like:",
        sprint(20,self.kernel[0])        
        # solve Ka = y here. get alpha. 
        self.solver(labda=labda)
        print "solver done: first elements of alpha:"
        sprint(5,self.alpha)

        # To test the training set:
        if args.anatrain:
            print "TRAINING SET RESULTS:"
            trainingresults = self.testnew(self.coulombs)
            mae = MAE(self.training_y, trainingresults)
            rmse = RMSE(self.training_y, trainingresults)
            print "mean absolute Error is:", mae
            print "RMS Error is:", rmse

        #now predictions 
        self.coulombs_t = tuple( self.coulomb(item) for item in self.test_xyzs)
        # testnews uses the kernel to get new predicted values of ys for the testset. 
        # uses also self.coulombs_t 
        # the results go to self.y_t
        if args.timer:
            with Timer() as t:
                self.y_t = self.testnew(self.coulombs_t) 
            print "=> elapsed test new data: %s s" % t.secs
        else:
            self.y_t = self.testnew(self.coulombs_t) 
 
        #now there is a self.y_t = predict en self.test_y is real value
        print "test done:"
        #sprint(5,self.y_t)
        print "TESTING RESULTS:"
        mae = MAE(self.y_t,self.test_y)
        rmse = RMSE(self.y_t,self.test_y)
        #print(self.y_t)
        scores = np.array([ i - j for i,j in zip(self.y_t,self.test_y) ])
        X_plot = np.linspace(-5, 10, 1000)
        print "\nscores:"
        sprint(10, scores)

        if True:
            import matplotlib.pyplot as plt
            from sklearn.neighbors import KernelDensity
            X_plot = np.linspace(-5, 5, 1000)[:, np.newaxis]
            fig, ax = plt.subplots()
            kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(scores[:,np.newaxis])
            log_dens = kde.score_samples(X_plot)
            print "logdens:"
            sprint(10,log_dens)
            ax.plot(X_plot, np.exp(log_dens), 'c-')
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
        if True:
            outname = 'logdens' + str(fraction*10) + str(sigma) + str(labda) + '.b'
            import pickle
            with open(outname,'wb') as fid:
                pickle.dump(log_dens,fid)
                #pickle.dump(


        print "mean absolute Error is:", mae
        print "RMS Error is:", rmse
        print "compare a few"
        sprint(10, self.y_t, self.test_y)

        if args.plot:
            import matplotlib.pyplot as plt
            # training = x-as / test is y 
            if args.anatrain:
                plt.plot( self.training_y, trainingresults, 'ro')
                plt.plot( self.y_t, self.test_y, 'bo',alpha=0.5)
            else:
                plt.plot( self.y_t, self.test_y, 'bo')
            plt.title( ' real values vs tested values with ML ' )
            plt.show()

        return rmse, mae

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

def indtoint(index, array):
    conf = zcon.indtocon(index)
    intl = zcon.contoint(conf,array)
    return intl

def get_table(inputfile='table_unbiased'):
    with open(inputfile,'rb') as f:
        table = pickle.load(f)
    return table

def distance_int(x1, x2):
    l = len(x1)
    #laplacian
    d = sum( [ abs(x1[i] - x2[i]) for i in range(l) ] )
    return d

def get_kernel_int(X, sigma=1):
    l = len(X)
    K = np.zeros([l,l])
    for i in range(l):
        for j in range(i+1): #for i in range(0) gives []
            if not i==j:
                K[i][j] = np.exp( - distance_int(X[i], X[j]) / sigma )
                K[j][i] = K[i][j]
            else:
                K[i][j] = 1.0
        print "#",
    return K

def kernel_callable(x1,x2, sigma=1):
    kij = np.exp( - distance_int( x1, x2 ) / sigma )
    return kij

def solver_int(K,Y,labda=1):
    l = len(K)
    I = np.identity(l)
    Ka = K + labda * I
    alg = 7
    if alg == 1:
        alpha = np.dot( np.linalg.inv(Ka) , (Y.T) )
    elif alg == 7:
        from scipy import linalg
        lu = linalg.lu_factor(Ka)
        alpha = linalg.lu_solve(lu, Y)
    print "alphashape:", alpha.shape
    return alpha

def testnew_int(alpha,X,newX,sigma):
    outtest = []
    for new in newX:
        ans = 0
        for i in xrange(len(X)):
            ans += alpha[i] * np.exp( - distance_int( X[i], new ) / sigma )
        outtest.append(ans)
        #print "$",
    return outtest

def plot(data1,data2,*args,**kwargs):
    import matplotlib.pyplot as plt
    import seaborn
    plt.plot(data1,data2,*args,**kwargs)
    plt.show()
    return 

def new_procedure(sigma=1, labda=1):
    import INDES
    table = get_table()[:15]
    param, array = INDES.read_input('INPUTBC')
    X = [ indtoint(item[0],array) for item in table ]
    print "X:",
    sprint(10,X)
    Y = np.asarray([ item[2] for item in table ])
    print "Y:",
    sprint(10,Y)
    K = get_kernel_int(X,sigma = sigma )
    print "K[1]:", K[1]
    #sprint(10,K[1])
    alpha = solver_int(K,Y, labda = labda)
    print "alpha:", alpha
    newX = X[:]
    newY  = testnew_int(alpha, X, newX, sigma=sigma)
    print "newY:",
    sprint(10,newY)
    R = RMSE(newY,Y)
    print "R=", R
    plot(newY,Y,'or')
    return

##### INDES CALL #####
@log_io()
def learn_int_procedure(table, indices, array, sigma=1e2, labda=1e-4, **kwargs):
    X = [ indtoint(item[0],array) for item in table ]
    print "X:",
    sprint(10,X)
    Y = np.asarray([ item[2] for item in table ])
    print "Y:",
    sprint(10,Y)
    K = get_kernel_int(X,sigma = sigma )
    print "K[1]:", K[1]
    #sprint(10,K[1])
    alpha = solver_int(K,Y, labda = labda)
    print "alpha:", alpha
    if debug:
        print "indices:", indices
    newX = [ indtoint(item, array) for item in indices ]
    newY  = testnew_int(alpha, X, newX, sigma=sigma)
    print "newY", newY
    #raise SystemExit('stop')
    return newY

@log_io()
def learn_int_skl_procedure(table, indices, array, sigma=1e4, labda=1e-4, **kwargs):
    X = [ indtoint(item[0],array) for item in table ]
    print "X:",
    sprint(10,X)
    Y = np.asarray([ item[2] for item in table ])
    print "Y:",
    sprint(10,Y)

    newX = [ indtoint(item, array) for item in indices ]

    from sklearn.kernel_ridge import KernelRidge
    if False:
        newY = skl_example(X,Y,newX)
        #print "svY:", svY
    elif False:
        gamma = 1. / ( 2 * sigma**2 )
        clf = KernelRidge(alpha = labda, kernel='rbf', gamma= gamma)
        clf_out = clf.fit(X,Y)
        print "clf_out:", clf_out
        newY = clf.predict(newX)
    else:
        clf = KernelRidge(alpha = labda, kernel = kernel_callable, kernel_params = { 'sigma':sigma } )
        clf_out = clf.fit(X,Y)
        score   = clf.score(X,Y)
        print "KRR score:", score
        newY = clf.predict(newX)
    if False:
        from sklearn.svm import SVR
        from sklearn import svm
        from sklearn import preprocessing
        from sklearn.model_selection import GridSearchCV
        scaler = preprocessing.StandardScaler().fit(X)
        X_scaled = scaler.transform(X)
        newX_scaled = scaler.transform(newX)
        #svr = SVR(        C     = 10   , kernel = kernel_callable, epsilon = 0.01 )
        svr = svm.SVR(kernel = 'rbf' ,
                      gamma  =   1e-8,
                      C      =   1.0 ,
                      epsilon=   0.1 )
        svr = GridSearchCV( SVR( kernel='rbf', gamma=0.1, epsilon = 0.1 ), cv=5,
                            param_grid = { 'C' : np.logspace(-10,10,5),
                                           'gamma': np.logspace(-20,20,5) ,
                                           'epsilon': [ 1, 0.1] } )
        svr_out = svr.fit(X_scaled,Y)
        print "best SVR params:", svr.best_params_
        svrY = svr.predict(newX_scaled)
        print "scores SVM:", svr.score(X_scaled,Y)
        print "svrY", svrY
        raise SystemExit('stop')
    print "newY", newY
    return newY

def skl_example(X,Y,newX):
    from sklearn.svm import SVR
    from sklearn.model_selection import GridSearchCV
    from sklearn.model_selection import learning_curve
    from sklearn.kernel_ridge import KernelRidge
    sigma = 1.
    if True:
        #svr = GridSearchCV( SVR( kernel= 'rbf', gamma= 0.1) ,
        #                    cv    = 5,
        #                    param_grid ={'C'     : np.logspace(-4, 4,10) ,
        #                                 'gamma' : np.logspace(-2, 2, 5) ,
        #                                 'epsilon':np.logspace(-2, 2, 5) })
        # unfortunately GridSearchCV cannot use a custom kernel
        krr = GridSearchCV( KernelRidge( kernel= kernel_callable  ),
                            cv    = 5,
                            param_grid = {'alpha' : [1e4, 1e2, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5,1e-6 ],
                                          'gamma' : np.logspace(-6, 6, 10) } )
    else:
        svr = SVR(         kernel = 'poly',degree=3, gamma = 0.01  , C = 10, epsilon = 0.01 )
        krr = KernelRidge( kernel = 'rbf', gamma = 0.0005, alpha = 1e-2 )
    #svr.fit(X, Y)
    krr.fit(X, Y)
    #svr_ratio = svr.best_estimator_
    #print "sv_ratio:", svr_ratio
    #krr_ratio = krr.best_estimator_
    #print "krr_ratio:", krr_ratio
    print "krr:", krr
    #print "svr:", svr

    #svY = svr.predict(newX)
    krY = krr.predict(newX)
    return  krY


################# RUN AS __MAIN___ ########################

def learn_int_skl_procedure(table, indices, array, sigma=1e4, labda=1e-4, **kwargs):
    X = [ indtoint(item[0],array) for item in table ]
    print "X:",
    sprint(10,X)
    Y = np.asarray([ item[2] for item in table ])
    print "Y:",
    sprint(10,Y)

    newX = [ indtoint(item, array) for item in indices ]

    from sklearn.kernel_ridge import KernelRidge
    if False:
        gamma = 1. / ( 2 * sigma**2 )
        clf = KernelRidge(alpha = labda, kernel='rbf', gamma= gamma)
        clf_out = clf.fit(X,Y)
        print "clf_out:", clf_out
        newY = clf.predict(newX)
    else:
        clf = KernelRidge(alpha = labda, kernel = kernel_callable, kernel_params = { 'sigma':sigma } )
        clf_out = clf.fit(X,Y)
        score   = clf.score(X,Y)
        print "KRR score:", score
        newY = clf.predict(newX)
    if False:
        from sklearn.svm import SVR
        from sklearn import svm
        from sklearn import preprocessing
        from sklearn.model_selection import GridSearchCV
        scaler = preprocessing.StandardScaler().fit(X)
        X_scaled = scaler.transform(X)
        newX_scaled = scaler.transform(newX)
        #svr = SVR(        C     = 10   , kernel = kernel_callable, epsilon = 0.01 )
        svr = svm.SVR(kernel = 'rbf' ,
                      gamma  =   1e-8,
                      C      =   1.0 ,
                      epsilon=   0.1 )
        svr = GridSearchCV( SVR( kernel='rbf', gamma=0.1, epsilon = 0.1 ), cv=5,
                            param_grid = { 'C' : np.logspace(-10,10,5),
                                           'gamma': np.logspace(-20,20,5) ,
                                           'epsilon': [ 1, 0.1] } )
        svr_out = svr.fit(X_scaled,Y)
        print "best SVR params:", svr.best_params_
        svrY = svr.predict(newX_scaled)
        print "scores SVM:", svr.score(X_scaled,Y)
        print "svrY", svrY
        raise SystemExit('stop')
    print "newY", newY
    return newY

######################################################################################


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
    parser.add_argument("-T","--timer",action="store_true",help="perform some time analyses")
    parser.add_argument("-k","--kernel",action="store",type = str,default='gaussian',help="which kernel to use: (laplacian, gaussian)")
    #parser.add_argument("-K","--kernel2",action="store",nargs=1,default='gaussian', type=str, help="kerneltype")
    parser.add_argument("-r","--random",action="store_true",help="use a randomly selected test set and training set")
    parser.add_argument("-s","--sigma",action="store",nargs='?',type=float,default=1.e2,const=1e2,help="do a sigma default 1e2 KRR")
    parser.add_argument("-c","--cutoff",nargs=2, type = float, help="cutoff values min max")
    parser.add_argument("-l","--labda",action="store",nargs='?',type=float,default=1.e-5,const=1e-5,help="do a labda default 1e-5 KRR")
    parser.add_argument("-f","--fraction",action="store",nargs='?',type=float,default=1,const=1,help="between 0-1 use this fraction as training set")
    parser.add_argument('file',help="a pickled tablebin file")
    args=parser.parse_args()

    import pickle
    with open(args.file, 'rb') as f:
        table = pickle.load(f)



    new_procedure(sigma=args.sigma, labda=args.labda)

