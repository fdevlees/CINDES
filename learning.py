#!/bin/env python
'''module for machine learning for CINDES2.py program'''

#import pyximport; pyximport.install()
#import cython_sum

from writings import log_io, sprint, print_title
import logging
import sys
import pickle
import math as m
import gc
import numpy as np

#for a progres bar:
import time
import progressbar
bar = progressbar.ProgressBar()

from copy import deepcopy
from pprint import pprint
from converter import Converter
import construction as zcon

class MachineLearning(object):
    ''' Class for making training set / kernel / coulomb / predictions etc. '''
    def __init__(self,name,type='normal',kerneltype='laplacian'):
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
            print "$",
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
            from sklearn.cross_validation import train_test_split
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

    def ML2(self,table,core,active,passive,converter, **kwargs):
        sigma = 1e2
        labda = 1e-5
        print "sigma: ", str(sigma), ' ', 'labda: ', str(labda)
        # get input from inputfile table
        self.y =  np.fromiter((item[1] for item in table ),np.float)
        mats =  tuple( contozma(zcon.indtocon(item[0]),core,active,passive) for item in table )
        self.xyzs = [ zmatoxyz(converter,item) for item in mats ]
        # calculate all coulomb matrices
        self.coulombs = tuple( self.coulomb(item) for item in self.xyzs)
        print "the first entries of the first two Coulomb matrices: "
        if self.type == 'norm3':
            for n in [0,1]:
                l = len( self.coulombs[0] )
                k = 0
                i = 0
                while True:
                    if i>20: break
                    for j in range(i+1):
                        try:
                            print self.coulombs[n][k],
                            k+=1
                        except IndexError:
                            break
                    else:
                        print
                        i+=1
                        continue
                    break
        else:
            sprint(1,self.coulombs)
        print "end"
        # calculate kernel
        self.get_kernel(sigma=sigma)
        # calculate alpha coefficients
        self.solver(labda=labda)
        return 
     

def symsort(mat):
    indexlist = np.argsort(np.linalg.norm(mat,axis=1))[::-1]
    return mat[indexlist][:,indexlist]

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
    plt.show()
    return   

def zmatoxyz(a,mat):
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

######################################################################################
@log_io()
def machinelearning(indices,table=[],**kwargs):
    '''this function will be called by CINDES2.py'''
    # 1. The table has to be converged to a table bin. This is in the generator functions in CINDES2.py
    from converter import Converter
    converter = Converter()
    #generate a table.xyz file were all the xyzs of the table are generated. 
    generate1(converter,table,**kwargs)
    print "a table.xyz file is generated"
    #this table.xyz is used to make the alpha vector via machinelearning.
    my_ML = MachineLearning('name')
    my_ML.ML()
    #now there is an alpha argument of my_ML
    print "alpha coefficients are calculated"
    # now the indices has to be converted to xyz coordinates to. 
    conffile='confs.xyz'
    print "indices:", indices
    generate_xyz(indices=indices,converter=converter,outputfile=conffile,**kwargs)
    # there is a confs.xyz.
    new_y = my_ML.predict(inputfile=conffile)
    print new_y
    return

def machinelearning3(*args,**kwargs):
    import concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        result = executor.submit(machinelearning2, *args, **kwargs).result()
    return result


@log_io()
def machinelearning2(indices,table,printlevel=1,**kwargs):
    ''' or this function will be called by CINDES'''
    from converter import Converter
    converter = Converter()
    kwargs['converter'] = converter

    ## NB: the type here:
    my_ML = MachineLearning('name',type='norm3')

    #now alpha and kernel are constructed
    my_ML.ML2(table,**kwargs)
    if printlevel==1:
        print "alpha coefficients are calculated"
        sprint(5,my_ML.alpha)
        print "kernel[0:1]:"
        sprint(2,my_ML.kernel)
    # now the indices has to be converted to xyz coordinates to. 
    new_y = my_ML.predict2(indices,**kwargs)
    zzz = gc.collect()
    print "gc.collect():", zzz
    gc.DEBUG_LEAK
    return new_y

def MC_init(table=[],**kwargs):
    from converter import Converter
    converter = Converter()
    generate1(converter,table,**kwargs)
    print "a table.xyz file is generated"
    my_ML = MachineLearning('name','norm3')
    my_ML.ML()
    return my_ML

def MC_test_ind(ml, indices,**kwargs):
    from converter import Converter
    converter = Converter()
    conffile='confs.xyz'
    print "indices:", indices
    generate_xyz(indices=indices,converter=converter,outputfile=conffile,**kwargs)
    # there is a confs.xyz.
    new_y = ml.predict(inputfile=conffile)
    print new_y
    return new_y

def machinelearning4():
    '''old version. try machinelearning5'''
    from converter import Converter
    converter = Converter()
    kwargs=dict()
    kwargs['converter'] = converter

    ## NB: the type here:
    mlin = MachineLearning('name',type='norm2')
    mlin.get_input(converter=converter)
    #now we have my mlin.xyzs and mlin.y
    print "some y values of total set"
    sprint(5, mlin.y)
    # now we need to split it in two parts. say 10 % training set and 90 % testset.
    a = 0.8
    mlin.ML0(a)   
    print "DONE ML"
    return

def machinelearning5(sigma=1e4,labda=0,fraction=0.2,kerneltype='laplacian'):
    from converter import Converter
    converter = Converter()
    kwargs = dict()
    kwargs['converter'] = converter
    print "KERNEL:", kerneltype
    mlin = MachineLearning('name',type='norm3',kerneltype=kerneltype)
    mlin.get_input(converter=converter)
    #now we have my mlin.xyzs and mlin.y
    print "some y values of total set"
    sprint(5, mlin.y)
    rmse, mae = mlin.ML0(fraction = fraction,sigma=sigma,labda=labda)
    print "DONE ML"
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
    parser.add_argument("-T","--timer",action="store_true",help="perform some time analyses")
    parser.add_argument("-k","--kernel",action="store",type = str,default='gaussian',help="which kernel to use: (laplacian, gaussian)")
    #parser.add_argument("-K","--kernel2",action="store",nargs=1,default='gaussian', type=str, help="kerneltype")
    parser.add_argument("-r","--random",action="store_true",help="use a randomly selected test set and training set")
    parser.add_argument("-s","--sigma",action="store",nargs='?',type=float,default=1.e2,const=1e2,help="do a sigma default 1e2 KRR")
    parser.add_argument("-c","--cutoff",nargs=2, type = float, help="cutoff values min max")
    parser.add_argument("-l","--labda",action="store",nargs='?',type=float,default=1.e-5,const=1e-5,help="do a labda default 1e-5 KRR")
    parser.add_argument("-f","--fraction",action="store",nargs='?',type=float,default=1,const=1,help="between 0-1 use this fraction as training set")
    args=parser.parse_args()
    # get a test c,a,p
    print args.cutoff 

    if args.interactive:
        pass
    else:
        if args.timer:
            from timer import Timer
            with Timer() as t:
                machinelearning5(sigma=args.sigma,labda=args.labda, fraction=args.fraction)
            print "=> elapsed learning5: %s s" % t.secs
        else:
            #print "args.kernel",args.kernel2
            machinelearning5( sigma = args.sigma,
                              labda = args.labda, 
                           fraction = args.fraction,
                         kerneltype = args.kernel )
    print "DONE LEARNING.PY"
    # load table.xyz


