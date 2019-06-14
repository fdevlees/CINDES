#!/bin/env python
ttert = False # keep a distinction between tert apical and termedial. this will give double amount of coefficients
debug=True


########################
#####   IMPORTS    #####
########################
import pickle
from pprint import pprint
from re import findall, split
import numpy as np
from sklearn import linear_model, metrics
from sklearn.model_selection import train_test_split
from abc import ABCMeta, abstractmethod
import pandas as pd
from CINDES.utils.writings import log_io, sprint

from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

#################################
#####    BASIC FUNCTIONS    #####
#################################

def indtocon(index):
    #return [list(item) for item in index.split('_')]
    #return  [ findall('[A-Z][^A-Z]*',item) for item in index.split('_') ]
    return index.split('_')

def contoind(conf):
    return '_'.join([''.join(item) for item in conf])

def indtoint(index, array):
    import CINDES.INDES.construction as zcon
    conf = zcon.indtocon(index)
    intl = contoint(conf,array)
    return intl

def contoint(conf,array):
    ''' makes an integer list representation of conf '''
    inconf = [ site.index(group) for site,group in zip(array,conf) ]
    return inconf

def get_X_int(indices, array):
    X = [ indtoint(index, array) for index in indices ]
    #print "in get_X_int:", X
    return np.asarray(X)

#++++++++++++++++++++++++++
#++++     GLOBALS    ++++++
#++++++++++++++++++++++++++

eV=27.2113838

funcs = {'CCFFF': '$C-CF_3$',
         'CCHHH': '$C-CH_3$',
         'CCN': '$C-C\\equiv N$',
         'CCl': '$C-Cl$',
         'CF': '$C-F$',
         'CH': '$C-H$',
         'CNHH': '$C-NH_2$',
         'CNOO': '$C-NO_2$',
         'CCOOH': '$C-COOH$',
         'CO': '$C=O$',
         'COH': '$C-OH$',
         'CSH': '$C-SH$',
         'N': '$N$',
         'O': '$O$',
         'S': '$S$',
         'B': '$B$'}
seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S', 'B']
fseq= [ funcs[func] for func in seq ] #whahaha :)

###########################
####      CLASSES     #####
###########################

class Dataset(): #abstract data class

    def extract(self, table=None, **kwargs):
        #confs, data = self.readfile(cutoff=20.0)
        if table==[]:
            print "read table from tablebin"
            column=args.column
            print("column: ", column)
            self.confs, data = self.readfile(column=column,**kwargs)
            #print "self.confs[20:22]:", self.confs[19:22]
            self.Y = self.extractY(data,column=1)
        else:
            print "read table from call"
            inds = [ item[0] for item in table ]
            self.confs = [ indtocon(item) for item in inds]
            data = [ item[1:] for item in table ]
            column = args.column-1
            self.Y = self.extractY(data,column=column)
            #self.Y = 27.2113838*self.Y

        # here symmetry operations: input: args.symmetry | self.confs | self.Y | self.syms
        if args.symmetry==True:
            datadict = dict() #dictionary to avoid duplicates
            for conf, ytje in zip(self.confs, self.Y): #loop over configurations + values
                conf_np = np.array(conf) #convert to np.array for boolean indexing
                if debug: print conf_np
                for isomer in self.syms: #for each symmetrical isomer
                    newconf = conf_np[isomer] #create a new conf. 
                    newind = contoind(newconf) #make the conf an index to be used in the dict
                    datadict[newind]=ytje #update the dict. if already present nothing will happen. 
            self.inds, self.Y = zip( *datadict.items() ) #convert the dictionary back to lists. 
            self.confs = [ indtocon(index) for index in self.inds ] #convert the indices back to configurations. 
            if debug:
                print "len confs:", len(self.confs)
                print "len Y:", len(self.Y)
                print "len dict", len(datadict)
        # output: self.confs(updated) | self.Y(updated) | self.inds | 

        self.X = self.extractX(self.confs)
        #print "extraction succesfull"
        if args.verbose<0:
            print "X:",self.X.shape
            print "Y:",self.Y.shape

        test= 0
        if test==1:
            n = min(2,len(self.Y))
            sprint(n,self.X,self.Y)
        return

    def extract2(self):
        self.X2 = self.extract2DX(self.confs)
        #PHENALENE:
        #hits = np.sum(self.X2, axis=0).reshape([11,11]).astype(int)
        # DIAMONDOIDS:
        if ttert: hits = np.sum(self.X2, axis=0).reshape([12,15*2]).astype(int)
        else: hits = np.sum(self.X2, axis=0).reshape([12,15]).astype(int)
        return hits

    def analyze(self,data):
        '''analyzes the structure of the data file'''
        print "data analyzation"
        print "verbosity:", args.verbose
        ldat = len(data)
        print "data consists of {} elements".format(ldat)
        for i in range(ldat):
            datel = data[i]
            dcolumn = [ item[args.column] for item in datel ]
            print "element {} has {} datapoints".format(i,len(datel)),
            avg = np.mean(dcolumn)
            std = np.std(dcolumn)
            print "with average of {} and std of {}".format(avg,std)
        if args.verbose>2:
            for i,item in enumerate(data[0]):
                print i, item
        return

    def openfile(self,file='tablebin'):
        ruwdata = []
        with open(file,'rb') as fid:
          if True:
            while True:
                try:
                    ruwdata.append(pickle.load(fid))
                except EOFError:
                    break
        if args.analyseinput:
            self.analyze(ruwdata)
        if args.verbose>1:
            print "first 10 of ruwdata:"
            for i in range(10):
                print ruwdata[0][i]
        return ruwdata

    def readfile(self,cutoff=False,column=1,**kwargs):
        '''extracts the data from filename: name '''
        ruwdata = self.openfile(**kwargs)
        datar=ruwdata[-1][:] # last entry in tablebin
        #########
        if False:
            data = datar[:488] + datar[1616:]
        else:
            data = datar
        if args.verbose>1:
            print "first 10 of data:"
            sprint(10,data)
        ########
        if cutoff:
            print "cutoff applied 15 eV"
            cutoff = 15.0
            data = [ [ item[0], float(item[column])] for item in datar if float(item[1])<cutoff ] 
        if False:
            data = [ [ item[0], float(item[column])] for item in datar if not 'B' in item[0] ]
        self.data = data
        confs = [ indtocon(item[0]) for item in data ]
        indices = [ item[0] for item in data ]
        self.indices = indices
        self.confs = confs
        return confs,data   

    def extractY(self,data,column=0):
        Y = np.array([ item[column] for item in data ])
        return Y

    def write_X(self):
        '''print the coefficient matrix. can only be called when
        the extractX of the child is already used''' 
        with open('Amatrix','w') as fid:
           for item in self.X:
               fid.write(str(item))
               fid.write('\n')
        return

class Adamantane(Dataset):
    '''  TYPE OF MOLECULE SUBCLASS OF DATASET   '''

    #specific adamantane parameters
    ngps = (12,12,12,12,15,15,15,15,15,15) #for every instance this is same

    def __init__(self,*rgs,**kwargs):
        self.seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S']
        self.bonds = ( (0,9),(0,5),(0,7),
                     (3,9),(3,4),(3,8),
                     (2,6),(2,8),(2,7),
                     (1,4),(1,5),(1,6) )
        if args.symmetry:
            Adasym = [
             [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10],
             [ 1, 3, 4, 2, 7, 8, 9,10, 5, 6],
             [ 1, 4, 2, 3, 9,10, 5, 6, 7, 8],
             [ 3, 2, 4, 1, 6, 7, 5, 9,10, 8],
             [ 4, 2, 1, 3, 7, 5, 6,10, 8, 9],
             [ 4, 1, 3, 2, 6,10, 8, 9, 7, 5],
             [ 2, 4, 3, 1,10, 5, 9, 7, 8, 6],
             [ 3, 1, 2, 4,10, 8, 6, 7, 5, 9],
             [ 2, 3, 1, 4, 9, 7, 8, 6,10, 5],
             [ 4, 3, 2, 1, 8, 9, 7, 5, 6,10],
             [ 3, 4, 1, 2, 5, 9,10, 8, 6, 7],
             [ 2, 1, 4, 3, 8, 6,10, 5, 9, 7] ]
            self.syms = np.array(Adasym) -1
            #ADA example
            #conf = ['a','b','c','d','e','f','g','h','i','j']
            #conf_new = np.array(conf)[Asyms[3]] #for example
        return

    def extractX(self,confs):
        nC= len(confs)
        if args.verbose>2:
            for i in range(20):
                print confs[i]
            print "seq:", seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        site0 = np.zeros([nC,12])
        site1 = np.zeros([nC,12])
        site2 = np.zeros([nC,12])
        site3 = np.zeros([nC,12])
        site4 = np.zeros([nC,15])
        site5 = np.zeros([nC,15])
        site6 = np.zeros([nC,15])
        site7 = np.zeros([nC,15])
        site8 = np.zeros([nC,15])
        site9 = np.zeros([nC,15])
        if args.verbose>1:
            print "dim site0", site0.shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                if group=='CNOO60':
                    group='CNOO'
                j = seq.index(group) #find the index of the group of that sequence
                if i == 0:
                    site0[k,j]=1 #index is number of configuration , number of group
                if i == 1:
                    site1[k,j]=1
                if i == 2:
                    site2[k,j]=1
                if i == 3:
                    site3[k,j]=1
                if i == 4:
                    site4[k,j]=1
                if i == 5:
                    site5[k,j]=1
                if i == 6:
                    site6[k,j]=1
                if i == 7:
                    site7[k,j]=1
                if i == 8:
                    site8[k,j]=1
                if i == 9:
                    site9[k,j]=1
                if i == 10:
                    print 'too long'
                if i in [0,1,2,3] and j > 11:
                    print 'too far index'
                if i in [4,5] and j > 14:
                    print 'too far index'
        if args.equalsites:
            tertiair = site0 + site1 + site2 + site3
            secondair = site4 + site5 + site6 + site7 + site8 + site9
            X = np.concatenate((tertiair,secondair),axis=1)
        else:
            X = np.concatenate((site0,site1,site2,site3,site4,site5,site6,site7,site8,site9),axis=1)
        return X

    def extract2DX(self,confs):
        nC= len(confs)
        nter = 12
        nsec = 15
        X = np.zeros( [ nC , nter*nsec ] )
        #bonds = ( (0,4),(2,4),(1,5),(3,5) )
        for k in range(len(confs)): # for all the configurations:
            B = np.zeros([12,15]) #so tertiary * secondairy sites
            #B = ( B04, B24, B15, B35 )
            #loop over all combinations 
            for combi in self.bonds:
                i1,i2 = combi
                group1 = confs[k][i1]
                group2 = confs[k][i2]
                if group1=='CNOO60':
                    group1='CNOO'
                if group2=='CNOO60':
                    group2='CNOO'
                gr1 = seq.index(group1)
                gr2 = seq.index(group2)
                B[ gr1, gr2 ] += 1
            if True: # so make one total matrix were all combis are combined
                #Btotal = np.sum( B , axis = 0)
                Bflatten = B.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        return X

    def test_conf(confs):
        results = []
        for conf in confs:
            for i1,i2 in zip(*self.bonds):
                if i1=='CSH' and i2=='CF':
                    results.append(conf)
        return results

############################ END CLASS ADAMANTANE

############################ START CLASS DIAMANTANE
class Diamantane(Dataset):
    '''Diamantane class'''
    ngps = (12,12,12,12,15,15) #for every instance this is same

    def __init__(self,*rgs,**kwargs):
        self.seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S','B']
        if args.symmetry:
            if False:
                self.syms = [ [ 0,1,2,3,4,5 ], [1,0,3,2,5,4] ,
                              [ 2,1,0,3,4,5 ], [1,2,3,0,5,4] ,
                              [ 0,3,2,1,4,5 ], [3,0,1,2,5,4] ,
                              [ 2,3,0,1,4,5 ], [3,2,1,0,5,4] ]
            else:
                self.syms = [ [ 0,1,2,3,4,5 ], [1,0,3,2,5,4] ]
        return

    def extract2DX(self,confs):
        nC= len(confs)
        nter = 12
        nsec = 15
        if ttert: X = np.zeros( [ nC, nter*nsec*2 ] )
        else: X = np.zeros( [ nC , nter*nsec ] )
        bonds = ( (0,4),(2,4),(1,5),(3,5) )  # bonds 0,4 and 1,5 are symmetrically similar and 2,4 and 3,5 idem dito
        
        for k in range(len(confs)): # for all the configurations:
            B04 = np.zeros([12,15]) #so tertiary * secondairy sites
            B24 = np.zeros([12,15]) #so tertiary * secondairy sites
            B15 = np.zeros([12,15]) #so tertiary * secondairy sites
            B35 = np.zeros([12,15]) #so tertiary * secondairy sites
            B = ( B04, B24, B15, B35 )
            #loop over all combinations 
            for combi,bmatrix in zip(bonds,B):
                i1,i2 = combi
                group1 = confs[k][i1]
                group2 = confs[k][i2]
                gr1 = seq.index(group1)
                gr2 = seq.index(group2)
                bmatrix[ gr1, gr2 ] = 1
            if not ttert: # so make one total matrix were all combis are combined
                Btotal = np.sum( B , axis = 0)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
            elif ttert: # keep a distinction between tert apical and tert medial sites. 
                Btertapical = np.sum( (B04,B15), axis=0)
                Btertmedial = np.sum( (B24,B35), axis=0)
                Btotal = np.concatenate( (Btertapical,Btertmedial),axis=1)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        print "Btotal constructed; shape X2:", np.shape(Btotal)
        return X

    def extractX(self,confs):
        nC= len(confs)
        #global seq
        #global newseq
        #newseq=[]
        #for item in seq:
        #    for conf in confs:
        #        if item in conf:
        #            newseq.append(item)
        #            break
        #print "newseq:", newseq
        #seq = newseq
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "seq:", seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        site0 = np.zeros([nC,12])
        site1 = np.zeros([nC,12])
        site2 = np.zeros([nC,12])
        site3 = np.zeros([nC,12])
        site4 = np.zeros([nC,15])
        site5 = np.zeros([nC,15])
        if args.verbose>1: print "dim site0", site0.shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                j = seq.index(confs[k][i]) #find the index of the group of that sequence
                if i == 0:
                    site0[k,j]=1 #index is number of configuration , number of group
                if i == 1:
                    site1[k,j]=1
                if i == 2:
                    site2[k,j]=1
                if i == 3:
                    site3[k,j]=1
                if i == 4:
                    site4[k,j]=1
                if i == 5:
                    site5[k,j]=1
                if i == 6:
                    print 'too long'
                if j > 14:
                    print 'too far index'
        X = np.concatenate((site0,site1,site2,site3,site4,site5),axis=1)
        if False: # try to make sites equal
            a1 = site0 + site1
            a2 = site2 + site3
            a3 = site4 + site5
            X = np.concatenate((a1,a2,a3),axis=1)
        return X

class Phenalene(Dataset):
    def __init__(self,nsites,*args,**kwargs):
        self.ngps = nsites * (12,)
        self.seq = [ 'N', 'CH', 'CF', 'CCFFF', 'CCHHH', 'COH', 'CNOO', 'CNHH', 'CCOOH', 'COCHHH', 'CNHCHHH' ]
        global seq
        seq = self.seq
        return

    def extractX(self,confs):
        nC= len(confs)
        nsites = len(self.ngps)
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "self.seq:", self.seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        LoS = []
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        if args.verbose>1: print "dim site0", LoS[0].shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                cleangroup = split('(\d+)',group)[0]
                j = self.seq.index(cleangroup) #find the index of the group of that sequence
                LoS[i][k,j] = 1
        X = np.concatenate(LoS,axis=1)
        return X

    def extract2DX(self,confs):
        nC= len(confs)
        X = np.zeros( [ nC , 121 ] )
        bonds = ( (0,6),(5,6),(4,8),(3,8),(1,7),(2,7),
                  (0,1),(2,3),(4,5) )

        for k in range(len(confs)): # for all the configurations:
            B = []
            for combi in bonds:
                B.append( np.zeros( [11,11] ) )
            #loop over all combinations 
            for combi,bmatrix in zip(bonds,B):
                i1,i2 = combi
                group1 = confs[k][i1]
                cleangroup1 = split('(\d+)',group1)[0]
                group2 = confs[k][i2]
                cleangroup2 = split('(\d+)',group2)[0]
                gr1 = self.seq.index(cleangroup1)
                gr2 = self.seq.index(cleangroup2)
                bmatrix[ gr1, gr2 ] = 1
                bmatrix[ gr2, gr1 ] = 1
            if True: # so make one total matrix were all combis are combined
                Btotal = np.sum( B , axis = 0)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        return X

class Propane(Dataset):
    def __init__(self, ngps=(3,3,3), seq=['CCHHH','CNHH','COH'] ):
        self.seq=seq
        self.ngps=ngps
        return

    def extractX(self, confs):
        nC = len(confs)
        LoS = [] # list of sites
        nsites = len(self.ngps)
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]

                # mv CH N O to CCHHH, CNHH, COH
                t = { 'CH': 'CCHHH', 'N':'CNHH', 'O':'COH' }
                if group in t:
                    cleangroup = t[group]
                else:
                    cleangroup = group
                j = self.seq.index(cleangroup) #find the index of the group of that sequence
                LoS[i][k,j] = 1
        X = np.concatenate(LoS,axis=1)

        print "X:", X
        print "X.shape:", X.shape
        return X

class Thiadiazinyl(Dataset):
    def __init__(self,*args,**kwargs):
        #self.seq = [ 'N', 'CH', 'CF', 'CCFFF', 'CCHHH', 'COH', 'CNOO', 'CNHH', 'CCOOH', 'COCHHH', 'CNHCHHH' ]
        #self.seq = [ 'CNHCH3', 'CSOCH3', 'COCH3', 'CSCH3', 'CSO3H', 'N', 'CCOOH', 'CCF3', 'CCH3', 'CHCO', 'CFCO', 'COOH','CSOH',
        #             'CNH2'  , 'COH', 'CSH', 'CCN', 'CH', 'CF', 'CBr', 'CCl']
        self.seq = [ 'CNHCHHH', 'CSOCHHH', 'COCHHH', 'CSCHHH', 'CSOOOH', 'N', 'CCOOH', 'CCFFF', 'CCHHH', 'CCHO', 'CCFO', 'COOH','CSOH',
                     'CNHH'  , 'COH', 'CSH', 'CCN', 'CH', 'CF', 'CBr', 'CCl', 'CSHO']
        global seq
        seq = self.seq
        self.ngps = 5 * (len(self.seq),)
        return

    def extractX(self,confs):
        print "confs[0]:", confs[0]
        nC= len(confs)
        nsites = len(self.ngps)
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "self.seq:", self.seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        LoS = []
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        if args.verbose>1: print "dim site0", LoS[0].shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                #cleangroup = split('(\d+)',group)[0]
                if group=='CSO2OH':
                    cleangroup= 'CSO3H'
                elif group=='CHNCH3':
                    cleangroup= 'CNHCH3'
                else:
                    cleangroup = group
                j = self.seq.index(cleangroup) #find the index of the group of that sequence
                LoS[i][k,j] = 1
        X = np.concatenate(LoS,axis=1)
        return X


class Pentacene(Dataset):
    def __init__(self,
            nsites=None,
            *args,**kwargs):
        self.seq = ['CH','N','CCHHH','CNHH','COH','CPh','CTh','CCl','CF','CSH','CNOO','CCCH','CCN','CNC','CCHO','CSOOOH','CCFFF']
        if not nsites is None:
            self.ngps = nsites * (len(self.seq),)
        else:
            self.ngps = None
        return

    def extractX(self,confs):
        if self.ngps is None:
            nsites = len(confs[0])
            print "nsites:", nsites
            self.ngps = nsites * (len(self.seq),)
        print "ngps:", self.ngps
        nC= len(confs)
        nsites = len(self.ngps)
        if args.verbose>=2:
            for i in range(10):
                print confs[i]
            print "self.seq:", self.seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        LoS = []
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        if args.verbose>1: print "dim site0", LoS[0].shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                j = self.seq.index(group) #find the index of the group of that sequence
                try:
                    LoS[i][k,j] = 1
                except Exception as e:
                    print "LoS[i].shape", LoS[i]
                    print "ijk", i, j, k
                    raise e
        X = np.concatenate(LoS,axis=1)
        return X

    def extract2DX(self,confs):
        seq = self.seq
        nC= len(confs)
        X = np.zeros( [ nC , len(seq)**2 ] )
        bonds = ( (0,1),(1,2),(2,3),(3,4),(4,5),(5,6))

        for k in range(len(confs)): # for all the configurations:
            B = []
            if False:
                for combi in bonds:
                    B.append( np.zeros( [len(seq),len(seq)] ) )
            else:
                from itertools import combinations
                for combi in combinations(range(7),2):
                    B.append( np.zeros( [len(seq),len(seq)] ) )
            #loop over all combinations 
            for combi,bmatrix in zip(bonds,B):
                i1,i2 = combi
                group1 = confs[k][i1]
                cleangroup1 = split('(\d+)',group1)[0]
                group2 = confs[k][i2]
                cleangroup2 = split('(\d+)',group2)[0]
                gr1 = self.seq.index(cleangroup1)
                gr2 = self.seq.index(cleangroup2)
                bmatrix[ gr1, gr2 ] = 1
                bmatrix[ gr2, gr1 ] = 1
            if True: # so make one total matrix were all combis are combined
                Btotal = np.sum( B , axis = 0)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        return X
#####################################
#####       MAIN PROGRAM       ######
#####################################

def supermain(args):
    if False:
        resultlist = []
        for column in [1,2,3]:
            args.column=column
            run=main(args)
            resultlist.append(run.coeft)
        print "resultlist:"
        print resultlist
        plots_lin(resultlist)
    else:
        main(args)
    return 

def main(args):
    identify = args.identify
    if 'ada' in identify or any( prop in identify for prop in ['GAP', 'HOMO', 'LUMO']):
        myrun = Adamantane(args.file)
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane(args.file)
    elif 'penta' in identify:
        myrun = Pentacene(name=args.file,nsites=9)
    elif 'phe' in identify:
        myrun = Phenalene(name=args.file,nsites=9)
    elif any(item in identify for item in ['fre','thia']):
        myrun = Thiadiazinyl(args.file)
    elif any(item in identify for item in ['fre','naph']):
        pass
        #myrun = Naphtol(args.file)
    else:
        print "no identify identified but I take adamantane as default!"
        myrun = Adamantane(args.file)
        #raise SystemExit('no identify was identified')
    myrun.extract(file=args.file)
    return myrun


#####################################
#####     END MAIN PROGRAM     ######
#####################################


#####################################
#####   START CALL FROM CINDES   ####
#####################################
#for default args:
class defaults(object):
   def __init__(self):
       self.verbose=2
       self.plot=0
       self.twosite=True
       self.xyplot=0
       self.symmetry=False
       self.column=1
       self.analyse=False
       self.equalsites=False
       self.fraction=None
       self.intersect=False
args = defaults()


def get_X_1D(indices, identify, descriptor='1DL',column=2, **kwargs):
    global args
    args.column=column

    if any(item in identify for item in ['GAP', 'HOMO', 'LUMO', 'ada', 'adhoma','apm6']):
        myrun = Adamantane('ada')
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane('dia')
    elif 'penta' in identify:
        myrun = Pentacene()
    elif 'pro' in identify:
        myrun = Propane(ngps=(3,3,3))
    elif 'thi' in identify:
        myrun = Thiadiazinyl()
    else:
        print "identify:", identify
        raise SystemExit('No identify_ identified')

    confs = [ indtocon(index) for index in indices ]
    if descriptor=='1DL':
        return myrun.extractX(confs=confs)
    elif descriptor=='2DL':
        return myrun.extract2DX(confs=confs)
    else:
        raise SystemExit('no valid descriptor given')


#####################################
#####    END CALL FROM CINDES   #####
#####################################

