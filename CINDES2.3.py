#!/bin/env python 
#
#
#   THIS VERSION WAS TAKEN FROM ~/INDES/CINDES2.2.py 
#   goal of this version is to include computational reduction by prescreaning via ML
#
#
# import libraries
from writings import log_io, print_title
from inspect import stack
import shutil #module to copy files
from platform import node
import pprint # pretty printer for printing lists
import os # for getting window width and testing existence of files
import re
from re import findall # now only needed in construction.py
import sys # for getting command line input
import glob # for testing existence of files matching a pattern
import random # for obtaining random geometry
#import numpy as np # for using np.array although not used yet
import time # for getting time/date and time delays
import pickle # for saving and getting the tablebin
import logging # instead of the large amount of print statements not using it at the moment
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
from itertools import izip
from math import exp #exp(x) returns e^x
from copy import deepcopy # for keeping matrices while changing others
# import my own modules
import inputreader as inr
import construction as zcon #all functions needed for constructing new geometries
import reader as r # this reads the zmatrix in gaussian format
import submitter as subm
import datareader

class Unbuffered(object):
    def __init__(self,stream):
        self.stream = stream
    def write(self,data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self,attr):
        return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

class Run(object):
    "This is the main object for all the parameters used during the process"
    def __init__(self,**entries):
        #self.name=name
        self.script = stack()[0][1]
        self.node = node()
        self.starttime = time.time()
        self.directory = os.getcwd()
        self.pid = os.getpid()
        self.ppid = os.getppid()
        #print(time.strftime('%c'))
        self.__dict__.update(entries) #here all the key/value pairs in entries are converted to attributes.

    def __str__(self):
        sb=[]
        for key in self.__dict__:
            sb.append("{key:20}='{value}'".format(key=key, value=self.__dict__[key]))
        return '\n'.join(sb)

    def __repr__(self):
        return self.__str__()

    def currenttime(self):
        return "Current time %s" % str(time.time() - self.starttime)

    def set_calculation_properties(self):
        if self.program in ['Gaussian','gaussian']:
            self.runspecs_gaussian()
        elif self.program in ['ORCA','orca','Orca']:
            runspecs_orca(param)
        else:
            raise SystemExit('PROGRAM NOT RECOGNIZED')
        return

    def predictor(self,table,allindices,data,count):
        '''makes the predictions using KRR(ML) / RR(LS) / DIF(MC) '''
        if self.ml==1 and not table==[] and not allindices==[]:
            import learning as ml
            #preds_ml = ml.machinelearning3(allindices,table,**TZmat)
            preds_ml = ml.machinelearning2(allindices,table,**TZmat)
            #print "AllIndices & MACHINE LEARNING PREDICTIONS:"
            #for index, pred_ml in zip(allindices, preds_ml):
            #    print index, pred_ml
        elif self.ml==2 and not table==[] and not allindices==[]:
            #reduce indices | not implemented!
            pass
        else: preds_ml=[]
        if self.regression==1 and not table==[] and not allindices==[]:
            import fitter
            self.printlevel=1
            param = self.__dict__
            preds_ls = fitter.regression(table,allindices,**param)
            print "AllIndices & PREDICTIONS:"
            for index, prediction in zip(allindices, preds_ls):
                print index, prediction
        else: preds_ls=[]
        if self.difmodel==1 and not table==[] and not allindices==[] and (self.restart>2 or count>1):
            import fitter
            instance = fitter.get_instance()
            preds_dif = fitter.dif_predict(allindices,maximum,instance)
            print "AllIndices & PREDICTIONS DIFMODEL:"
            for index, prediction in zip(allindices, preds_dif):
                print index, prediction
        else: preds_dif = []
        if preds_ml==[] and preds_ls==[]:
            print "all indices: "
            pprint.pprint(allindices)
        predictions = [preds_ml, preds_ls, preds_dif ]

        if self.ml==2: #prescrean calculate only the best 50 %
            raise SystemExit('prescreaning not implemented')
            pred_data = zip(allindices, preds_ml) #get indices and predictions in same list
            pred_ml_sorted = sorted(pred_data,key=lambda x:x[1] ) #sort them based on prediction
            if self.optimum=='maximum': pred_ml_sorted = pred_ml_sorted[::-1] #when not optimum minimum reverse the list
            for index in allindices:
                if not item in data:
                    pass
                    # add worst 50 % to data_nocal
                    # add best 50 % to indices_tocal
        elif self.regression==2:
            raise SystemExit('prescreaning not implemented')
            pass
        else:
            data_nocal = data
            indices_tocal = allindices
        return predictions, data_nocal, indices_tocal

    def runspecs_gaussian(self):
        param=self.__dict__
        if param['stab']==1:
            # extra parameters needed:
            #param['positions'] = (2,6,7,9,11,12) # HARD CODING positions to add a Hydrogen
            self.gasconstant = 8.3144621
            self.bde_a = -12.68 #kJ/mol/eV^2
            self.bde_b = -218.1 #kJ/mol
            self.stab_h = 235.8 #kJ/mol
            self.Dw_h = 0.063 #eV
            self.chi_h = 2.20 
            self.chi_c = 2.60
            self.chi_n = 3.05
            self.H_h = -0.516817233 #a.u.
            self.avtc = -28.1290706 #kJ/mol #average thermal correction for 5 random structures kJ/mol
            if param['semiempirical'] == 1:
                self.gaussianline1 =  '# opt am1 \n'
                self.gaussianline2 =  '# geom=check am1\n' #also for 456
                self.gaussianline3 =  param['gausline2']
            else:
                self.gaussianline1 =  '# opt ub3lyp/6-31g(d) pop=npa \n'
                self.gaussianline2 =  '# geom=check guess=read b3lyp/6-311+G(d,p) scf=xqc\n' #also for 456
                self.gaussianline3 =  '# geom=check guess=read b3p86/6-311+G(d,p) scf=xqc\n' #also for 7
        elif param['polar']==1:
            if param['volume']==1:
                self.gaussianline = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianlinepolar = '#p geom=allcheck guess=read polar volume=tight '+param['functional']+'/'+param['basisset']+'\n'
            else:
                self.gaussianline = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianlinepolar = '#p geom=allcheck guess=read polar '+param['functional']+'/'+param['basisset']+'\n'
        elif param['ip']==1 or param['ea']==1:
            self.gaussianline = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['twojob']=1
            #param['multiplejobs']=1
            self.gaussianline2 = '# geom=check guess=read ' + param['functional'] +'/'+ param['basisset'] +'\n'
        else:
            if param['property']=='dipole':
                logging.warning('NO Geometry optimization will be performed!!!')
                self.gaussianline = '# ' + param['functional'] +'/'+ param['basisset'] +'\n'
            else:
                self.gaussianline = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            if param['twojob'] == 1: 
                self.gaussianline2 = '# geom=check guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
        for key in ['ip','ea','polar']:
            if param[key]==1:
                self.multiplejobs +=1
        return
#END CLASS RUN -----------------------------------------

def skipper(indices,data):
    ''' generate random data '''
    print "submit is skipped! random data is generated"
    import string
    for item in indices:
        propx= len(item.replace('_',''))
        if 'bcprop' in param:
            propy= sum([ string.uppercase.index(itempje)+1 for itempje in list(item.replace('_',''))]) 
            data.append([item,propx,propy])
        else:
            data.append([item,propx])
    return data

def submittingprocedure(confs,indices,data,myrun,**kwargs):
    global once
    # here submitting thing knows at least the path
    fileparameters = myrun.__dict__
    if myrun.program in ['ORCA','orca','Orca']:
        import orcafunctions
        data = orcafunctions.submittingprocedure(confs,indices,data,fileparameters,**kwargs)
    elif myrun.program in ['Gaussian','gaussian']:
        import gaussianfunctions as gausf
        data = gausf.procedure(myrun,confs,indices,data,kwargs)
    else:
        print "program not recognized!:", myrun.program
        raise SystemExit('no program recognized')
    return data

def randomconf(subarray,maxconf,nrandsites=0): #version 20/01/2016
    '''this function makes a random configuration. choosing one sub for each site'''
    arlen = len(subarray) #is length of subarray
    #print "arraylength:", arlen
    while True:
        if nrandsites == 0: #then choose a whole new configuration
            conf = []
            for i in range(arlen):
                conf.append(random.choice(array[i]))
        else: #only change nrandsites
            conf = maxconf[:] #start from same conf
            #this gives an error because range(arlen) seems to an integer.
            sitenumbers = range(arlen)
            rands = random.sample(sitenumbers,nrandsites) #choose nrandsites
            for i in rands:
                conf[i]=random.choice(subarray[i])
        if not conf==maxconf: break
    return conf

def generate(p):
    '''this generates zero or one on a probability of p'''
    return random.random() <= p

@log_io()
def montecarloprocedure(fileparameters, subarray, maxi, table,**kwargs): #version 4/10/2015
    #MONTE CARLO PROCEDURE. 
               #maxsite = montecarloprocedure(beta, array, maximum, table)
    # INPUT: beta - maximum - table - array
    # OUTPUT: maxsite
    from operator import mul
    T = fileparameters['montecarlo']
    kb = 8.6e-5 #boltzmann constant
    beta = 1.0 / ( kb * T )
    print "Monte Carlo switched on!"
    print "intial temperature is: ",T
    cmaximum = zcon.indtocon(maxi[0]) # maximum is in index format. change to confformat
    Dtable = dict([ item[:2] for item in table] )
    Tcount = 0 #temperature counter. to zero after increased.
    Rcount = 0 #number of random confs tested
    Tcountmax = int ( 10** ( float( 1 + fileparameters['nrandsites'] )/ 2 ) ) 
    print "Number of tested configurations per temperature:", Tcountmax

    #FOR ML
    if fileparameters['ml']==1:
        import learning
        ml_instance = learning.MC_init(table, **kwargs)

    while True:
        rconf = randomconf(subarray,cmaximum,fileparameters['nrandsites']) # make a total random configuration
        rind = zcon.contoind(rconf)
        deltaetje = 0
        for i in range(len(rconf)): # now we want to have a value erandom for this configuration and test it with a certain probability
            if not rconf[i] == cmaximum[i]:
                confje = cmaximum[0:i] + [rconf[i]] + cmaximum[i+1:]
                indje= zcon.contoind(confje)
                #print "Dtable[indje]:", Dtable[indje]
                deltaetje += Dtable[indje] - maxi[1]
                #print "deltaetje:", deltaetje
        erandom = float (maxi[1] + deltaetje)

        # FOR MACHINE LEARNING
        if fileparameters['ml']==1:
            print "ml_instance:", ml_instance
            print "indje:", indje
            erandom_ML = learning.MC_test_ind(ml=ml_instance, indices=[indje],**kwargs) 
            print "erandom_ML:", erandom_ML

        # calculate the gradient energy. > resulttry
        p = exp(- beta * (abs( erandom - maxi[1] )))
        acceptance = generate(p)
        Rcount +=1
        if rind == maxi[0]:
            print "rconf similar to maxconf. not accepted"
            acceptance = 0
        if acceptance == 1:
            print maxi[0],'\n',rind
            print "configuration accepted"
            print Rcount, " configurations tested"
            print "Random Conf:" , rind
            print "property_random:", erandom
            print "chance of acceptance:",p
            print "final temperature while acceptance:", T
            break
        else:
            #print "configuration not accepted"
            Tcount +=1
            if Tcount >= Tcountmax: #once in hundred configurations, the temperature is increased by 10%
                Tcount = 0
                T = T * 1.1
                beta = 1.0 / ( kb * T )
    return [ rind, erandom ]

def runspecs_orca(param):
    if param['stab']==1:
        # extra parameters needed:
        #param['positions'] = (2,6,7,9,11,12) # HARD CODING positions to add a Hydrogen
        gasconstant = 8.3144621
        bde_a = -12.68 #kJ/mol/eV^2
        bde_b = -218.1 #kJ/mol
        stab_h = 235.8 #kJ/mol
        Dw_h = 0.063 #eV
        chi_h = 2.20 
        chi_c = 2.60
        chi_n = 3.05
        H_h = -0.516817233 #a.u.
        avtc = -28.1290706 #kJ/mol #average thermal correction for 5 random structures kJ/mol
        if param['semiempirical'] == 1:
            pass
        else:
            param['orcaline1'] = '! opt'
            param['orcaline2'] = '! dft'
            param['functional1'] = 'B3LYP'
            param['functional2'] = 'B3P86'
            param['basisset1'] = '6-31g(d)'
            param['basisset2'] = '6-311+G(d,p)'
    else:
        if param['ip']==1 or param['ea']==1:
            param['orcaline'] = '! opt\n'
        else: #band gap optimization
            param['orcaline'] = '! opt\n'
        # here sum up how many extra jobs there are for dataanalysis. 
        for key in ['ip','ea','polar','IP','EA']:
            if param[key]==1:
                param['multiplejobs']+=1
    return

def geometry(file,param):
    '''reads the zmat from a file and splits it'''
    zmat,fileid = r.zmatread(file)
    zmatdic = r.zmatvalues(fileid)
    logging.debug(pprint.pformat(zmatdic))
    fileid.close()
    # FORMATTING AND SPLITTING OF ZMATRIX
    zmat = r.zmatprinter(zmat,zmatdic)
    logging.debug("zmat:\n" + pprint.pformat(zmat))
    (coremat, activemat, passivemat) = r.sitesplitter(zmat, param['ncore'], param['line1'], param['nch3'])
    # now i save here the matrices for later use, and then the others are allowed to change for each molecule
    logging.info('coremat:' + pprint.pformat(coremat))
    logging.info('activemat:' + pprint.pformat(activemat))
    logging.info('passivemat:' + pprint.pformat(passivemat))
    logging.info("----- END FORMATTING & SPLITTING -----")
    Total_Zmat = { 'core':coremat, 'active':activemat, 'passive':passivemat }
    return Total_Zmat

def testmax(myrun, data, bcok):
    param = myrun.__dict__
    if 'bcprop' in param:
        if param['bcoptimum'] in ['min','Min','MIN']:
        #test if BC fullfilled. 
            try:
                voldoende = [ it for it in data if float(it[2]) < float(param['bcval']) ]
            except TypeError: pass
        else:
            assert param['bcoptimum'] in ['max','Max','MAX']
            try:
                voldoende = [ it for it in data if float(it[2]) > float(param['bcval']) ]
            except TypeError: pass

        print "voldoende:\n", pprint.pprint(voldoende)
        print "the BC condition is bc<:", param['bcval']
        if voldoende==[]: #so if there is at least one fullfilling BC
            bcok=0
            print "BC not fullfilled:"
            #maxsite = min(data,key = lambda x:x[2])
            maxsite = min(data,key = lambda x:abs( float(x[2]) - float(param['bcval']) ) )
        else: #BC nog niet
            bcok=1
            print "BC fullfilled; voldoende is not empty:", pprint.pprint(voldoende)
            #maxsite = max(voldoende,key = lambda x:x[1])
            if param['optimum']== 'minimum':
                maxsite = min(voldoende,key = lambda x:x[1])
            else:
                maxsite = max(voldoende,key = lambda x:x[1])
    else:
        bcok=1 #no BC but need this variable to test later on
        # maximum of the list or MINIMUM
        if param['optimum']== 'minimum':
            if param['cutoff'] == 0:
                maxsite = min(data,key = lambda x:x[1])
            else:
                testdata = [ item for item in data if abs(item[1]) > param['cutoff'] ]
                maxsite = min(testdata,key = lambda x:x[1])
        else:
            maxsite = max(data,key = lambda x:x[1])
    logging.warning('maxisite:' + pprint.pformat(maxsite))
    return maxsite,bcok

def runtest(param, maximum, maxsite, count, bcok,mctable=[],**kwargs):
    converged=0
    # test if this is same as previous maximum. if so then converged and break
    if (count > 1 and bcok) or param['restart']>=3: #BCOK is a test of the boundary condition is already fullfilled
        if maximum[1] == maxsite[1]: 
            print "maximum is the same!"
            print "converged to a maximum configuration!"
            if param['montecarlo'] == 0:
               converged = 1
            else:
               if param['ml']==0:
                   maxsite = montecarloprocedure(param, array, maximum, mctable)
               else:
                   maxsite = montecarloprocedure(param, array, maximum, mctable, **kwargs)
               print "maxsite:",maxsite
        else:
            print "maximum and maxsite are not the same yet"
            print "maximum:" ,maximum
            print "maxsite:" ,maxsite
    else: #except NameError:
        print "NameError no maximum or BC not yet fullfilled."
        #pass
    maximum = maxsite[:]
    return maximum, maxsite, converged

@log_io()
def loggings(data,table,count,k,l,predictions=[]):
    def formatitem(item):
        index = '{:50s}'.format(item[0])
        datas = ' '.join(( '{:15.8f}'.format(datatje) for datatje in item[1:] ) )
        itemstring = index + datas
        return itemstring

    #--- LOGGINGS: CYCLESINFO
    with open('cyclesinfo','a') as cfid:
        filedata=deepcopy(data[:])
        for item in filedata:
            item.extend([count,k,l])
            cfid.write(' '.join(pprint.pformat(i) for i in item)+'\n')
    del filedata
    #---- LOGGINGS: TABLEBIN
    # here move the new data to table except duplicates
    for item in data:
        if not item[0] in [tja[0] for tja in table]: table.append(item)
    with open('tablebin','wb') as tfid: # write the table to a file
        pickle.dump(table,tfid)

    #---- LOGGINGS: to screen
    if predictions in ([],[[],[],[]]):
        print "data:"
        for item in data:
            print '{}'.format(formatitem(item))
    elif predictions[1]==[]:
        print "data + normal regression predictions"
        for item, pred in zip(data,predictions):
            print '{} {:10.4}'.format(formatitem(item), pred)
    elif predictions[0]==[]:
        print "data + machine learing predictions"
        for item, pred_ml in zip(data,preds_ml):
            print '{} {:10.4}'.format(formatitem(item), pred_ml)
    else:
        print "data + normal regression + machine learning"
        for item, pred, pred_ml in zip(data,predictions,preds_ml):
            print '{} {:10.4} {:10.4}'.format(formatitem(item), pred, pred_ml)

    return table



# ------------------------ #
# PROGRAM MAIN STARTS HERE #
# ------------------------ #

# initial global variables
once = 0
#(rows, columns) = os.popen('stty size', 'r').read().split() # get window width
pp = pprint.PrettyPrinter(indent=4, width=100)
kb = 8.6e-5 #boltzmann constant # FOR MC
#specific initial parameters.
file = "ZMAT"
siteinput = "INPUTBC"
logging.info("        name of zmatfile:" + file)
logging.info("name of substituent-file:" + siteinput)

@log_io()
def read_input(siteinput):
    subinp = inr.openfile(siteinput) #this is the fileID
    param = inr.readfile(subinp) #inputline is a tuple with all kind of input variables
    #print "PARAMETERS:", pp.pprint(param)
    #there are defaults given in the inputreader module
    param['nsites'] = len(param['line1'])
    if param['procedure'] in [ 'genconf' ]:
        print "Generate Configuration Procedure Active"
        array = []
    else:
        array = inr.substireader(param['nsites'],subinp)
        if not param['procedure'] in ['getrandom', 'genrandom']:
            print "ARRAY:", 
            pp.pprint(array)
    if not param['procedure'] in ['getrandom', 'genrandom']:
        logging.info("INPUT PARAMETERS:")
        for key,value in param.iteritems():
            logging.info(key + ' : ' + str(value))
        logging.info("-----END INPUT READING-----")
    return param, array

@log_io()
def setup_filesystem(param):
    if not param['nosub']==1:
        param['workdir'] = os.getcwd()
        path = param['workdir']  + '/databc'
        param["path"]=str(path)
        logging.info("PATH:"+str(path))
        if not os.path.exists(param['path']):
            os.makedirs(param['path'])
        if param['program'] in ['Gaussian','gaussian']:
            shutil.copy(os.getcwd()+'/ID_gauss',param['path'])
        elif param['program'] in ['ORCA','Orca','orca']:
            shutil.copy(os.getcwd() + '/ID_orca',param['path'])
        else:
            raise SystemExit('ERROR: No valid program specified')
    return param, path

def get_startconf(param,array):
    logging.info("random start molecule: ")
    startconf = []
    if param['restart'] >= 2:
        logging.info("I am restarting from this configuration:")
        if 'startconf' in param:
            logging.info("read from input file:")
            startconf = zcon.indtocon(param['startconf'])
        else:
            logging.info("read hard coded in main file")
            startconf = [['C', 'N', 'H', 'H'], ['C', 'C', 'O', 'O', 'H'], ['C', 'C', 'N'], ['C', 'O', 'H'], ['S'], ['C', 'O']]
    else:
        if not param['startind'] == '':
            startconf = zcon.indtocon(param['startind'])
            logging.info("read startconf from input")
        else:  
            for i in range(len(array)):
                startconf.append(random.choice(array[i]))
            logging.info("constructed random start configuration")
    logging.warning("startconf:"+pprint.pformat(startconf))
    return startconf

def genrandom(param,array):
    ''' generate random structures and print '''
    #set restart and startind to default:
    print
    for _ in xrange(param['nrandom']):
            conf = []
            for i in range(len(array)):
                conf.append( random.choice(array[i]))
            print zcon.contoind(conf)
    print
    return

def generate2(core,active,passive,converter):
    from itertools import product
    confs=[]
    for item in product(*array):
        confs.append(item)
    printindices=1
    # confs to xyz:
    with open('chemspace.xyz','a') as fid:

        for i in range(len(confs)):
            c = deepcopy(core)
            a = deepcopy(active)
            p = deepcopy(passive)
            logging.debug("i=" + str(i))
            mat = zcon.constructor2(confs[i],c,a,p)
            xyz = zmatoxyz(converter,mat)
            if printindices==1:
                fid.write('{:05d} {:4d} {:s}\n'.format(i+1,len(xyz), zcon.contoind(confs[i])))
            else:
                fid.write('{:05d} {:4d}\n'.format(i+1,len(xyz)))
            #fid.write('{:12.8f}\n'.format(Y[i])) #no property for chemspace
            for item in xyz:
                fid.write('{:3s} {:10.4f} {:10.4f} {:10.4f}\n'.format(item[0],item[1][0],item[1][1],item[1][2]))
            fid.write('\n')
            # now for each mat transform to xyz
    return confs

def get_sequence(count, myrun):
    # START set sequence INPUT: param
    param = myrun.__dict__
    if 'sequences' in param:
        try:
            sequence = param['sequences'][count-1] #accounting for the fact count starts counting at 1
            return sequence
        except IndexError:
            sequence=random.sample(range(param['nsites']),param['nsites'])
            return sequence
    if count==1 and param['restart']==4:
        sequence=param['sequence']
    else:
        if param['norandom']==1:
            sequence=range(param['nsites'])
        elif param['sequence'] == []:
            sequence=random.sample(range(param['nsites']),param['nsites'])
        else:
            print "sequence read from file"
            sequence=param['sequence']
    ##output sequence
    logging.warning(str(sequence))
    return sequence

def generate_procedure(param,array):
    '''generate all structures and print in format'''
    from converter import Converter
    converter = Converter()

    #get structure
    TZmat = geometry(file,param)
    #to get an xyz file with the data from tablebin do generate1()
    import learning
    with open('tablebin','rb') as f:
        table = pickle.load(f)
    import writings
    print table[508:510]
    learning.generate1(converter=converter,table=table,**TZmat)

    #to get an xyz file with all the possible structures possible:
    #generate2(core,active,passive,converter)

    #we have to generate all possible iterations from the array

    #get table
    return

def genconf(param):
    myrun = Run(**param)
    myrun.set_calculation_properties()
    #if param['program'] in ['Gaussian','gaussian']:
    #    runspecs_gaussian(param)
    #elif param['program'] in ['ORCA','orca','Orca']:
    #    runspecs_orca(param)
    param = myrun.__dict__
    TZmat = geometry(file,param)
    conf = zcon.indtocon(param['startind'])
    print "in GENCONF: conf is:", conf
    param['workdir'] = os.getcwd()
    path = param['workdir']
    param["path"]=str(path)

    c = deepcopy(TZmat['core'])
    a = deepcopy(TZmat['active'])
    p = deepcopy(TZmat['passive'])
    mat = zcon.constructor2(conf,c,a,p)
    zcon.filewriter2(mat,param['startind'],**param)
    return

def standard_procedure(param,array):
    param, path = setup_filesystem(param)
    startconf = get_startconf(param,array)
    main(param,array,startconf)
    return

def set_table(myrun):
    if myrun.restart>0:
        with open('tablebin','rb') as f:
            table = pickle.load(f)
    else:
        table = []
        open('tablebin','wb').close()
    return table

def set_maximum(myrun,table):
    if myrun.restart >= 3:
        tabledict = dict( [ item[0:2] for item in table ] )
        maximum = [ myrun.startconf, tabledict[myrun.startconf] ]
        logging.info("maximum:"+ str(maximum))
    else:
        maximum = 0
    return maximum

def main(param,array,startconf):
    bcok=0 #TO REMOVE LATER
    param['bcok']=0
    myrun = Run(**param)
    print(myrun) #this should print all the class elements via the __str__ function

    # the table with all the results of all calculated configs
    table = set_table(myrun)

    #zmatrix reading and splitting
    TZmat = geometry(file,param)

    # set maximum
    maximum = set_maximum(myrun,table)

    #set calculation properties
    myrun.set_calculation_properties()
            # ------------------------------------- # 
            # --- HERE THE MAIN LOOP STARTS --- --- #
            # ------------------------------------- #
    count = 1 # so we start counting at 1!
    while True:
        #logging.info("-----------------------\n  COUNT: "+ str(count) + "\n-----------------------")
        print_title("COUNT: " + str(count),outline='l',signator="-")

        ### set site order in sequence INPUT: param, count
        sequence = get_sequence(count, myrun)

        ### for each site in sequence:
        for l in range(len(sequence)):
            k = sequence[l]
            print_title("k(site)= " + str(k) + " l(nsite)= "+ str(l),outline='l',signator='=')
            #print "k=",k  , "l=",l
            if not l == 0 or count > 1: #define new startconfiguration if not first cycle
                # define new starting geometry
                print "maxsite[0]",maxsite[0]
                del startconf
                startconf = zcon.indtocon(maxsite[0])
                print "newconf: ", startconf

            #get indicesall and the indices that still need to be calculated
            # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
            indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table)
            print "----- END random start configurations -----"
            print "indices_todo:",indices_todo, "indices_all:", indices_all
            print "data_nodo:", data_nodo #all item[1]==1 in data_nodo 

            # perform prescreaning in a predictions. 
            predictions,data_nocal,indices_tocal = myrun.predictor(table, indices_todo,data_nodo, count)

            # START OF SUBMITTING PART
            if not myrun.nosub==1:
                data_all = submittingprocedure(configurations,indices_tocal,
                                           data_nocal,
                                           myrun,
                                           **TZmat) # here call submitting procedure
            else: data_all = skipper(indices,data)

            # sort data in same order as allindices:
            data_all = sorted(data_all, key=lambda x:indices_all.index(x[0]))
            # logs new elements in data to table and tablebin and whole data to cyclesinfo
            table = loggings(data_all,table,count,k,l,predictions)

            # decide what the maximum site is and if the bc if fullfilled
            print "BCOK:", bcok
            maxsite, bcok = testmax(myrun, data_all, bcok)

            print("--- %s seconds ---" % (time.time() - myrun.starttime))
            print(myrun.currenttime())
        # END LOOP OVER SITES

        #get maximum and test convergence
        maximum, maxsite,converged = runtest(param, maximum, maxsite, count, bcok, mctable=table,**TZmat)
        if converged==1: break
        count +=1
        if count > param['maxiter']:
            print "maxiterations is reached"
            print "maximum is: ", maximum
            break
    # ---------------------------- #
    # ------ END OF LOOPING ------ # 
    # ---------------------------- #
    print "DONE"
def testrun(param,array):
    pass
if __name__ == "__main__":
    print_title("C I N D E S\nAn Inverse Molecular Design Program\nwritten by Jos L. Teunissen", newlines=True)
    # INPUT READING
    param, array = read_input(siteinput)
    # END INPUT READING

    if param['procedure'] == 'standard':
        standard_procedure(param,array)
    elif param['procedure'] == 'test':
        testrun(param,array)
    elif param['procedure'] == 'generate':
        generate_procedure(param,array)
    elif param['procedure'] == 'genconf':
        genconf(param)
    elif param['procedure'] in [ 'getrandom' ,'genrandom']:
        genrandom(param,array)
    else:
        logging.warning('proceduretype not recognized')
    print "bla"
