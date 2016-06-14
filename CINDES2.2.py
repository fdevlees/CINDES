#!/bin/env python 
#
#
#   THIS VERSION WAS TAKEN FROM ~/PYTHON/PHENID/BDE*/BDE_no_freq.py 13-10-2015
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
#import subprocess # for submitting jobs
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
#import mailer
 
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
    def __init__(self,name):
        self.name=name
        self.script = stack()[0][1]
        self.node = node()
        self.starttime = time.time()
        self.directory = os.getcwd()
        self.pid = os.getpid()
        self.ppid = os.getppid()
        print(time.strftime('%c'))
    def currenttime(self):
        return "Current time %s" % time.time() - self.starttime()


#
# most important FUNCTIONS 
# other functions are located in modules
#
def skipper(indices,data):
    ''' generate random data '''
    for item in indices:
        propx= len(item.replace('_',''))
        if 'bcprop' in param:
            propy= sum([ string.uppercase.index(itempje)+1 for itempje in list(item.replace('_',''))]) 
            data.append([item,propx,propy])
        else:
            data.append([item,propx])
    return data

@log_io()
def filemaker(confs,indices,fileparameters,passive,active,core,**kwargs): #----- dict with info for filewriter has to pass here
    path = fileparameters['path']
    for i in range(len(confs)):
	c = deepcopy(core)
	a = deepcopy(active)
	p = deepcopy(passive)
	logging.debug("i=" + str(i))
	mat = zcon.constructor2(confs[i],c,a,p)
        #print 
        #print mat
        logging.debug("in filemaker kwargs:")
        logging.debug(pprint.pformat(fileparameters))
        #now use this variable to test if we have to make all kind of extra files or not
         

        if not fileparameters['stab']==1:
	    zcon.filewriter2(mat,indices[i],**fileparameters) #------------------------------------------------HERE IS THE FILEWRITER CALL
        else:
            zcon.filewriterA(mat,indices[i],**fileparameters) #here we have to use makers to construct the AH files
            #we make a folder with the indexname in /data/indeces[i]
            if not os.path.exists(path + '/' + indices[i]): #path is $WORKDIR/data
                os.makedirs(path + '/' + indices[i])
                # and make sure ID_gauss is in the folder!
                shutil.copy(path +'/ID_gauss',path+'/'+indices[i])
            filename = fileparameters['path'] + '/' + fileparameters['identify'] + str(indices[i]) + ".com" #same line as in filewriter. open it again.
# ----part of changer
            logging.debug("filename: " + filename)
            multcharge = re.compile('^\-?[0-9]\s[0-9]') #only set the compiler
            with open(filename) as fid: #again open as fid
                for line in fid:
                    if multcharge.match(line): #from where there is a match it reads the subsequant lines as the zmat
                        #print "match!"
                        #print multcharge.match(line).group()
                        zmat=[]
                        line=next(fid)
                        while not line == '\n': #until empty line
                            zmat.append(line.split())
                            line=next(fid)
                        break #so that only the first match is used. after the other matches there is no zmat
            for pos in fileparameters['positions']:
                zmat2 = deepcopy(zmat)
                hornot = maker1(zmat2,pos,indices[i],**fileparameters) #returns a value indicating if there is already a hydrogen (or a nitrogen)
                # FOR NOW ONLY DO ONE POSSIBILITY THIS IS EASIER BECAUSE WE KNOW EXACTLY HOW MANY JOBS THERE HAVE TO BE SUBMITTED
                #if not hornot == 1: #if not there are two ways to place the hydrogen.
                    #maker2(zmat,pos,indices[i],**fileparameters)
            print "DONE"
    return
# ---- end part of changer
# --- maker here. old one in construction.py

def maker1(zmat,pos,index,**fileparameters):
    """makes new file with hydrogen attached on first dihedral"""
    zmatnew = deepcopy(zmat)
    spos = str(pos) #spos is string of pos. pos = position
    h=0
    #print "zmat[spos-1]:",zmat[pos-1]
    #print "pos:",spos
    #print "fileparamters ncore:", fileparameters['ncore']
    logging.debug(pprint.pformat(zmat))
    #print "zmat[fileparameters['ncore']:]:"
    #pp.pprint(zmat[fileparameters['ncore']:])
    if zmat[pos-1][0] == 'N': 
        item = zmat[pos-1]
        h=1
        if len(item)==1: #when pos is 1 so first index of a zmat
            hline = ['H',1,0.9,2,109.5,3,126.0]
        elif len(item)==3: #when pos is 2 so second index of a zmat
            hline = ['H',2,0.9,3,109.5,4,126.0]
        else:
            bondindex=item[1] #or if item only has length 1
            dihedralindex=item[3]
            hline= ['H',spos,0.9,bondindex,109.5,item[3],126.0]#LOOK AT THIS
    else:
        for item in zmatnew[fileparameters['ncore']:]:
            #print "in loop", "spos:",spos,"str(item[1]",str(item[1])
            if str(item[1]) == spos:
                if item[0] == 'H': h=1
                hline=item[:] 
                item[6] = '126.0'
                hline[6] = '234.0'
                hline[0] = 'H'
                #print "hline:",hline
    zmatnew.append(hline)      
    zcon.filewriterAH(zmatnew,spos,index,**fileparameters) #now it is important where this will be written.
    return h
# ----- end of maker
def maker2(zmat,pos,index,**fileparameters):
    '''makes new file with hydrogen attached on second dihedral. 
    this is only necessary when there is not already another hydrogen on the compound
    or that the site is nitrogen or possibly sulfur doped. '''
    zmatnew = zmat[:]
    spos = str(pos)
    h=0
    #print "pos:",spos
    for item in zmatnew[fileparameters['ncore']:]:
        if str(item[1]) == spos:
            hline=item[:] 
            item[6] = '234.0'
            hline[6] = '126.0'
            hline[0] = 'H'
    zmatnew.append(hline)      
    zcon.filewriterAH(zmatnew,spos + '_2',index,**fileparameters)
    return
# ----- end of maker2
                
def submittingprocedure(confs,indices,path,data,fileparameters,**kwargs):
    global once
    # here submitting thing knows at least the path
    if fileparameters['program'] in ['ORCA','orca','Orca']:
        import orcafunctions
        data = orcafunctions.submittingprocedure(confs,indices,data,fileparameters,**kwargs)
        #---           -------------           ---------------           --------------           ---#
        #   SYSTEM EXIT             SYSTEM EXIT               SYSTEM EXIT              SYSTEM EXIT
        #---           -------------           ---------------           --------------           ---#
        #raise SytemExit('Exit')
    else:
        if fileparameters['no1sub']==1 and once==0:
            once = 1
            print " "
        else:
            #pp.pprint(kwargs)
            
            filemaker(confs,indices,fileparameters,**kwargs) #----------------------------------HERE IS THE FILEWRITER CALL
            print "----- END making of the files -------------"
        # now the jobs have to be submitted 
        jobids = []
        indicesall = deepcopy(indices) #here i copy the indices. The indices are submitted. The indicesall are not all submitted but are all read out.

        if fileparameters['stab']==1:#then submit also the jobs in folders
            if fileparameters['try_ready']==1:
                indices,paths = jobtester3(indices,path,fileparameters,returnpath=True)
                print "try_ready activated"
                print "indices:", indices
                #for path in paths:
                #    print path
                #raise SystemExit('stopped')
            for item in indices:
                name1 = item + '.com'
                jobid = subm.submit(path,name1,fileparameters['identify']).strip()
                jobids.append(jobid)
                for pos in fileparameters['positions']:
                    path2 = path + '/' + item
                    name2 = item + '_' + str(pos) + '.com'
                    jobid = subm.submit(path2,name2,fileparameters['identify']).strip()
                    jobids.append(jobid)
        elif once==1 and fileparameters['no1sub']==1:
            print "submit skipped"
            once = 2
        else:
            if fileparameters['try_ready']==1:
                indices = jobtester3(indicesall,path,fileparameters)
            for item in indices:
                name = item + '.com'
                if fileparameters['nosub']==2:
                    time.sleep(1)
                    jobid = subm.nosubmit(path,item,fileparameters['identify'])
                    print item + 'submitted'
                else:
                    jobid = subm.submit(path,name,fileparameters['identify']).strip()
                jobids.append(jobid)
        logging.info("----- END all jobs are submitted ----------")
        #status = subm.jobstatus(jobid) #CANNOT WORK WITH STAB
        # test of all jobs are ready | later change to two minutes or so. 
        jobtester2(indices,jobids,path,fileparameters)
        print "indicesall:", indicesall
        data = datareader.datareader(indicesall,jobids,path,data,fileparameters)
    return data

def jobtester3(indices,path,fileparameters,returnpath=False):
    paths = [] #here we are going to make a list of paths of the jobs
    if 'positions' in fileparameters: positions = fileparameters['positions']
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'][:-1] + '*_' + indices[i] + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
                paths.append(path2)
    indicescopy = deepcopy(indices)
    #print "paths:", paths
    newpaths = paths[:]
    if fileparameters['stab']==1: #test if all necessary A and AH calculations are performed
        k=0
        for i in range(len(indices)): # all indices
            l=0
            if glob.glob(paths[k]): # test A
                print "already calculated:", paths[k]
                newpaths.remove(paths[k])
                l+=1
            for j in range(len(positions)): # test all AH
                k +=1
                if glob.glob(paths[k]):
                    print "already calculated:", paths[k]
                    newpaths.remove(paths[k])
                    l+=1
            print "len(positions):", len(positions)
            print "l:", l
            if l == len(positions) + 1: #if all AH and A then remove from indices
                indices.remove(indicescopy[i])
            k+=1
    else:
        for i in range(len(paths)):
            if glob.glob(paths[i]):
                print "already calculated:", indicescopy[i]
                indices.remove(indicescopy[i])
    if returnpath:
        return indices,newpaths
    else:
        return indices

@log_io(signator='=')
def jobtester2(indices,jobids,path,fileparameters):
    tijdje = 0
    paths = [] #here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
                paths.append(path2)
    while True: # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje>fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy= paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths==[]:
            break
        print "time/h:", tijdje/3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje+=fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime']) #just wait for the files to write back before opening them
    return           

def randomconf_old(subarray):
    '''this function makes a random configuration. choosing one sub for each site'''
    conf = []
    for i in range(len(array)):
        conf.append(random.choice(array[i]))
    return conf

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

def runspecs(param):
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
            param['gaussianline1'] =  '# opt am1 \n'
            param['gaussianline2'] =  '# geom=check am1\n' #also for 456
            param['gaussianline3'] =  param['gausline2']
        else:
            param['gaussianline1'] =  '# opt ub3lyp/6-31g(d) pop=npa \n'
            param['gaussianline2'] =  '# geom=check guess=read b3lyp/6-311+G(d,p) scf=xqc\n' #also for 456
            param['gaussianline3'] =  '# geom=check guess=read b3p86/6-311+G(d,p) scf=xqc\n' #also for 7
    elif param['polar']==1:
        if param['volume']==1:
            param['gaussianline'] = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['gaussianlinepolar'] = '#p geom=allcheck guess=read polar volume=tight '+param['functional']+'/'+param['basisset']+'\n'
        else:
            param['gaussianline'] = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['gaussianlinepolar'] = '#p geom=allcheck guess=read polar '+param['functional']+'/'+param['basisset']+'\n'
    elif param['ip']==1 or param['ea']==1:
        param['gaussianline'] = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
        param['twojob']=1
        #param['multiplejobs']=1
        param['gaussianline2'] = '# geom=check guess=read ' + param['functional'] +'/'+ param['basisset'] +'\n'
    else:
        if param['property']=='dipole':
            logging.warning('NO Geometry optimization will be performed!!!')
            param['gaussianline'] = '# ' + param['functional'] +'/'+ param['basisset'] +'\n'
        else:
            param['gaussianline'] = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
        if param['twojob'] == 1: 
            param['gaussianline2'] = '# geom=check guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
    
    for key in ['ip','ea','polar']:
        if param[key]==1:
            param['multiplejobs']+=1
    return

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


def testmax(param, data, bcok):
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

#        try:
#    	voldoende = [ it for it in data if int(it[2]) < int(param['bcval']) ]
#        except TypeError: pass
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
def loggings(data,table,count,k,l,predictions=[],preds_ml=[]):
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
    if predictions==[] and preds_ml==[]:
        print "data:"
        for item in data:
            print '{}'.format(formatitem(item))
    elif preds_ml==[]:
        print "data + normal regression predictions"
        for item, pred in zip(data,predictions):
            print '{} {:10.4}'.format(formatitem(item), pred)
    elif predictions==[]:
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
# FOR MC
kb = 8.6e-5 #boltzmann constant
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
 
def get_sequence(count, param):
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
    if param['program'] in ['Gaussian','gaussian']:
        runspecs(param)
    elif param['program'] in ['ORCA','orca','Orca']:
        runspecs_orca(param)

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

def main(param,array,startconf):
    bcok=0
    maximum = 0 
    myrun = Run('myrun')
    logging.info(pprint.pformat('starttime' + str(myrun.starttime)))
    logging.info(pprint.pformat('script' + myrun.script))
    logging.info(myrun.directory)
    logging.info(myrun.node)
    logging.info(pprint.pformat('pid: ' + str(myrun.pid)))
    logging.info(pprint.pformat('ppid: ' + str(myrun.ppid)))

    # the table with all the results of all calculated configs
    if param['restart']>0:
        with open('tablebin','rb') as f:
            table = pickle.load(f)
    else:
        table = []
        open('tablebin','wb').close()
    
    #zmatrix reading and splitting
    TZmat = geometry(file,param)
    
    # START set maximum INPUT: param, table
    if param['restart'] >= 3:
        tabledict = dict( [ item[0:2] for item in table ] )
        maximum = [ param['startconf'] , tabledict[param['startconf'] ] ]
        logging.info("maximum:"+ str(maximum))
    # END set maximum OUTPUT maximum
    
    # START set sequence INPUT: param
    if param['restart'] == 4:
        sequence = param['sequence']
    # END set sequence OUTPUT sequence
    
    if param['program'] in ['Gaussian','gaussian']:
        runspecs(param)
    elif param['program'] in ['ORCA','orca','Orca']:
        runspecs_orca(param)
            # ------------------------------------- # 
            # --- HERE THE MAIN LOOP STARTS --- --- #
            # ------------------------------------- #
    count = 1 # so we start counting at 1!
    while True:
        #logging.info("-----------------------\n  COUNT: "+ str(count) + "\n-----------------------")
        print_title("COUNT: " + str(count),outline='l',signator="-")
        ### set site order in sequence INPUT: param, count
        sequence = get_sequence(count, param)
        ##output sequence
        logging.warning(str(sequence))
    
        ### for each site in sequence:
        for l in range(len(sequence)):
            k = sequence[l]
            print_title("k(site)= " + str(k) + " l(nsite)= "+ str(l),outline='l',signator='=')
            #print "k=",k  , "l=",l
            if not l == 0 or count > 1:
                # define new starting geometry
                print "maxsite[0]",maxsite[0]
                startconf = zcon.indtocon(maxsite[0])
                print "newconf: ", startconf
    
            data=[] # here the data is stored for one cycle
    
            # MAKE HERE CONFIGURATIONS OF ALL SUBSTITUENTS FOR SITE k
            configurations =  [ startconf[0:k] + [array[k][i]] + startconf[k+1:] for i in range(len(array[k]))]
    
            logging.debug(pprint.pformat(configurations))
    
            indices,data,configurations,allindices = zcon.indexmaker(configurations,data,table)

            print "indices:",indices, "allindices:", allindices
            if param['ml']==1 and not table==[] and not allindices==[]:
                import learning as ml
                #preds_ml = ml.machinelearning3(allindices,table,**TZmat)            
                preds_ml = ml.machinelearning2(allindices,table,**TZmat)            
                #print "AllIndices & MACHINE LEARNING PREDICTIONS:"
                #for index, pred_ml in zip(allindices, preds_ml):
                #    print index, pred_ml
            elif param['ml']==2 and not table==[] and not allindices==[]:
                #reduce indices | not implemented!
                pass
            else: preds_ml=[]
            if param['regression']==1 and not table==[] and not allindices==[]:
                import fitter
                param['printlevel']=1
                predictions = fitter.regression(table,allindices,**param)
                print "AllIndices & PREDICTIONS:"
                for index, prediction in zip(allindices, predictions):
                    print index, prediction
            else: predictions=[]
            if param['difmodel']==1 and not table==[] and not allindices==[] and (param['restart']>2 or count>1):
                import fitter
                instance = fitter.get_instance()
                preds_dif = fitter.dif_predict(allindices,maximum,instance)
                print "AllIndices & PREDICTIONS DIFMODEL:"
                for index, prediction in zip(allindices, preds_dif):
                    print index, prediction
            else: preds_dif = []
            if preds_ml==[] and predictions==[]:    
                print "all indices: " 
                pprint.pprint(allindices)
            del startconf
            print "----- END random start configurations -----"
    
            # START OF SUBMITTING PART
            if not param['nosub']==1:
                data = submittingprocedure(configurations,indices,
                                           param['path'],
                                           data,
                                           param, **TZmat) # here call submitting procedure
            else:
                print "submit is skipped! random data is generated"
                import string
                data = skipper(indices,data)
            # sort data in same order as allindices:
            data = sorted(data, key=lambda x:allindices.index(x[0]))
            # logs new elements in data to table and tablebin and whole data to cyclesinfo
            table = loggings(data,table,count,k,l,predictions,preds_ml)
    
            # decide what the maximum site is and if the bc if fullfilled
            print "BCOK:", bcok
            maxsite, bcok = testmax(param, data, bcok)
            
        print("--- %s seconds ---" % (time.time() - myrun.starttime))
        # END LOOP OVER SITES
        
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
   
    #----
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
