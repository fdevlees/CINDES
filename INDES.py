#!/bin/env python 
#
#
#   THIS VERSION WAS TAKEN FROM ~/PYTHON/PHENID/BDE*/BDE_no_freq.py 13-10-2015
#
#
# import libraries
positions = (2,6,7,9,11,12) # HARD CODING positions to add a Hydrogen
from inspect import stack
import shutil #module to copy files
print stack()[0][1] # print name of python file
from platform import node
print node() # print node. so you can see from where it is submitted
import pprint # pretty printer for printing lists
import os # for getting window width and testing existence of files
import re
from re import findall # now only needed in construction.py
import sys # for getting command line input
import glob # for testing existence of files matching a pattern
import subprocess # for submitting jobs
import random # for obtaining random geometry
#import numpy as np # for using np.array although not used yet
import time # for getting time/date and time delays
now = time.strftime("%c")
print ("Current time %s" % now)
start_time = time.time()
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

# initial global variables
#(rows, columns) = os.popen('stty size', 'r').read().split() # get window width
pp = pprint.PrettyPrinter(indent=4, width=100)

#
#
semiempirical = 0
#
#

# FOR MC
kb = 8.6e-5 #boltzmann constant

#specific initial parameters.
file = "ZMAT"
siteinput = "INPUTBC"
print "name of zmatfile", file
print "name of substituent-file:", siteinput

# make a folder for all the files that are made during the calculation
workdir = os.getcwd()
path = workdir + '/databc'
print "PATH:",path
if not os.path.exists(path):
    os.makedirs(path)
    shutil.copy(os.getcwd()+'/ID_gauss',path)

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
    return data

def filemaker(core,active,passive,confs,indices,**fileparameters): #----- dict with info for filewriter has to pass here
    for i in range(len(confs)):
	c = deepcopy(core)
	a = deepcopy(active)
	p = deepcopy(passive)
	logging.debug("i=" + str(i))
	mat = zcon.constructor2(confs[i],c,a,p)
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
                
def submittingprocedure(core,active,passive,confs,indices,path,data,fileparameters):
    # here submitting thing knows at least the path
    filemaker(core,active,passive,confs,indices,**fileparameters) #----------------------------------HERE IS THE FILEWRITER CALL
    print "----- END making of the files -------------"
    #---           -------------           ---------------           --------------           ---#
    #   SYSTEM EXIT             SYSTEM EXIT               SYSTEM EXIT              SYSTEM EXIT
    #---           -------------           ---------------           --------------           ---#
    #raise SytemExit('Exit')
    # now the jobs have to be submitted 
    jobids = []
    if fileparameters['stab']==1:#then submit also the jobs in folders
        for item in indices:
            name1= item + '.com'
            jobid = subm.submit(path,name1,fileparameters['identify']).strip()
            jobids.append(jobid)
            for pos in fileparameters['positions']:
                path2 = path + '/' + item
                name2 = item + '_' + str(pos) + '.com'
                jobid = subm.submit(path2,name2).strip()
                jobids.append(jobid)
    else:
        for item in indices:
            jobid = subm.submit(path,item,fileparameters['identify']).strip()
            jobids.append(jobid)
    print "----- END all jobs are submitted ----------"
    #status = subm.jobstatus(jobid) #CANNOT WORK WITH STAB
    # test of all jobs are ready | later change to two minutes or so. 
    jobtester2(indices,jobids,path,fileparameters)
    data = datareader.datareader(indices,jobids,path,data,fileparameters)
    return data

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
    print "All jobs are READY"
    time.sleep(180) #just wait for the files to write back before opening them
    return           

def randomconf_old(subarray):
    '''this function makes a random configuration. choosing one sub for each site'''
    conf = []
    for i in range(len(array)):
        conf.append(random.choice(array[i]))
    return conf

def randomconf(subarray,maxconf,nrandsites=0): #version 4/10/2015
    '''this function makes a random configuration. choosing one sub for each site'''
    arlen = len(subarray) 
    if nrandsites == 0: #then choose a whole new configuration
        conf = []
        for i in range(arlen):
            conf.append(random.choice(array[i]))
    else: #only change nrandsites
        conf = maxconf[:] #start from same conf
        rands = random.sample(range(arlen),nrandsites) #choose nrandsites
        for i in rands:
            conf[i]=random.choice(array[i])
    return conf

def generate(p):
    '''this generates zero or one on a probability of p'''
    return random.random() <= p

def montecarloprocedure(T, subarray, max, table): #version 4/10/2015
    #MONTE CARLO PROCEDURE. 
               #maxsite = montecarloprocedure(beta, array, maximum, table)
    # INPUT: beta - maximum - table - array
    # OUTPUT: maxsite
    kb = 8.6e-5 #boltzmann constant
    beta = 1.0 / ( kb * T )
    print "monte carlo"
    cmaximum = zcon.indtocon(max[0]) # maximum is index. change to confformat
    Dtable = dict(table)
    Tcount = 0 #temperature counter. to zero after increased.
    Rcount = 0 #number of random confs tested
    while True:
        rconf = randomconf(subarray,max[0],nrandsites=2) # make a total random configuration
        rind = zcon.contoind(rconf)
        deltaetje = 0
        for i in range(len(rconf)): # now we want to have a value erandom for this configuration and test it with a certain probability
            if not rconf[i] == cmaximum[i]:
                confje = cmaximum[0:i] + [rconf[i]] + cmaximum[i+1:]
                indje= zcon.contoind(confje)
                #print "Dtable[indje]:", Dtable[indje]
                deltaetje += Dtable[indje] - max[1]
                #print "deltaetje:", deltaetje
        erandom = float (max[1] + deltaetje)
        # calculate the gradient energy. > resulttry
        p = exp(- beta * (abs( erandom - max[1] )))
        acceptance = generate(p)
        Rcount +=1
        if acceptance == 1:
            print "configuration accepted"
            print Rcount, " configurations tested"
            print "Random Conf:" , rind
            print "erandom:", erandom
            print "chance of acceptance:",p
            break
        else:
            #print "configuration not accepted"
            Tcount +=1
            if Tcount >= 10:
                Tcount = 0
                T = T * 1.1
                beta = 1.0 / ( kb * T )
    return [ rind, erandom ]

def montecarloprocedure_old(T, subarray, max, table): #version 24/09/2015
    #MONTE CARLO PROCEDURE. 
               #maxsite = montecarloprocedure(beta, array, maximum, table)
    # INPUT: beta - maximum - table - array
    # OUTPUT: maxsite
    kb = 8.6e-5 #boltzmann constant
    beta = 1.0 / ( kb * T )
    print "monte carlo"
    cmaximum = zcon.indtocon(max[0]) # maximum is index. change to confformat
    Dtable = dict(table)
    Tcount = 0
    while True:
        rconf = randomconf(subarray) # make a total random configuration
        rind = zcon.contoind(rconf)
        print "Random Conf:" , rind
        deltaetje = 0
        for i in range(len(rconf)): # now we want to have a value erandom for this configuration and test it with a certain probability
            if not rconf[i] == cmaximum[i]:
                confje = cmaximum[0:i] + [rconf[i]] + cmaximum[i+1:]
                indje= zcon.contoind(confje)
                print "Dtable[indje]:", Dtable[indje]
                deltaetje += Dtable[indje] - max[1]
                print "deltaetje:", deltaetje
        erandom = float (max[1] + deltaetje)
        print "erandom:", erandom
        # calculate the gradient energy. > resulttry
        p = exp(- beta * (abs( erandom - max[1] )))
        print "chance of acceptance:",p
        acceptance = generate(p)
        if acceptance == 1:
            print "configuration accepted"
            break
        else:
            print "configuration not accepted"
            Tcount +=1
            if Tcount >= 10:
                Tcount = 0
                T = T * 1.1
                beta = 1.0 / ( kb * T )
    return [ rind, erandom ]

# ------------------------ #
# PROGRAM MAIN STARTS HERE #
# ------------------------ #

# INPUT READING 
subinp = inr.openfile(siteinput) #this is the fileID
param = inr.readfile(subinp) #inputline is a tuple with all kind of input variables
print "PARAMETERS:", pp.pprint(param)
param["path"]=str(path)
#there are defaults given in the inputreader module

nsites = len(param['line1'])
array = inr.substireader(nsites,subinp)
logging.info("INPUT PARAMETERS:")
for key,value in param.iteritems():
    logging.info(key + ' : ' + str(value))
print "-----END INPUT READING-----"

if param['stab']==1:
    # extra parameters needed:
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
        param['gaussianline1'] =  '# opt ub3lyp/6-31g(d) pop=npa\n'
        param['gaussianline2'] =  '# geom=check guess=read b3lyp/6-311+G(d,p)\n' #also for 456
        param['gaussianline3'] =  '# geom=check guess=read b3p86/6-311+G(d,p)\n' #also for 7
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
    param['gaussianline2'] = '# geom=check guess=read ' + param['functional'] +'/'+ param['basisset'] +'\n'
else:
    param['gaussianline'] = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
    if param['twojob'] == 1: 
        param['gaussianline2'] = '# geom=check guess=read ' + param['functional'] +'/'+ param['basisset'] +'\n'

# the table with all the results of all calculated configs
if param['restart']>0:
    with open('tablebin','rb') as f:
        table = pickle.load(f)
else:
    table = []
    open('tablebin','wb').close()

# ZMAT READING 
(zmat,fileid) = r.zmatread(file)
zmatdic = r.zmatvalues(fileid)
logging.debug(pprint.pformat(zmatdic))
fileid.close()
# FORMATTING AND SPLITTING OF ZMATRIX
zmat = r.zmatprinter(zmat,zmatdic)
logging.debug("zmat:\n" + pprint.pformat(zmat))
(coremat, activematrix, passivematrix) = r.sitesplitter(zmat, param['ncore'], param['line1'], param['nch3'])
# now i save here the matrices for later use, and then the others are allowed to change for each molecule
corematrix0 = coremat
activematrix0 = activematrix
passivematrix0 = passivematrix
logging.info('coremat:' + pprint.pformat(coremat))
logging.info('activemat:' + pprint.pformat(activematrix))
logging.info('passivemat:' + pprint.pformat(passivematrix))
print "----- END FORMATTING & SPLITTING -----"

# START GEOM DEFINER
print "random start molecule: "
startconf = []
if param['restart'] >= 2:
    print "I am restarting from this configuration:"
    if 'startconf' in param:
        print "read from input file:"
        startconf = zcon.indtocon(param['startconf'])
    else:
        print "read hard coded in main file"
        startconf = [['C', 'N', 'H', 'H'], ['C', 'C', 'O', 'O', 'H'], ['C', 'C', 'N'], ['C', 'O', 'H'], ['S'], ['C', 'O']]
    print startconf
else:
    for i in range(len(array)):
        startconf.append(random.choice(array[i]))
if param['restart'] >= 3:
    maximum = startconf[:]

if param['restart'] == 4:
    sequence = param['sequence']
        # ------------------------------------- # 
        # --- HERE THE MAIN LOOP STARTS --- --- #
        # ------------------------------------- #
count = 1 # so we start counting at 1!
while True:
    print "-----------------------\n  COUNT: ", count, "\n-----------------------"
    if count==1 and param['restart']==4:
        sequence=param['sequence']
    else:
        sequence=random.sample(range(nsites),nsites)
    logging.warning(str(sequence))
    for l in range(len(sequence)):
        k = sequence[l]
        print "k=",k  , "l=",l
        if not l == 0 or count > 1:
            # define new starting geometry
            print "maxsite[0]",maxsite[0]
            startconf = zcon.indtocon(maxsite[0])
            print "newconf: ", startconf

        data=[] # here the data is stored for one cycle
        # MAKE HERE CONFIGURATIONS OF ALL SUBSTITUENTS FOR SITE k
        configurations =  [ startconf[0:k] + [array[k][i]] + startconf[k+1:] for i in range(len(array[k]))]
        indices,data,configurations = zcon.indexmaker(configurations,data,table)
        print "indices: " 
        pprint.pprint(indices)
        del startconf
        print "----- END random start configurations -----"
        # START OF SUBMITTING PART
        if not param['nosub']==1:
            data = submittingprocedure(corematrix0,activematrix0,passivematrix0,
                                         configurations,indices,path,data,param) # here call submitting procedure
        else:
            print "submit is skipped! random data is generated"
            import string
            data = skipper(indices,data)    
        #---
        with open('cyclesinfo','a') as cfid:
            filedata=deepcopy(data[:])
            for item in filedata:
                item.extend([count,k,l])
                cfid.write(' '.join(pprint.pformat(i) for i in item)+'\n')
        del filedata
        #---        

        logging.info('data:' + pprint.pformat(data))
        # maximum of the list or MINIMUM
        if 'bcprop' in param:
            #test if BC fullfilled. 
	    try:
		voldoende = [ it for it in data if int(it[2]) < int(param['bcval']) ]
	    except TypeError: pass
            print "voldoende:\n", pprint.pprint(voldoende)
            print "the BC condition is bc<:", param['bcval']
            if voldoende==[]: #so if there is at least one fullfilling BC
                bcok=0
                print "BC not fullfilled:"
                maxsite = min(data,key = lambda x:x[2])
            else: #BC nog niet
                bcok=1
                print "BC fullfilled; voldoende is not empty:", pprint.pprint(voldoende)
                maxsite = max(voldoende,key = lambda x:x[1])
        else:
            bcok=1 #no BC but need this variable to test later on
            maxsite = min(data,key = lambda x:x[1])
        logging.warning('maxisite:' + pprint.pformat(maxsite))
            
        # here move the new data to table except duplicates
        for item in data:
            if not item[0] in [tja[0] for tja in table]: table.append(item)
        with open('tablebin','wb') as tfid: # write the table to a file
            pickle.dump(table,tfid)
    # test if this is same as previous maximum. if so then converged and break
    if count > 1 and bcok: #BCOK is a test of the boundary condition is already fullfilled
        if maximum[1] == maxsite[1]: 
            print "maximum is the same!"
            print "converged to a maximum configuration!"
            if param['montecarlo'] == 0:
               break
            else:
               maxsite = montecarloprocedure(param['montecarlo'], array, maximum, table)
               print "maxsite:",maxsite
        else:
            print "maximum and maxsite are not the same yet"
            print "maximum:" ,maximum
            print "maxsite:" ,maxsite
    else: #except NameError:
        print "NameError no maximum or BC not yet fullfilled."
        #pass
    maximum = maxsite[:]
    count +=1
    print("--- %s seconds ---" % (time.time() - start_time))
    
    if count > param['maxiter']:
        print "maxiterations is reached"
        print "maximum is: ", maximum
        break
# ---------------------------- #
# ------ END OF LOOPING ------ # 
# ---------------------------- #
print "DONE"
