# inputreader module

import logging
from pprint import pprint
import re
inrlog = logging.getLogger('substireader')
inrlog.setLevel(logging.INFO)
# set handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
# set formatter
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)

inrlog.addHandler(ch)

def openfile(filename):
    '''opens file in reading mode and returns fileid'''
    subinp = open(filename,'r')   
    return subinp

def readfile(subinp):
    '''this method reads all the inputkeywords'''
    #default values
    paras={'program':'gaussian',
           'procedure':'standard',
           'optimum':'minimum',
           'ml':0,
           'nosub':0, 
           'cutoff':0,
           'no1sub':0,
           'try_ready':0,
           'regression':0,
           'difmodel':0,
           'nprocs':2,
           'debug':False,
           'norandom':0,
           'sequence':[],
           'restart':0,
           'semiempirical':0,
           'twojob':0,
           'multiplejobs':0,
           'stab':0,
           'semiempirical':0,
           'polar':0,
           'ip':0,
           'ea':0,
           'startind': '',
           'extrawaittime': 2,
           'timelimit':250000,
           'timestep':300,
           'maxiter':10,
           'basisset':'6-31G',
           'functional':'b3lyp',
           'ncore':10,
           'nch3':16,
           'identify':'unspecified_',
           'property':'gap',
           'charge':0,
           'mult':1,
           'montecarlo':0,   #Temperature at start
           'nrandsites':2}   #n random sites changed. for all choose 0
    #scans all the lines until if will find the END keyword
    #this is a bit tricky because keywords can appear everywere in the file before END 
    #but the value of the keywords is always on the second and further position
    # i can if i want change all lines with:
    # if 'mystring' in line:
    # substitute by
    # if 'mystr' in line.split()[0]:
    while True:
        line = subinp.readline()
        if not line: break
        if 'bc' in line:
            paras['bcprop'] = line.split()[1]
            paras['bcval'] = line.split()[2] 
            try:
                paras['bcoptimum'] = line.split()[3]
            except IndexError:
                paras['bcoptimum'] = 'min'
            assert paras['bcoptimum'] in ['min','Min','max','Max','MIN','MAX']
        if 'nosub' in line: 
            try:
                paras['nosub'] = int( line.split()[1] )
            except IndexError:
                paras['nosub'] = 1
        if 'restart' in line: 
            paras['restart'] = int(line.split()[1])
            if paras['restart'] > 1:
                try:
                    paras['startconf'] = line.split()[2]
                except IndexError:
                    logging.warning("no startconf given while expected!")
            if paras['restart'] > 3: #sequence to do in first run
                line = subinp.readline()
                # line=next(subinp)
                paras['sequence']= [ int(item) for item in line.split() ]
        #example:
        # restart 4 CH_COH_CCHHH_N_CCOOH
        # 3 4 5
        if 'ml'==line[:2]:
            paras['ml']=1

        if 'semiempirical' in line: paras['semiempirical'] = 1
        if 'nprocs' in line: paras['nprocs'] = int(line.split()[1])
        if 'debug' in line: paras['debug'] = True
        if 'norandom' in line: paras['norandom'] = 1
        if 'no1sub' in line: paras['no1sub'] = 1
        if 'regression' in line: paras['regression'] = 1
        if 'difmodel' in line: paras['difmodel'] = 1
        if 'optimum' in line: paras['optimum'] = line.split()[1]
        if 'sequence' in line: 
            nsequences = int(line.split()[1])
            sequences = []
            for _ in range(nsequences):
                line = subinp.readline()
                sequence  = [ int(item) for item in line.split() ]
                sequences.append(sequence)
            paras['sequence'] = sequences[0]
            paras['sequences'] = sequences
        if 'simple' in line: paras['simple'] = 1
        if 'try_ready' in line: paras['try_ready'] = 1
        if 'twojob' in line: paras['twojob'] = 1
        if 'startind' in line: paras['startind'] = line.split()[1]
        if 'procedure' in line: 
                paras['procedure'] = line.split()[1]
                if paras['procedure'] in ['genrandom', 'getrandom']:
                    try:
                        paras['nrandom'] = int(line.split()[2])
                    except IndexError:
                        raise SystemExit("NO number of random structures specified!")
        if 'timelimit' in line: paras['timelimit'] = int(line.split()[1])
        if 'timestep' in line: paras['timestep'] = int(line.split()[1])
        if 'maxiter' in line: paras['maxiter'] = int(line.split()[1])
        if 'mult' in line: paras['mult'] = int(line.split()[1])
        if 'charge' in line: paras['charge'] = int(line.split()[1])
        if 'basisset' in line: paras['basisset'] = line.split()[1]
        if 'functional' in line: paras['functional'] = line.split()[1]
        if 'identify' in line: paras['identify'] = line.split()[1]
        if any(item in line for item in ('ncore','natomscore')): paras['ncore'] = int(line.split()[1])
        if 'nch3' in line: paras['nch3'] = int(line.split()[1])
        if 'montecarlo' in line: 
            paras['montecarlo'] = float(line.split()[1])
            try:
                paras['nrandsites']= int(line.split()[2])
            except IndexError: pass
        if 'extrawaittime' in line: paras['extrawaittime'] = float(line.split()[1])
        if 'cutoff' in line: paras['cutoff'] = float(line.split()[1])
        if 'property' in line: paras['property'] = line.split()[1]
        if 'sites' in line: paras['line1'] = [ int(item) for item in line.split()[1:] ]       
        if 'positions' in line: paras['positions'] = [ int(item) for item in line.split()[1:] ]       
        if any(item in line.split()[0] for item in ('program','AI','Program','prog','Prog','programma')): paras['program']= line.split()[1]
        if 'END' in line: break
    # it turns out to be helpful to have a flag to know if the stab or polar property has to be calculated so:
    # NOTE THAT HERE it is not possible to use polar and stab simultaneously
    if 'bcprop' in paras: #test if we use a BC
        if paras['bcprop']=='stab': 
            paras['stab']=1 #test if maybe that one is stab. if such stab=1
        elif paras['property']=='stab':
            paras['stab']=1
        elif paras['bcprop']=='polar':
            paras['polar']=1
        elif paras['property']=='polar':
            paras['polar']=1
        if paras['bcprop'] in ['ip','IP'] or paras['property'] in ['ip','IP']: paras['ip']=1
        if paras['bcprop'] in ['ea','EA'] or paras['property'] in ['ea','EA']: paras['ea']=1
    else:
        if paras['property']=='stab':
            paras['stab']=1
        elif paras['property']=='polar':
            paras['polar']=1
        elif paras['property'] in ['ip','IP']: paras['ip']=1
        elif paras['property'] in ['ea','EA']: paras['ea']=1
    #print "in inputreader paras:", paras
    return paras

def substireader(nsit,subinp):
    '''this one reads all the different substituents for all different positions'''
    substiarray = []
    #print "nsit,subinp", nsit, subinp
    for i in range(nsit):
        # read number of substituents
        try:
            line2 = subinp.readline().split()
            nsubsit = int(line2[1])
        except IndexError:
            print "no functional groups present or wrong formatted"
            break
        inrlog.debug("site number: " + str(i+1))
        inrlog.debug("nsubsit: " + str(nsubsit))
        #print "site number: ", i+1 
        #print "nsubsit: ", nsubsit
        site = [ subinp.readline().split() for line in range(nsubsit)]
        # for j in range(nsubsit):
        #    atoms = subinp.readline().split()
        substiarray.append(site)
        # read each substituent in an array or so
    return(substiarray)

if __name__ == "__main__":
    import sys
    fid=openfile(sys.argv[1])
    line1,line2 = readfile(fid)
    print "inputline: ", line1
    print "tweede inputline: ", line2
    nsit=len(line1[-1])
    subs= substireader(nsit,fid)
    print "subs:",
    pprint(subs)
