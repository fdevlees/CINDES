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

from CINDES4.utils.writings import log_io

# MAIN FUNCTION
@log_io()
def read_input(siteinput):
    subinp = openfile(siteinput) #this is the fileID
    param = readfile(subinp) #inputline is a tuple with all kind of input variables
    #here for a new link feature. nsites is len(line1) - nlinks
    if 'nlinks' in param:
        param['nsites'] = len(param['line1']) - param['nlinks']
    else:
        param['nsites'] = len(param['line1'])
    #####################
    if param['procedure'] in [ 'genconf' ]:
        print "Generate Configuration Procedure Active"
        array = []
    else:
        array = substireader(param['nsites'],subinp)
        if not param['procedure'] in ['getrandom', 'genrandom']:
            print "ARRAY:",
            pprint(array)
    if not param['procedure'] in ['getrandom', 'genrandom']:
        logging.info("INPUT PARAMETERS:")
        for key,value in param.iteritems():
            logging.info(key + ' : ' + str(value))
    return param, array

def openfile(filename):
    '''opens file in reading mode and returns fileid'''
    subinp = open(filename,'r')
    return subinp

def get_preds(subinp, line):
    ''' for future development a more extensible format for giving which predictions are tried
    it returns a list of dictionaries with each dictionary having one obligatory type key '''
    # default predicition types: 
    defaults = { 'ML' : { 'type': 'ML', 'descriptor':'coulomb'},
                 'iML': { 'type':'iML' },
                 'NN' : { 'type': 'NN', 'descriptor':'BoB'},
                 '1D' : { 'type': '1D' },
                 '2D' : { 'type': '2D' }
               }

    npredictions = int(line.split()[1])
    preds = []
    for _ in range(npredictions):
        line = subinp.readline()
        ptype= line.split()[0]
        pred = defaults[ptype]
        try:
            pred['descriptor']=line.split()[1]
        except IndexError:
            pass
        preds.append(pred)
    return subinp, preds

def get_genalg_params(subinp, line):
    pass
    return subinp, genalg_params

def readfile(subinp):
    '''this method reads all the inputkeywords'''
    #default values
    paras={'program':'gaussian',
           'procedure':'standard',
           'optimum':'minimum',
           'ml':0,
           'nosub':0,
           'nosub_file:':'',
           'cutoff':0,
           'no1sub':0,
           'try_ready':0,
           'test_ready':2,
           'regression':0,
           'tdregression':0,
           'difmodel':0,
           'nprocs':2,
           'debug':False,
           'norandom':0,
           'sequence':[],
           'symlinks':[],
           'restart':0,
           'semiempirical':0,
           'twojob':0,
           'multiplejobs':0,
           'stab':0,
           'polar':0,
           'ip':0,
           'ea':0,
           'startind': '',
           'extrawaittime': 2,
           'timelimit':250000,
           'timestep':300,
           'predictions':[],
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
        if line[0]=='#':continue
        elif 'bc' in line:
            paras['bcprop'] = line.split()[1]
            paras['bcval'] = line.split()[2]
            try:
                paras['bcoptimum'] = line.split()[3]
            except IndexError:
                paras['bcoptimum'] = 'min'
            assert paras['bcoptimum'] in ['min','Min','max','Max','MIN','MAX']
        elif 'nosub' in line:
            splitted = line.split()
            try:
                paras['nosub'] = int( splitted[1] )
            except IndexError:
                paras['nosub'] = 1
            if paras['nosub']==3:
                try:
                    paras['nosub_file'] = splitted[2]
                    logging.info( "nosub3. external file is used for data!: " + paras['nosub_file'])
                except IndexError:
                    logging.warning( "no file found. nosub downgraded to 1" )
                    paras['nosub'] = 1
        elif 'restart' in line:
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
        elif 'ml'==line[:2]:
            paras['ml']= int(line.split()[1])
        elif 'semiempirical' in line: paras['semiempirical'] = 1
        elif 'nprocs' in line: paras['nprocs'] = int(line.split()[1])
        elif 'debug' in line: paras['debug'] = True
        elif 'norandom' in line: paras['norandom'] = 1
        elif 'no1sub' in line: paras['no1sub'] = 1
        elif 'regression' in line: paras['regression'] = 1
        elif 'twodimreg' in line: paras['tdregression'] = 1
        elif 'difmodel' in line: paras['difmodel'] = 1
        elif 'optimum' in line: paras['optimum'] = line.split()[1]
        elif 'predictions' in line:
            subinp, paras['predictions'] = get_preds( subinp, line)
        elif 'sequence' in line:
            nsequences = int(line.split()[1])
            sequences = []
            for _ in range(nsequences):
                line = subinp.readline()
                sequence  = [ int(item) for item in line.split() ]
                sequences.append(sequence)
            paras['sequence'] = sequences[0]
            paras['sequences'] = sequences
        elif 'symlinks' in line: #NEW FEATURE! - not yet fully implemented
            nlinks = int(line.split()[1])
            links = []
            for _ in range(nlinks):
                line = subinp.readline()
                link = [ int(item) for item in line.split() ]
                assert len(link)>=2, 'link of len 1 is no link' # links larger than two could be allowed. 1 3 4 = 1-3 3-4 1-4
                links.append(link)
            paras['nlinks']   = nlinks
            paras['symlinks'] = links
            print "SYMMETRY ACTIVATED!"
        elif 'simple' in line: paras['simple'] = 1
        elif 'try_ready' in line: paras['try_ready'] = 1
        elif 'test_ready' in line: paras['test_ready'] = int(line.split()[1])
        elif 'twojob' in line: paras['twojob'] = 1
        elif 'startind' in line: paras['startind'] = line.split()[1]
        elif 'procedure' in line: 
                paras['procedure'] = line.split()[1]
                if paras['procedure'] in ['genrandom', 'getrandom']:
                    try:
                        paras['nrandom'] = int(line.split()[2])
                    except IndexError:
                        raise SystemExit("NO number of random structures specified!")
                elif paras['procedure'] in [ 'GA', 'genalg' ]:
                    #subinp, paras['genalg_params'] = get_genalg_params( subinp, line)
                    pass
        elif 'timelimit' in line: paras['timelimit'] = int(line.split()[1])
        elif 'timestep' in line: paras['timestep'] = int(line.split()[1])
        elif 'maxiter' in line: paras['maxiter'] = int(line.split()[1])
        elif 'mult' in line: paras['mult'] = int(line.split()[1])
        elif 'charge' in line: paras['charge'] = int(line.split()[1])
        elif 'basisset' in line: paras['basisset'] = line.split()[1]
        elif 'functional' in line: paras['functional'] = line.split()[1]
        elif 'identify' in line: paras['identify'] = line.split()[1]
        elif any(item in line for item in ('ncore','natomscore')): paras['ncore'] = int(line.split()[1])
        elif 'nch3' in line: paras['nch3'] = int(line.split()[1])
        elif 'montecarlo' in line:
            paras['montecarlo'] = float(line.split()[1])
            try:
                paras['nrandsites']= int(line.split()[2])
            except IndexError: pass
        elif 'extrawaittime' in line: paras['extrawaittime'] = float(line.split()[1])
        elif 'cutoff' in line: paras['cutoff'] = float(line.split()[1])
        elif 'property' in line: paras['property'] = line.split()[1]
        elif 'sites' in line: paras['line1'] = [ int(item) for item in line.split()[1:] ]
        elif 'positions' in line: paras['positions'] = [ int(item) for item in line.split()[1:] ]

        elif any(item in line.split()[0] for item in ('program','AI','Program','prog','Prog','programma')):
            if line.split()[1] in ['gaussian','g09','Gaussian','G09']:
                paras['program']= 'gaussian'
            elif line.split()[1] in ['orca','ORCA','Orca']:
                paras['program'] = 'orca'
            elif line.split()[1] in ['MOLPRO','Molpro','molpro']:
                paras['program'] = 'molpro'
                raise SystemExit('Molpro not yet implemented')
            else:
                raise SystemExit('program not recognized')
        elif 'END' in line: break
        else:
            print "line is not interpreted!", line

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
