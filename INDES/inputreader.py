

# inputreader module

import logging
import numpy as np
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
    if param['nlinks']:
        param['nsites'] = len(param['sites']) - param['nlinks']
    else:
        param['nsites'] = len(param['sites'])
    #param['nsites'] = len(param['line1'])
    #####################
    if param['procedure'] in [ 'genconf' ]:
        print "Generate Configuration Procedure Active"
        array = []
    else:
        array = substireader(param['nsites'],subinp)
        if not param['procedure'] in ['getrandom', 'genrandom','testpred']:
            print "ARRAY:"
            for i,item in enumerate(array):
                print "site{:>2d}:  |".format(i),
                for sub in item:
                    print " {} ".format("".join(sub)),
                    print "|",
                print
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
    defaults = { 'ml' : { 'type': 'ml', 'descriptor':'coulomb'},
                 'iml': { 'type':'iml' },
                 'nn' : { 'type': 'nn', 'descriptor':'bob'},
                 '1d' : { 'type': '1d', 'descriptor':'1DL', 'subtype':'ridge', 'intercept':True },
                 '2d' : { 'type': '2d' },
                 'knn': { 'type': 'knn'},
                 'gp' : { 'type': 'gp' },
                 'svr': { 'type': 'svr'},
                 'krr': { 'type': 'krr', 'kernel':'rbf'},
                 'qml': { 'type': 'qml', 'kernel':'rbf'}
               }
    # set n_folds default for each experiment:
    for experiment in defaults.values(): experiment.update( {'n_folds':5 ,
                                                                 'pca':False,
                                              'n_principal_components':100,
                                                               'plots':[],
                                                          'tableindex':1 } )

    npredictions = int(line.split()[1])
    preds = [] # this becomes a list of predictions to make

    for _ in range(npredictions):
        while True:
            line = subinp.readline()
            if not '#' in line: break
        pname= line.split()[0]   # the first word is a unique prediction identifier (just a name which has to be unique)
        ptype= line.split()[1]   # the second word indicates the prediction type
        pred = defaults[ptype].copy()   # the defaults for that prediction type are then loaded in pred
        pred['name'] = pname     # set the name 

        # the rest of the line is than interpreted: add new arguments or change default arguments
        try:
            # these lines:
            #    - splits the rest of the line in keyword
            #    - adds apostrophs around the keys
            #    - join the key:value pairs with comma's
            splitted = [ item for item in line.split()[2:] ]
            formatted= [ '\'{}\':{}'.format(*item.split(':')) for item in splitted ]
            options = ','.join(formatted)
        except IndexError:
            pass
        if options:
            #print "options:", options
            pred.update( eval( '{{{}}}'.format(options) ) )
            #print "prediction keywords are changed:", pred

        # the fully declared prediction type is than saved to the prediction list
        preds.append(pred)
    return subinp, preds

def get_prop_function(subinp, line):
    line = subinp.readline()
    splitted = line.split()
    #print "functional property:", line
    import re

    # we need to find the properties going into the function. properties only contain [a-zA-Z]
    word = re.compile('(^[a-zA-Z_]*$)')

    # the properties in that line are: #set because one property can occur multiple times in function
    props = { item for item in splitted if word.match(item) and not item in ['if', 'else' ] }
    #print "properties:", props

    # props need to be separated by a comma
    arguments = ','.join(props)
    #print "arguments:", arguments
    func = eval('lambda {}:{}'.format(arguments, line))

    return subinp, func, props

def get_jobs(subinp, line):
    njobs = int(line.split()[1])
    jobs=[]
    for _ in range(njobs):
        job=dict()
        # read propline
        line = subinp.readline().split()
        job['info']=set(line)
        # read mult/charge/hotline
        line = subinp.readline().split()
        job['charge'], job['mult'], job['hotline'] = (line[0], line[1], ' '.join(line[2:]))
        jobs.append(job)
    return jobs

def get_genalg_params(subinp, line):
    defaults = { 'ngenerations' : 20,
                 'npopulation'  : 20,
                 'CXP'          : 0.8, #crossover probability
                 'MUP'          : 0.2, #mutation probability
                 'elitism'      : True,
                 'optimum'      : 'maximum',
                 'nelitism'     : 1,
                 'scaling'      : 'sigmatrunc',
                 'db_identify'  : 'ex' + str(np.random.randint(0,90)),
                 'freq_stats'   : 10,
                 'seed'         : 0,
                 'selector'     : 'RouletteWheel'
                 }
    try:
        n_extra_lines = int(line.split()[2])
    except IndexError:
        print "WARNING: all default values for the genetic algorithms will be used:"
    else: # execute only when no exception is thrown
        for _ in range(n_extra_lines):
            line = subinp.readline()
            key = line.split()[0]
            if not key in defaults:
                print "keyword in GA section not recognized:", key
                raise SystemExit('program stopped')
            value_type = type(defaults[key])
            defaults[key] = value_type( line.split()[1] )
            #if key in [ 'ngenerations', 'npopulation', 'nelitism', 'freq_stats' ]:
            #    defaults[key] = int(line.split()[1])
            #elif key in ['CXP', 'MUP']:
            #    defaults[key] = float(line.split()[1])
            #elif key in ['elitism']:
            #    defaults[key] = bool(line.split()[1])
            #elif key in [ 'scaling', 'db_identify', 'optimum' ]:
            #    defaults[key] = line.split()[1]
            #else:
            #    print "line not interpreted:", line
        print " defaults of genetic algorithm are changed. new values:"
    print defaults
    return subinp, defaults

def readfile(subinp):
    '''this method reads all the inputkeywords'''
    #default values
    randomseed = np.random.randint(0,100)
    paras={
           'adjust_dihedrals':False,
           #'aea':0,
           #'aip':0,
           #'basisset':'6-31G',
           'bc': False,
           'charge':0,
           'cutoff':0,
           'debug':False,
           'difmodel':0,
           #'ea':0,
           'extrajobs':[],
           'extra_props': [],
           'extrawaittime': 2,
           #'functional':'b3lyp',
           'function': lambda x:x,
           'identify':'unspecified_',
           #'ip':0,
           'jobs':[],
           'maxiter':10,
           'maxcycles':'100',
           #'ml':0,
           'montecarlo':0,   #Temperature at start
           'mult':1,
           #'multiplejobs':0,
           'nch3':16,
           'ncore':10,
           'nlinks':False,
           'no1sub':0,
           'norandom':0,
           'nosub':0,
           'nosub_file:':'',
           'nprocs':2,
           'nrandsites':2,
           'optimum':'minimum',
           'optga':False,
           #'polar':0,
           'predictions':[],
           'procedure':'standard',
           'program':'gaussian',
           'property':'gap',
           'regression':0,
           'restart':0,
           'restrictions':[],
           'seed':randomseed,
           #'semiempirical':0,
           'sequence':[],
           #'solv':False,
           'stab':0,
           'stabjobs':[],
           'startind': '',
           'symlinks':[],
           'tablename':'table.json',
           'tdregression':0,
           'test_ready':2,
           'timelimit':250000,
           'timestep':300,
           'try_ready':0
           #'twojob':0,
           #'volume':False
           }   #n random sites changed. for all choose 0
    #scans all the lines until if will find the END keyword
    #this is a bit tricky because keywords can appear everywere in the file before END 
    #but the value of the keywords is always on the second and further position
    # i can if i want change all lines with:
    # if 'mystring' in line:
    # substitute by
    # if 'mystr' in line.split()[0]:
    while True:
        line = subinp.readline()
        if line == '\n':continue
        if line[0]=='#':continue
        # 1. some capital sensitive keywords:
        elif 'tablename' in line:
            paras['tablename'] = line.split()[1]
            continue
        elif 'startind' in line:
            paras['startind'] = line.split()[1]
            continue
        elif 'identify' in line:
            paras['identify'] = line.split()[1]
            continue
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
            continue
        elif 'END' in line: break

        # 2. capital insensitive keywords:
        line = line.split('#')[0].lower()
	if 'adjust_dihedrals' in line: paras['adjust_dihedrals']=True
        elif 'bc' in line:
            paras['bc'] = True
            paras['bcprop'] = line.split()[1]
            paras['bcval'] = line.split()[2]
            try:
                paras['bcoptimum'] = line.split()[3]
            except IndexError:
                paras['bcoptimum'] = 'min'
            assert paras['bcoptimum'] in ['min','max']
        elif 'basisset' in line: paras['basisset'] = line.split()[1]
        elif 'cutoff' in line: paras['cutoff'] = float(line.split()[1])
        elif 'charge' in line: paras['charge'] = int(line.split()[1])
        elif 'debug' in line: paras['debug'] = True
        elif 'difmodel' in line: paras['difmodel'] = 1
        elif 'extra_props' in line: paras['extra_props'] = line.split()[1:]
        elif 'extrajobs' in line:
            paras['extrajobs']=get_jobs(subinp, line)
        elif 'extrawaittime' in line: paras['extrawaittime'] = float(line.split()[1])
        elif 'functional' in line: paras['functional'] = line.split()[1]
	elif 'stabjobs' in line:
            # this code has to come before 'jobs' because also 'jobs' in 'stabjobs'
            njobs = int(line.split()[1])
            jobs=[]
            for _ in range(njobs):
                job=dict()
                # read propline
                line = subinp.readline().split()
                job['info']=set(line)
                # read mult/charge/hotline
                line = subinp.readline().split()
                job['charge'], job['mult'], job['hotline'] = (line[0], line[1], ' '.join(line[2:]))
                jobs.append(job)
                paras['stabjobs']=jobs
        elif 'jobs' in line:
            paras['jobs']=get_jobs(subinp, line)
        elif 'maxcycles' in line: paras['maxcycles']= str(int(line.split()[1]))
        elif 'maxiter' in line: paras['maxiter'] = int(line.split()[1])
        elif 'montecarlo' in line:
            paras['montecarlo'] = float(line.split()[1])
            try:
                paras['nrandsites']= int(line.split()[2])
            except IndexError: pass

        #-----
        # do not interchange the following two keywords!!!"
        elif 'multiplejobs' in line:
            print "set multiple jobs:"
            paras['multiplejobs'] = int(line.split()[1])
        elif 'mult' in line: paras['mult'] = int(line.split()[1])
        #-----

        elif 'ml'==line[:2]:
            paras['ml']= int(line.split()[1])
        elif 'norandom' in line: paras['norandom'] = 1
        elif 'no1sub' in line: paras['no1sub'] = 1
        elif any(item in line for item in ('ncore','natomscore')): paras['ncore'] = int(line.split()[1])
        elif 'nch3' in line: paras['nch3'] = int(line.split()[1])
        elif 'nprocs' in line: paras['nprocs'] = int(line.split()[1])
        elif 'optimum' in line:
            if 'max' in line.split()[1]:
                print "changed optimization to maximum instead of minimum!"
                paras['optimum'] = 'maximum'
	elif 'optga' in line: paras['optga']=True
        elif 'positions' in line: paras['positions'] = [ int(item) for item in line.split()[1:] ]
        elif any(item in line.split()[0] for item in ('program','ai','program','prog','programma')):
            if line.split()[1] in ['gaussian','g09']:
                paras['program']= 'gaussian'
            elif line.split()[1] in ['orca']:
                paras['program'] = 'orca'
            elif line.split()[1] in ['molpro']:
                paras['program'] = 'molpro'
                raise SystemExit('Molpro not yet implemented')
            elif line.split()[1] in ['nwchem']:
                paras['program'] = 'nwchem'
            else:
                raise SystemExit('program not recognized')
        elif 'predictions' in line:
            subinp, paras['predictions'] = get_preds( subinp, line)
        elif 'property' in line:
            prop = line.split()[1]
            if 'func' in prop:
                paras['property']='func'
                subinp, paras['function'], paras['func_args'] = get_prop_function( subinp, line )
            else:
                paras['property'] = prop
        elif 'procedure' in line:
                paras['procedure'] = line.split()[1]
                if paras['procedure'] in ['genrandom', 'getrandom']:
                    try:
                        paras['nrandom'] = int(line.split()[2])
                    except IndexError:
                        raise SystemExit("NO number of random structures specified!")
                elif paras['procedure'] in [ 'ga', 'genalg' ]:
                    subinp, paras['genalg'] = get_genalg_params( subinp, line)
                    pass
                elif paras['procedure'] in [ 'testpred', 'makepred' ]:
                    try:
                        paras['datacolumn'] = int(line.split()[2])
                    except IndexError:
                        print "first column of table is taken as datacolumn. (default)"
                        paras['datacolumn'] = 1
                    paras['restart'] = 1
                elif paras['procedure']=='getdivers':
                    paras['divers_nmax'] = int(line.split()[2])
                    paras['divers_batchsize'] = int(line.split()[3])
                    paras['divers_divindex'] = int(line.split()[4])
        elif 'regression' in line: paras['regression'] = 1
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
        elif 'restrictions'  in line: paras['restrictions'] = [ int(item) for item in line.split()[1:] ]
        elif 'seed' in line: paras['seed'] = int(line.split()[1])
        elif 'semiempirical' in line: paras['semiempirical'] = 1
        elif 'sequence' in line:
            nsequences = int(line.split()[1])
            sequences = []
            for _ in range(nsequences):
                line = subinp.readline().translate(None,'[,]') # that is: read line and replace the chars in '[,]' by None. 
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
        elif 'sites' in line: paras['sites'] = [ int(item) for item in line.split()[1:] ]
        elif 'twodimreg' in line: paras['tdregression'] = 1
        elif 'try_ready' in line: paras['try_ready'] = 1
        elif 'test_ready' in line: paras['test_ready'] = int(line.split()[1])
        elif 'twojob' in line:
            try:
                paras['twojob'] = int(line.split()[1])
            except IndexError:
                paras['twojob'] = 1
        elif 'procedure' in line:
                paras['procedure'] = line.split()[1]
                if paras['procedure'] in ['genrandom', 'getrandom']:
                    try:
                        paras['nrandom'] = int(line.split()[2])
                    except IndexError:
                        raise SystemExit("NO number of random structures specified!")
                elif paras['procedure'] in [ 'ga', 'genalg' ]:
                    subinp, paras['genalg'] = get_genalg_params( subinp, line)
                    pass
        elif 'timelimit' in line: paras['timelimit'] = int(line.split()[1])
        elif 'timestep' in line: paras['timestep'] = int(line.split()[1])
        else:
            if line.strip(): # so if not just an empty line:
                print "line: \"{}\" is not interpreted".format(line.strip('\n')),
                raise SystemExit('program stopped')

    # get a list of all properties that need to be calculated:
    if paras['property']=='func':
        props = []
        props.extend( paras['func_args'] )
    else:
        props = [ paras['property'] ]
    try:
        props.append(paras['bcprop'])
    except KeyError:
        pass
    props.extend( paras['extra_props'] )
    paras['props'] = set(props)
    #for prop in ['stab', 'polar', 'ip', 'aip', 'ea', 'aea', 'solv' ]:
    for prop in ['stab']:
        if prop in paras['props']:
            paras[prop]=True

    # extra sanity checks on input
    # sanity check 1: optimum in ga and bfs input similar
    if hasattr(paras,'genalg'):
        if not paras['genalg']['optimum'] == paras['optimum']:
            print "WARNING OPTIMUM KEYWORDS ARE NOT THE SAME. \n    Please check carefully if program is working correctly!"
            print "set optimum to genalg.optimum", paras['genalg']['optimum']
            paras['optimum'] = paras['genalg']['optimum']

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
