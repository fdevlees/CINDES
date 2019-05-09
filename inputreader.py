
# inputreader module

import logging
import numpy as np
import random
from pprint import pprint
import re
from CINDES import run

from CINDES.utils.writings import log_io

# MAIN FUNCTION


@log_io()
def read_input(inputfilename='INPUT'):
    subinp = openfile(inputfilename)  # this is the fileID
    param = readfile(subinp)  # inputline is a tuple with all kind of input variables


    # here we set some extra parameters:
    # 1. if there is symmetry the real number of sites is smaller than the number of changeable sites
    if param['nlinks']:
        param['nsites'] = len(param['sites']) - param['nlinks']
    else:
        try:
            param['nsites'] = len(param['sites'])
        except KeyError:
            param['nsites'] = 0

    # here we set the fragment library per site, called array
    if param['procedure'] in ['genconf']:
        print "Generate Configuration Procedure Active"
        array = []
    else:
        array = substireader(param['nsites'], subinp)
        if not param['procedure'] in ['getrandom', 'genrandom', 'testpred']:
            print "ARRAY:"
            for i, item in enumerate(array):
                print "site{:>2d}:  |".format(i),
                for sub in item:
                    print " {} ".format("".join(sub)),
                    print "|",
                print
    if not param['procedure'] in ['getrandom', 'genrandom']:
        logging.debug("INPUT PARAMETERS:")
        for key, value in param.iteritems():
            logging.debug(key + ' : ' + str(value))

    # This is new and not yet fully functional
    param['array'] = array

    # even newer. There is only one common input-object. i.e. a Run object
    # so this can be initiated here:
    if param['procedure'] is None:
        myrun = run.BaseRun(**param)
    else:
        myrun = run.FrameRun(**param)

    return myrun


def openfile(filename):
    '''opens file in reading mode and returns fileid'''
    subinp = open(filename, 'r')
    return subinp


def get_preds(subinp, line):
    ''' for future development a more extensible format for giving which predictions are tried
    it returns a list of dictionaries with each dictionary having one obligatory type key '''
    # default predicition types:
    defaults = {'ml': {'type': 'ml', 'descriptor': 'coulomb'},
                'iml': {'type': 'iml'},
                'nn': {'type': 'nn', 'descriptor': 'bob'},
                '1d': {'type': '1d', 'descriptor': '1DL', 'subtype': 'ridge', 'intercept': True},
                '2d': {'type': '2d'},
                'knn': {'type': 'knn'},
                'gp': {'type': 'gp'},
                'svr': {'type': 'svr'},
                'krr': {'type': 'krr', 'kernel': 'rbf'},
                'qml': {'type': 'qml', 'kernel': 'rbf'}
                }
    # set n_folds default for each experiment:
    for experiment in defaults.values():
        experiment.update({'n_folds': 5,
                           'pca': False,
                           'n_principal_components': 100,
                           'plots': [],
                           'tableindex': 1})

    npredictions = int(line.split()[1])
    preds = []  # this becomes a list of predictions to make

    for _ in range(npredictions):
        while True:
            line = subinp.readline()
            if not '#' in line:
                break
        # the first word is a unique prediction identifier (just a name which has to be unique)
        pname = line.split()[0]
        ptype = line.split()[1]   # the second word indicates the prediction type
        pred = defaults[ptype].copy()   # the defaults for that prediction type are then loaded in pred
        pred['name'] = pname     # set the name

        # the rest of the line is than interpreted: add new arguments or change default arguments
        try:
            # these lines:
            #    - splits the rest of the line in keyword
            #    - adds apostrophs around the keys
            #    - join the key:value pairs with comma's
            splitted = [item for item in line.split()[2:]]
            formatted = ['\'{}\':{}'.format(*item.split(':')) for item in splitted]
            options = ','.join(formatted)
        except IndexError:
            pass
        if options:
            #print "options:", options
            pred.update(eval('{{{}}}'.format(options)))
            #print "prediction keywords are changed:", pred

        # the fully declared prediction type is than saved to the prediction list
        preds.append(pred)
    return subinp, preds


def get_prop_function(subinp, line):
    line = subinp.readline()
    isword = re.compile('(^[a-zA-Z_]+[0-9]?\(?$)')
    splitted1 = re.split('(\[|\])', line)
    splitted2 = map(lambda x:x.split(), splitted1)

    props = set()
    inkey = False
    for _split in splitted2:
        for split in _split:
            if split=='[': inkey=True
            elif split==']':inkey=False
            if inkey:
                continue
            if isword.match(split) and (not split in ['if', 'else', 'abs(']):
                props.add(split)
    arguments = ",".join(props)
    funcstr = "lambda {}:{}".format(arguments, line)
    print "funcstr:", funcstr
    func = eval(funcstr)
    func.__doc__ = 'lambda function: {}'.format(line.strip())

    return subinp, func, props

def get_calcs(subinp, line):
    supercalcs = []
    calcs = []
    i, j=1, 1
    line = subinp.readline()
    while True:
        assert line.strip()=='{:d}.{:d}'.format(i,j)
        line = subinp.readline()
        calcs.append(get_jobs(subinp, line, index=0))
        # try to see if yet another job is given
        line = subinp.readline().strip()
        if line=='endcalcs' or line.split('.')[0]==str(i+1):
            if len(calcs)==1:
                supercalcs.append(calcs[0])
            else:
                supercalcs.append(calcs)
            calcs = []
        if line=='endcalcs':
            break
        i, j = map(int, line.split('.'))
    return supercalcs


def get_jobs(subinp, line, index=1):
    def get_extra_line(line):
        key, value = line.strip().split(None, 1)
        assert key in ['identify', 'nosub', 'program', 'nprocs', 'mem',
                'geom', 'script', 'positions', 'fafoom', 'rdfreq']
        if key in ['nosub', 'nprocs', 'fafoom', 'mem']:
            value = int(value)
        elif key in ['positions']:
            value = map(int, value.split())
        calc[key] = value
        return
    calc = dict()
    splitted = line.split()
    njobs = int(splitted[index])
    if len(splitted) == index+2:
        n_extra_lines = int(splitted[index+1])
        for _ in range(n_extra_lines):
            line = subinp.readline()
            get_extra_line(line)
    jobs = []
    for _ in range(njobs):
        job = dict()
        # read propline
        line = subinp.readline().split()
        job['info'] = set(line)
        # read mult/charge/hotline
        line = subinp.readline().split()
        job['charge'], job['mult'], job['hotline'] = (line[0], line[1], ' '.join(line[2:]))
        jobs.append(job)
    calc['jobs'] = jobs
    return calc


def get_genalg_params(subinp, line):
    defaults = {'ngenerations': 20,
                'npopulation': 20,
                'CXP': 0.0,  # crossover probability
                'MUP': 0.2,  # mutation probability
                'nelitism': 1,
                'scaling': 'sigmatrunc',
                'db_identify': 'ex' + str(np.random.randint(0, 90)),
                'freq_stats': 10,
                'restart':False,
                'seed': 0,
                'selector': 'RouletteWheel',
                'tournamentPoolsize':2,
                'weights': (1.0,),
                'algorithm':'GA'
                }
    try:
        n_extra_lines = int(line.split()[2])
    except IndexError:
        print "WARNING: all default values for the genetic algorithms will be used:"
    else:  # execute only when no exception is thrown
        for _ in range(n_extra_lines):
            line = subinp.readline()
            key = line.split()[0]
            if not key in defaults:
                print "keyword in GA section not recognized:", key
                raise SystemExit('program stopped')
            value_type = type(defaults[key])
            if value_type is bool:
                defaults[key] = True
            elif value_type is tuple:
                defaults[key] = tuple(map(float, line.split()[1:]))
            else:
                defaults[key] = value_type(line.split()[1])
        print " defaults of genetic algorithm are changed. new values:"
    return subinp, defaults

def get_pso_params(subinp, line):
    defaults = {'ngenerations': 20,
                'npopulation': 20,
                'w1': 1.0,  # local optimum factor
                'w2': 1.0,  # global optimum factor
                'w':0.8,
                'c1': 0.2,  # random factor1
                'c2': 0.2,  # random factor2
                'epsilon': 0.7,  # random factor2
                'db_identify': 'ex' + str(random.randint(0, 90)),
                'freq_stats': 10,
                'type':'concrete'
                }
    try:
        n_extra_lines = int(line.split()[2])
    except IndexError:
        print "WARNING: all default values for the genetic algorithms will be used:"
    else:  # execute only when no exception is thrown
        for _ in range(n_extra_lines):
            line = subinp.readline()
            key = line.split()[0]
            if not key in defaults:
                print "keyword in PSO section not recognized:", key
                raise SystemExit('program stopped')
            value_type = type(defaults[key])
            defaults[key] = value_type(line.split()[1])
        print " defaults for particle swarm optimization are changed. new values:"

    # change defaults of c1/c2 for PPSO
    if not defaults['type'] == 'concrete':
        for key in ['c1', 'c2']:
            if defaults[key]==0.2:
                defaults[key]=1.4

    return subinp, defaults

def get_gprf_params(subinp, line):
    defaults = {'algorithm': 'GP',
                'acq_func': 'gp_hedge',
                'acq_optimizer': 'sampling',
                'n_initial_points': 10,
                'seed':1287294368,
                'batchsize': 20,
                }
    try:
        n_extra_lines = int(line.split()[2])
    except IndexError:
        print "WARNING: all default values for the genetic algorithms will be used:"
    else:  # execute only when no exception is thrown
        for _ in range(n_extra_lines):
            line = subinp.readline()
            key = line.split()[0]
            if not key in defaults:
                print "keyword in GP/RF section not recognized:", key
                raise SystemExit('program stopped')
            value_type = type(defaults[key])
            defaults[key] = value_type(line.split()[1])
        print " defaults for GP/RF are changed. new values:", defaults
    return subinp, defaults


def readfile(subinp):
    '''this method reads all the inputkeywords'''
    # default values
    randomseed = np.random.randint(0, 100)
    paras = {

        #          GLOBAL RUN PARAMETERS
        'adjust_dihedrals': False,
        'batchsize': None,
        'bc': False,
        'calcs':[],
        'cutoff': 0,  # this cutoff has to apply to the final target property value only
        'debug': False,
        'defaultgroups':None,
        'difmodel': 0,
        'extrajobs': [],
        'extra_props': [],
        'extrawaittime': 2,
        'function': lambda x: x,
        'ignore': 0,
        'jobs': [],
        'logging':logging.INFO,
        'ml': 0,
        'maxiter': 10,
        'montecarlo': 0,  # Temperature at start
        'nch3': Ellipsis,
        'ncore': Ellipsis,
        'nlinks': False,
        'no1sub': 0,
        'norandom': 0,
        'nrandsites': 2,
        'optimum': 'minimum',
        'optga': False,
        'predictions': [],
        'prejobs': None,
        'procedure': None,
        'property': 'gap',
        'regression': 0,
        'restart': 0,
        'readtable':False,
        'restrictions': [],
        'seed': randomseed,
        'sequence': [],
        'startind': '',
        'symlinks': [],
        'tablename': 'table.json',
        'test_ready': 2,
        'threading': False,
        'timelimit': 250000,
        'timestep': 300,
        'try_ready': 0,
        'secret_file': '',
        'worker': False,
        'write':True,
        'TZmat': {},
        'zmatrixfile': 'ZMAT',

        #          LOCAL JOB PARAMETERS
        'geom2': None,
        'identify': 'unspecified_',
        'nprocs': 2,
        'program': 'gaussian',
        'stab': 0,
        'stabjobs': [],
        'nosub': 0
    }  # n random sites changed. for all choose 0
    # scans all the lines until if will find the END keyword
    # this is a bit tricky because keywords can appear everywere in the file before END
    # but the value of the keywords is always on the second and further position
    # i can if i want change all lines with:
    # if 'mystring' in line:
    # substitute by
    # if 'mystr' in line.split()[0]:
    while True:
        line = subinp.readline()
        splitted = line.split()
        if line == '\n':
            continue
        if line[0] == '#':
            continue
        key = splitted[0]
        # 1. some capital sensitive keywords:
        if 'tablename' in line:
            paras['tablename'] = splitted[1]
            continue
        elif 'startind' in line:
            paras['startind'] = splitted[1]
            continue
        elif 'nosub' in line:
            try:
                paras['nosub'] = int(splitted[1])
            except IndexError:
                paras['nosub'] = 1
            if paras['nosub'] == 3:
                raise SyntaxError('this functionality is renamed to: secret_file <filename>')
            continue
        elif 'secret_file' in line:
            paras['secret_file'] = splitted[1]
            logging.info("nosub3. external file is used for data!: " + paras['secret_file'])
            continue
        elif 'zmatrixfile' in line:
            if len(splitted[1:])==1:
                paras['zmatrixfile'] = splitted[1]
            else:
                print "multiple ZMAT files!"
                paras['zmatrixfile'] = splitted[1:]
            continue
        elif key=='defaultgroups':
            ndefaults = int(splitted[1])
            defaults = {}
            for _ in range(ndefaults):
                line = subinp.readline().split()
                defaults[int(line[0])] = line[1]
            paras['defaultgroups']=defaults
            continue
        elif 'END' in line:
            break

        # 2. capital insensitive keywords:
        line = line.split('#')[0].lower()
        if 'adjust_dihedrals' in line:
            paras['adjust_dihedrals'] = True
        elif 'batchsize' in line:
            paras['batchsize'] = int(line.split()[1])
        elif 'bc' in line:
            paras['bc'] = True
            paras['bcprop'] = line.split()[1]
            paras['bcval'] = line.split()[2]
            try:
                paras['bcoptimum'] = line.split()[3]
            except IndexError:
                paras['bcoptimum'] = 'minimum'
            if paras['bcoptimum'] in ['min', 'minimum']:
                paras['bcoptimum'] = 'min'
            else:
                paras['bcoptimum'] = 'max'
        elif 'startcalcs' in line:
            paras['calcs'] = get_calcs(subinp, line)
        # elif 'basisset' in line: paras['basisset'] = line.split()[1]
        elif 'cutoff' in line:
            paras['cutoff'] = float(line.split()[1])
        # elif 'charge' in line: paras['charge'] = int(line.split()[1])
        elif 'debug' in line:
            paras['debug'] = True
        elif 'difmodel' in line:
            paras['difmodel'] = 1
        elif 'extra_props' in line:
            paras['extra_props'] = line.split()[1:]
        elif 'extrajobs' in line:
            paras['extrajobs'] = get_jobs(subinp, line)
        elif 'extrawaittime' in line:
            paras['extrawaittime'] = float(line.split()[1])
        elif 'ignore' in line:
            try:
                paras['ignore'] = int(line.split()[1])
            except IndexError:
                paras['ignore'] = 3600
        elif 'stabjobs' in line:
            #raise NotImplementedError('this is not updated to the multiprogram scheme')
            # this code has to come before 'jobs' because also 'jobs' in 'stabjobs'
            njobs = int(line.split()[1])
            jobs = []
            for _ in range(njobs):
                job = dict()
                # read propline
                line = subinp.readline().split()
                job['info'] = set(line)
                # read mult/charge/hotline
                line = subinp.readline().split()
                job['charge'], job['mult'], job['hotline'] = (line[0], line[1], ' '.join(line[2:]))
                jobs.append(job)
                paras['stabjobs'] = jobs
        elif splitted[0] == 'jobs':
            paras['jobs'] = get_jobs(subinp, line)
        elif splitted[0] == 'logging':
            level = splitted[1].lower()
            if level == 'debug':
                logging.getLogger().setLevel(logging.DEBUG)
            elif level == 'info':
                logging.getLogger().setLevel(logging.INFO)
            elif level == 'warning':
                logging.getLogger().setLevel(logging.WARNING)
            else:
                raise SystemExit('logging level not recognized')
            logging.info('set logging level to:' + level)
        elif 'maxcycles' in line:
            paras['maxcycles'] = str(int(line.split()[1]))
        elif 'maxiter' in line:
            paras['maxiter'] = int(line.split()[1])
        elif 'montecarlo' in line:
            paras['montecarlo'] = float(line.split()[1])
            try:
                paras['nrandsites'] = int(line.split()[2])
            except IndexError:
                pass
        elif key=='mult':
            paras['mult'] = int(line.split()[1])
        elif 'ml' == line[:2]:
            paras['ml'] = int(line.split()[1])
        elif 'norandom' in line:
            paras['norandom'] = 1
        elif 'no1sub' in line:
            paras['no1sub'] = 1
        elif any(item in line for item in ('ncore', 'natomscore')):
            paras['ncore'] = int(line.split()[1])
        elif 'nch3' in line:
            paras['nch3'] = int(line.split()[1])
        elif 'nprocs' in line:
            paras['nprocs'] = int(line.split()[1])
        elif 'optimum' in line:
            if 'max' in line.split()[1]:
                paras['optimum'] = 'maximum'
            else:
                paras['optimum'] = 'minimum'
        elif 'optga' in line:
            paras['optga'] = True
        elif 'positions' in line:
            paras['positions'] = [int(item) for item in line.split()[1:]]
        elif 'predictions' in line:
            subinp, paras['predictions'] = get_preds(subinp, line)
        elif 'prejobs' in line:
            paras['prejobs'] = get_jobs(subinp, line)
        elif 'property' in line:
            prop = line.split()[1]
            try:
                functionfile = line.split()[2]
            except IndexError:
                functionfile = 'function'
            if 'load_func' in prop:
                paras['property'] = 'func'
                functionscript = __import__(functionfile)
                paras['function'] = functionscript.function
                assert callable(paras['function'])
                paras['func_args'] = functionscript.arguments
            elif 'func' in prop:
                paras['property'] = 'func'
                subinp, paras['function'], paras['func_args'] = get_prop_function(subinp, line)
            else:
                paras['property'] = prop
        elif 'procedure' in line:
            paras['procedure'] = line.split()[1]
            if paras['procedure'] in ['genrandom', 'getrandom']:
                try:
                    paras['nrandom'] = int(line.split()[2])
                except IndexError:
                    raise SystemExit("NO number of random structures specified!")
            elif paras['procedure'] in ['ga', 'genalg', 'deap']:
                subinp, paras['genalg'] = get_genalg_params(subinp, line)
            elif paras['procedure'] in ['pso', 'particleswarm']:
                subinp, paras['pso'] = get_pso_params(subinp, line)
            elif paras['procedure'] in ['gp', 'rf']:
                subinp, paras['gprf'] = get_gprf_params(subinp, line)
            elif paras['procedure'] in ['testpred', 'makepred']:
                try:
                    paras['datacolumn'] = int(line.split()[2])
                except IndexError:
                    print "first column of table is taken as datacolumn. (default)"
                    paras['datacolumn'] = 1
                paras['restart'] = 1
                paras['readtable'] = True
            elif paras['procedure'] == 'getdivers':
                paras['divers_nmax'] = int(line.split()[2])
                paras['divers_batchsize'] = int(line.split()[3])
                paras['divers_divindex'] = int(line.split()[4])
            elif paras['procedure'] in ['generate']:
                try:
                    paras['ngenerate'] = int(line.split()[2])
                except IndexError:
                    print "all possible molecules from site array will be calculated!"
                else: # do only when no error catched
                    indices = []
                    for _ in range(paras['ngenerate']):
                        line = subinp.readline()
                        indices.append(line.strip())
                    paras['generatemols'] = indices
                    print "read {:d} indices to generate".format(paras['ngenerate'])
        elif 'regression' in line:
            paras['regression'] = 1
        elif 'readtable' in line:
            paras['readtable'] = True
        elif 'restart' in line:
            print "RESTART KEYWORD IS DEPRECATED: USE READTABLE INSTEAD"
            paras['restart'] = int(line.split()[1])
            if paras['restart'] >= 1:
                paras['readtable'] = True
        elif 'restrictions' in line:
            paras['restrictions'] = [int(item) for item in line.split()[1:]]
        elif 'seed' in line:
            paras['seed'] = int(line.split()[1])
        # elif 'semiempirical' in line: paras['semiempirical'] = 1
        elif 'sequence' in line:
            nsequences = int(line.split()[1])
            sequences = []
            for _ in range(nsequences):
                # that is: read line and replace the chars in '[,]' by None.
                line = subinp.readline().translate(None, '[,]')
                sequence = [int(item) for item in line.split()]
                sequences.append(sequence)
            paras['sequence'] = sequences[0]
            paras['sequences'] = sequences
        elif 'symlinks' in line:  # NEW FEATURE! - not yet fully implemented
            nlinks = int(line.split()[1])
            links = []
            for _ in range(nlinks):
                line = subinp.readline()
                link = [int(item) for item in line.split()]
                # links larger than two could be allowed. 1 3 4 = 1-3 3-4 1-4
                assert len(link) >= 2, 'link of len 1 is no link'
                links.append(link)
            paras['nlinks'] = nlinks
            paras['symlinks'] = links
            print "SYMMETRY ACTIVATED!"
        elif 'simple' in line:
            paras['simple'] = 1
        elif 'sites' in line:
            paras['sites'] = [int(item) for item in line.split()[1:]]
        elif 'try_ready' in line:
            paras['try_ready'] = 1
        elif 'test_ready' in line:
            paras['test_ready'] = int(line.split()[1])
        elif 'threading' in line:
            paras['threading'] = True
        # elif 'twojob' in line:
        #    try:
        #        paras['twojob'] = int(line.split()[1])
        #    except IndexError:
        #        paras['twojob'] = 1
        elif 'procedure' in line:
            paras['procedure'] = line.split()[1]
            if paras['procedure'] in ['genrandom', 'getrandom']:
                try:
                    paras['nrandom'] = int(line.split()[2])
                except IndexError:
                    raise SystemExit("NO number of random structures specified!")
            elif paras['procedure'] in ['ga', 'genalg']:
                subinp, paras['genalg'] = get_genalg_params(subinp, line)
                pass
        elif 'timelimit' in line:
            paras['timelimit'] = int(line.split()[1])
        elif 'timestep' in line:
            paras['timestep'] = int(line.split()[1])

        # LOCAL FUTURE JOB KEYWORDS:
        elif 'geom2' in line:
            paras['geom2'] = line.split()[1]
        elif 'identify' in line:
            paras['identify'] = line.split()[1]
            if not paras['identify'].endswith('_'):
                print "underscore appended to identify"
                paras['identify'] += '_'
            continue
        elif key=='worker':
            paras['worker'] = True
        elif key=='write':
            paras['write'] = bool(int(line.split()[1]))
        elif any(item in line.split()[0] for item in ('program', 'ai', 'program', 'prog', 'programma')):
            if line.split()[1] in ['gaussian', 'g09']:
                paras['program'] = 'gaussian'
            elif line.split()[1] in ['orca']:
                paras['program'] = 'orca'
            elif line.split()[1] in ['molpro']:
                paras['program'] = 'molpro'
                raise SystemExit('Molpro not yet implemented')
            elif line.split()[1] in ['nwchem']:
                paras['program'] = 'nwchem'
            else:
                raise SystemExit('program not recognized')

        else:
            if line.strip():  # so if not just an empty line:
                print "line: \"{}\" is not interpreted".format(line.strip('\n')),
                raise SystemExit('program stopped')

    # get a list of all properties that need to be calculated:
    if paras['property'] == 'func':
        props = []
        props.extend(paras['func_args'])
    elif paras['property'] == '_':
        props = []
    else:
        props = [paras['property']]
    try:
        props.append(paras['bcprop'])
    except KeyError:
        pass
    props.extend(paras['extra_props'])
    paras['props'] = set(props)
    for prop in ['stab']:
        if prop in paras['props']:
            paras[prop] = True

    return paras


def substireader(nsit, subinp):
    '''this one reads all the different substituents for all different positions'''
    substiarray = []
    for i in range(nsit):
        # read number of substituents
        try:
            line2 = subinp.readline().split()
            nsubsit = int(line2[1])
        except IndexError:
            print "no functional groups present or wrong formatted"
            break
        logging.debug("site number: " + str(i + 1))
        logging.debug("nsubsit: " + str(nsubsit))
        site = [subinp.readline().split() for line in range(nsubsit)]
        # for j in range(nsubsit):
        #    atoms = subinp.readline().split()
        substiarray.append(site)
        # read each substituent in an array or so
    return(substiarray)


if __name__ == "__main__":
    import sys
    fid = openfile(sys.argv[1])
    line1, line2 = readfile(fid)
    print "inputline: ", line1
    print "tweede inputline: ", line2
    nsit = len(line1[-1])
    subs = substireader(nsit, fid)
    print "subs:",
    pprint(subs)
