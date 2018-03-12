#!/bin/env python 

# debug flag
debug=1

# import python libraries
import time # for getting time/date and time delays
start = time.clock()
from inspect import stack
import shutil #module to copy files
from platform import node
import pprint # pretty printer for printing lists
import os # for getting window width and testing existence of files
import re
from re import findall # now only needed in construction.py
import sys # for getting command line input
#import glob # for testing existence of files matching a pattern
import random # for obtaining random geometry
import json
import logging # instead of the large amount of print statements not using it at the moment
from copy import deepcopy # for keeping matrices while changing others

# import my own modules
import inputreader as inr
import construction as zcon #all functions needed for constructing new geometries
import reader as r # this reads the zmatrix in gaussian format
from predictions import predictor
from montecarlo import montecarloprocedure
from loggings import loggings
# import submitter as subm
# import datareader

# import utils 
# from CINDES.utils.molecule import Molecule
from CINDES.utils.writings import log_io, print_title, sprint, dump

# initial global variables
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
once = 0
zmatrixfile = "ZMAT"
#pp = pprint.PrettyPrinter(indent=4, width=100)

print "time for imports:", time.clock()-start
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class BaseRun(object):
    ''' This is the main object for all the parameters used during any process
    this object is initiated with a dictionary from the inputreader '''
    def __init__(self, **entries):
        self.__dict__.update(entries) #here all the key/value pairs in entries are converted to attributes.
        # set system variables 
        self.script = stack()[0][1]
        self.node = node()
        self.starttime = time.time()
        self.directory = os.getcwd()
        self.pid = os.getpid()
        self.ppid = os.getppid()
        # for self.setup_filesystem one needs to have: self.(-nosub / -program)
        self.setup_filesystem()
        self.set_calcs()
        return

    def __str__(self):
        sb=['BaseRun object with the following attributes:']
        empty_attributes=[]
        for key,value in sorted(self.__dict__.items()):
            try:
                if not value: # i.e. value is either None, False, zero, empty list/string
                    empty_attributes.append(key)
                    continue
            except ValueError: pass
            if key in ['predictions']:
                sb.append("{key:20}=".format(key=key))
                sb.append( dump( value ) )
            elif key in ['TZmat','genalg', 'adj', 'jobs', 'stabjobs', 'prejobs', 'extrajobs', 'calcs']:
                sb.append("{key:20}=".format(key=key))
                sb.append( pprint.pformat(value, width=150) )
            elif key=='function' and callable(value): #i.e. the value is a lambda function
                if value.__doc__:
                    sb.append("{key:20}={value}".format(key=key, value=value.__doc__))
                else:
                    sb.append("{key:20}= lambda function")
            elif key in ['adj']:
                sb.append("{key:20}=\n".format(key=key))
                f = lambda v:''.join([('0','1')[int(item)] for item in v ])
                sb.append('\n'.join(map(f,value)))
            else:
                sb.append("{key:20}='{value}'".format(key=key, value=value))
        sb.append("empty attributes    ={}".format(" ".join(empty_attributes)))
        return '\n'.join(sb)

    def __repr__(self):
        return self.__str__()

    # the the calculationskeyword:
    def set_calcs(self):
        def tocalc(paras, job):
            # path only set after setup_filesystem!
            calc=paras[job]
            tohavekeys=['program', 'nprocs', 'identify', 'path', 'nosub']
            for key in tohavekeys:
                if not key in calc:
                    calc[key]=paras[key]
            return calc
        def check(cal, i):
            if cal['identify'] in identifiers:
                cal['identify']="{}{}_".format(cal['identify'],str(i))
                assert not cal['identify'] in identifiers
                i+=1
            identifiers.append(cal['identify'])
            return cal, i

        paras=self.__dict__
        calcs=[tocalc(paras, 'jobs')]
        if paras['prejobs']: calcs.insert(0, tocalc(paras, 'prejobs'))
        for key in ['stabjobs', 'extrajobs']:
            if paras[key]: calcs[-1] = [calcs[-1], tocalc(paras, key)]

        # verify that there are not similar identifiers
        identifiers=[]
        i=1
        for calc in calcs:
            if isinstance(calc, list):
                for cal in calc:cal, i=check(cal, i)
            else: calc, i=check(calc, i)
        #print "calcs:", pprint.pprint(calcs)
        self.calcs=calcs
        return

    def setup_filesystem(self):
        param = self.__dict__
        if self.nosub==1:
            path=''
        else:
            #param['workdir'] = os.getcwd()
            self.workdir = os.getcwd()
            path = self.workdir  + '/CALC'
            param["path"]=str(path)
            if not os.path.exists(path):
                os.makedirs(path)
            if param['program'] == 'gaussian':
                self.script='ID_gauss'
                self.extension='.com'
                shutil.copy(os.getcwd()+'/ID_gauss',path)
            elif param['program'] == 'orca':
                self.script='ID_orca'
                shutil.copy(os.getcwd() + '/ID_orca',path)
            elif param['program'] == 'nwchem':
                self.script='ID_NWChem'
                self.extension=''
                shutil.copy(os.getcwd()+'/ID_NWChem',path)
            else:
                raise SystemExit('ERROR: No valid program specified')
        self.path = path
        return param, path

    def currenttime(self):
        return "Current time %s" % str(time.time() - self.starttime)

class FrameRun(BaseRun):
    ''' This inherites from BaseRun and is the main object for all BFS/SD molecular frame based 
    procedures.
    '''
    def __init__(self,**entries):
        super(FrameRun, self).__init__(**entries)
        # specific for FrameRun:
        self.TZmat = r.geometry(**entries)
        self.adj = self.set_adj(self.TZmat['core'], self.TZmat['active'])
        self.corresp = self.set_corresp( self.TZmat['active'], self.TZmat['passive'])
        return

    def set_adj(self, core, active):
        debug=0
        from CINDES.utils.converter import Converter
        import numpy as np
        conv = Converter()
        conv.read_zmalist(core)
        xyz = np.asarray( [ atom[1] for atom in conv.zmatrix_to_cartesian() ] )
        #print xyz
        ncore = len(xyz)
        adj = np.zeros([ ncore, ncore ])
        for i in range(ncore):
            for j in range(i,ncore):
                adj[i][j]= 0.1 < np.linalg.norm( xyz[i] - xyz[j] ) < 2.0
                adj[j][i]= adj[i][j]
        sites = [ int(methyl[0][1])-1 for methyl in active ]
        sites_adj = adj[sites][:,sites]
        if debug:
            print "adjacency matrix of core:", adj
            print "self.sites:", self.sites
            print "active: ", active
            print "sites: ", sites
            print "sites_adj:", sites_adj
        return sites_adj

    def set_corresp(self, active, passive):
        '''makes a dictionary that gives the correspondance of sites with position in core matrix'''
        corresp = dict()
        for site in active:
            corresp[ int(site[0][1]) ] = int(site[1][1])
        for site in passive:
            corresp[ int(site[0][1]) ] = int(site[1][1])
        return corresp

#END CLASS RUN -----------------------------------------

# during the RUN one has to set different variables based on current state and input:
# 1. startconfiguration
# 2. geometry (zmatrices of core / active sites / passive sites
# 3. site order
# 4. table (database)
# 5. optimum at the start of the run
# 6. optimum configuration within run
# 7. global optimum
# 8. logging of output # loggings module

# 1 startconfiguration
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

# 3 site order (sequence)
def get_sequence(count, myrun):
    # START set sequence INPUT: param
    param = myrun.__dict__
    #nsites = param['nsites'] - param['nlinks']
    nsites = param['nsites']
    if 'sequences' in param:
        try:
            sequence = param['sequences'][count-1] #accounting for the fact count starts counting at 1
        except IndexError:
            sequence=random.sample(range(nsites), nsites)
        finally:
            print "SEQUENCE: " , str(sequence)
            return sequence
    if count==1 and param['restart']==4:
        sequence=param['sequence']
    else:
        if param['norandom']==1:
            sequence=range(nsites)
        elif param['sequence'] == []:
            sequence=random.sample(range(nsites),nsites)
        else:
            print "sequence read from file"
            sequence=param['sequence']
    ##output sequence
    print "SEQUENCE: " , str(sequence)
    return sequence

# 4 table (database)
def set_table(myrun, array=[]):
    '''this function loads molecules from a given database it uses a few runattributes:
        - tablename (str)
        - restart (int)
        - props (set)

    '''
    #------- enclosed function 1
    def try_oldstyle(tablename):
        import pickle
        print tablename
        with open(tablename,'rb') as f:
            pickle_db = pickle.load(f)
        print "pickled table is loaded"
        print "pickle_db:", pickle_db
        json_db = dict()
        #tableprops=['mw','solv', 'e0_solv', 'e1_solv', 'lumo', 'solv']
        #tableprops=['omega']
        tableprops=['gap', 'lumo', 'homo']
        for item in pickle_db:
            key=item[0]
            value={prop:prop_value for prop,prop_value in zip(tableprops,item[1:])}
            json_db[key]=value
        print "an old_style formatted tablefile was loaded with props:", tableprops
        # touch new json file
        with open('{}.json'.format(tablename),'w') as f2:
            json.dump(json_db, f2, indent=-1)
        raise SystemExit('stop')
        return json_db
    #-------- enclosed function 2:
    def adjust_dihedrals(table, array):
        # 1. create a dictionary for each site that maps the group to group-dihedral
        #print "array:", array
        specific_dihedrals=[{} for _ in array]
        for groupssite, groupsdihedrals in zip(array, specific_dihedrals):
            for group in groupssite:
                if is_float(group[-1]):
                    dgroup=''.join(group)
                    group=''.join(group[:-1])
                    groupsdihedrals[group]=dgroup
        print "specificdihedrals:", specific_dihedrals
        ###
        # 2. convert each index in table to the correct dihedral 
        print table.keys()[:20]
        new_table={}
        for key in table:
            conf = key.split('_')
            new_conf=[]
            for group, groupsdihedrals in zip(conf, specific_dihedrals):
                # 1. remove any dihedrals from group
                try:
                    group=group.translate(None, '0123456789')
                except TypeError:
                    group=group.translate({ord(ch): None for ch in '0123456789'})
                dgroup=groupsdihedrals.get(group,group)
                new_conf.append(dgroup)
            new_key='_'.join(new_conf)
            new_table[new_key]=table[key]
        print new_table.keys()[:20]
        #raise SystemExit
        return new_table

    #-------- enclosed function 3
    def remove_dihedrals(table):
        # convert each index in table to the correct dihedral 
        new_table={}
        i=0
        for key in table:
            conf = key.split('_')
            new_conf=[]
            for group in conf:
                # 1. remove any dihedrals from group
                try:
                    group=group.translate(None, '0123456789')
                except TypeError:
                    group=group.translate({ord(ch): None for ch in '0123456789'})
                new_conf.append(group)
            new_key='_'.join(new_conf)
            new_table[new_key]=table[key]
        return new_table
    # -------
    print "------------"
    # try if tablename is given
    try:
        tablename = myrun.tablename
    except AttributeError:
        tablename = 'table.json'
    # look if extension is used otherwise set it automatically
    if not tablename[-5:]=='.json': tablename='{}.json'.format(tablename)

    if myrun.restart>0:
        # look if tablename is given in INPUT otherwise default
        try:
            with open(tablename,'rb') as f:
                db = json.load(f)
        except (IOError,ValueError):
            print "no json table"
            print "try to load as pickle {}".format(tablename[:-5])
            try:
                db = try_oldstyle(tablename[:-5])
            except IOError:
                print "also no correct pickled table"
                raise
        print "loaded json database with {} molecules".format(len(db))
        # myrun.props has to be a subset of value.viewkeys(): set operations <= means "is subset of"
        table = { key:value for key,value in db.iteritems() if myrun.props <= value.viewkeys() }
        print "made a table with {} molecules that have the required properties".format(len(table))

        # should the function value be included in the table? otherwise here is the place ;)
        if myrun.adjust_dihedrals:
            table = adjust_dihedrals(table, array)
        else:
            table = remove_dihedrals(table)

    else:
        table = dict()
        open(tablename,'wb').close()
    return table

# 5 optimum at the start of the run
def set_optimum(myrun,table):
    if myrun.restart >= 3:
        raise SystemExit('deprecated functionality')
        #tabledict = dict( [ item[0:2] for item in table ] )
        optimum = molecule(conf=myrun.startconf)
        optimum.props = table[molecule.index]
        # set Pvalue?
        #optimum = [ myrun.startconf, tabledict[myrun.startconf] ]
        logging.info("optimum:"+ str(optimum))
    else:
        optimum = None
    return optimum

# 6 optimum within the global iterations
def testmax(myrun, mols, bcok):
    ''' sets optsite
    multiple boundary conditions are not yet implemented

    '''
    param = myrun.__dict__
    # there are molecules that are predicted and are not calculated so they have molecule.Pvalue is None? 
    # NOT: no they should have a Pvalue but just there predicted attribute is set to True
    # so which molecule could possibly have a None Pvalue? 
    # neglected molecules will have a None value so indeed filter them out 
    # data = [ molecule.log() for molecule in mols if ( molecule.predicted == False and not molecule.Pvalue is None) ]
    # data = [ molecule.log() for molecule in mols if not ( molecule.predicted == False or not molecule.Pvalue is None) ]
    # INDEED mols that are only predicted are left out! Pred values are only used at the decision for tocal/nocal!
    mols = [ mol for mol in mols if mol.predicted == False and not mol.Pvalue is None ]

    if 'bcprop' in param:
        if param['bcoptimum'] in ['min','Min','MIN']:
        #test if BC fullfilled. 
            try:
                satisfactory = [ mol for mol in mols if mol.boundaries[0] < float(param['bcval']) ]
            except TypeError as e:
                print e
                pass
        else:
            assert param['bcoptimum'] in ['max','Max','MAX']
            try:
                satisfactory = [ mol for mol in mols if mol.boundaries[0] > float(param['bcval']) ]
            except TypeError as e:
                print e
                pass

        print "satisfactory:\n", pprint.pformat(satisfactory, width=100)
        print "the BC condition is bc<:", param['bcval']
        if satisfactory==[]: #so if there is at least one fullfilling BC
            bcok=0
            print "BC not fullfilled:"
            #optsite = min(data,key = lambda x:x[2])
            optsite = min(mols,key = lambda mol:abs( mol.boundaries[0] - float(param['bcval']) ) )
        else: #BC not yet
            bcok=1
            print "BC fullfilled; satisfactory is not empty:", pprint.pformat(satisfactory, width=100)
            #optsite = max(satisfactory,key = lambda x:x[1])
            if param['optimum'] in ['minimum', 'min']:
                optsite = min(satisfactory,key = lambda x:x[2])
            else:
                optsite = max(satisfactory,key = lambda x:x[2])
    else:
        bcok=1 #no BC but need this variable to test later on
        # optimum of the list or MINIMUM
        if param['optimum'] in ['minimum', 'min']:
            if param['cutoff'] == 0:
                optsite = min(mols,key = lambda x:x.Pvalue)
            else:
                testmols = [ mol for mol in mols if abs(mol.Pvalue) > param['cutoff'] ]
                optsite = min(testmols, key = lambda mol:mol.Pvalue)
        else:
            optsite = max(mols, key = lambda mol:mol.Pvalue)

    optsite.opt=True
    #logging.warning('optsite:' + pprint.pformat(optsite, width=100))
    return optsite, bcok

# 7 set global optimum and define convergence and redirect to Monte Carlo component
def runtest(run, optimum, optsite, count, bcok,mctable=[], array=[]):
    param = run.__dict__
    TZmat= run.TZmat
    converged=0
    print

    #raise SystemExit('optimum and optsite should be Molecule instances now')
    if (count > 1 and bcok) or param['restart']>=3: #BCOK is a test of the boundary condition is already fullfilled
        if optimum==optsite:  #test the property value! not 1 anymore!
            print "optimum is the same!"
            print "converged to a optimum configuration!"
            if param['montecarlo'] == 0:
               converged = 1
            else:
               if param['ml']==0:
                   optsite = montecarloprocedure(run, array, optimum, mctable)
               else:
                   optsite = montecarloprocedure(run, array, optimum, mctable, **TZmat)
               print "optimal_after_this_site:", pprint.pformat( optsite, width=100 )
        else:
            print "Global_Iteration_optimum and optimum_after_this_site are not the same yet"
            print "gi_optimum:" ,optimum
            print "current optimum:" ,optsite
    else: #except NameError:
        print "NameError no optimal structure or BC not yet fullfilled."
        #pass
    optimum = optsite.copy()
    return optimum, optsite, converged

def get_property_table(table, myrun):
    '''set a dict with {'index1':prop1, etc. } to use for montecarlo and prediction making '''
    db=dict()
    for key,value in table.iteritems():
        if myrun.property=='func':
            kwargs = { prop:value[prop] for prop in myrun.func_args }
            Pvalue = myrun.function(**kwargs)
            db[key]=Pvalue
        else:
            Pvalue = value[ myrun.property ]
            db[key]=Pvalue
    return db

# DATA GETTING:
# A: fake data for testing (skipper)
def skipper(mols_tocal,mols_nocal,iprint=True):
    ''' generate random data '''
    if iprint: print "submit is skipped! random data is generated"
    import string
    for molecule in mols_tocal:
        item = molecule.index
        #propx= sum([ string.uppercase.index(itempje)+1 for itempje in list(item.replace('_',''))]) 
        output = 0
        replaced = item.replace('_','')
        replaced = filter(lambda x:x.isalpha(), replaced)
        for i in replaced:
            try:
                output += string.uppercase.index(i)
            except ValueError:
                output += string.lowercase.index(i)
        propx = float(output)
        try:
            if 'bcprop' in param:
                propy= len(item.replace('_',''))
                molecule.boundaries = [ float(propy) ]
        except NameError:
            pass

        molecule.Pvalue = propx
        molecule.predicted = False

    mols_all = mols_tocal + mols_nocal
    return mols_all

def restriction1(mols_todo, mols_nodo, run):
    ''' test if not B-B A or N-N bond present in molecules '''
    def has_forbidden_combination(conf, adj):
        for i in range(len(conf)):
            for j in range(i,len(conf)):
                if conf[i]==conf[j] and adj[i][j]==1.0 and conf[i] in [ ['N'], ['B'] ]:
                    print "forbidden combination: ", i, conf[i], j, conf[j], adj[i]
                    return True
        return False

    from itertools import combinations, ifilterfalse

    print "nmol:", len(mols_todo)

    for molecule in mols_todo:
        print "mol.conf:", molecule.conf
        if has_forbidden_combination(molecule.conf, run.adj):
            mols_todo.remove(molecule)


    print "nmol:", len(mols_todo)

    return mols_todo, mols_nodo





# B: getting the real data by submitting 
def submittingprocedure(mols_tocal,mols_nocal,myrun,**kwargs):
    fileparameters = myrun.__dict__
    import calculator
    data = calculator.procedure(myrun,mols_tocal,mols_nocal,kwargs)
    return data

# THERE ARE DIFFERENT GLOBAL PROGRAM FLOW PROCEDURES:
# 1: STANDARD PROCEDURE: Best First Search: BFS()
# 2: Generate 1 Configuration input file: genconf
# 3: Generate total chemical space defined by the sites and functionalisations: generate
# 4: Generate a number of random structures and print them to screen: genrandom
# 5: A testrun. Not implemented. a helper function for the test functions in ./tests/tests.py: testrun
# 6: Steepest Descent algorithm. Looks like BFS but there is no loop over sites
# 7: Generate database based on farthest point selection. (based on diversity index)

# 1: standard BFS
def BFS(param,array):
    bcok=0 #TO REMOVE LATER
    param['bcok']=0
    startconf = get_startconf(param,array)

    #SET MYRUN CLASS and assign all necessary attributes
    myrun = FrameRun(**param)
    print(myrun) #this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun, array)
    property_table = get_property_table(table, myrun)
    # set optimum
    optimum = set_optimum(myrun,table)
    #set calculation properties
    #END MYRUN CLASS assignments. from now myrun should contain all the necessary information to work with during the whole program run.


            # ------------------------------------- # 
            # --- HERE THE MAIN LOOP STARTS --- --- #
            # ------------------------------------- #
    count = 1 # so we start counting at 1!
    while True:
        print_title("COUNT: " + str(count),outline='l',signator="-")

        ### set site order in sequence INPUT: param, count
        sequence = get_sequence(count, myrun)

        ### for each site in sequence:
        for l in range(len(sequence)):
            k = sequence[l]
            print_title("k(site)= {} l(nsite)= {} (c={})".format(k,l,count),outline='l',signator='=')
            if not l == 0 or count > 1: #define new startconfiguration if not first cycle
                # define new starting geometry
                print "optsite:",optsite
                del startconf
                startconf = zcon.indtocon(optsite.index)

            # STEP 1: INDEXMAKER
            #get indices_all and the indices that still need to be calculated
            # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
            #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table )
            mols_todo, mols_nodo = zcon.classmaker2(startconf,array,k,table, myrun )
            if 1 in myrun.restrictions: # this are actually filters!
                mols_todo, mols_nodo = restriction1(mols_todo, mols_nodo, myrun )
            print "|      NEW POPULATION CONSTRUCTED:"
            print "|   mols_todo:"
            if mols_todo:
                for mol in mols_todo: print "|      {}".format(mol)
            else: print "|      -"
            print "|   mols_nodo:"
            if mols_nodo:
                for mol in mols_nodo: print "|      {}".format(mol)
            else: print "|      -"

            # STEP 2: PREDICTOR
            # perform prescreaning in a predictions.
            mols_nocal, mols_tocal, made_pred = predictor(
                    myrun,
                    property_table,
                    mols_todo,mols_nodo,
                    count,
                    array=array,
                    nsite=l
                    )

            # STEP 3: SUBMITTING PART
            if not myrun.nosub==1:
                mols_all = submittingprocedure(mols_tocal,
                                               mols_nocal,
                                               myrun,
                                             **myrun.TZmat     ) # here call submitting procedure
            else: mols_all = skipper(mols_tocal,mols_nocal)

            # STEP 4: UPDATE OPTIMUM STRUCTURE
            # decide what the optimum site is and if the bc if fullfilled
            print "BCOK:", bcok
            optsite, bcok = testmax(myrun, mols_all, bcok)

            # STEP 5: UPDATE DATABASE and LOG results of microiteration
            # logs new elements in data to table and tablebin and whole data to cyclesinfo
            table = loggings(mols_all,
                    table,
                    count,
                    k,l,
                    made_pred,
                    tablename = myrun.tablename)
            property_table = get_property_table(table, myrun)

            print("--- %s seconds ---" % (time.time() - myrun.starttime))
            print(myrun.currenttime())
        # END LOOP OVER SITES

        #get optimum and test convergence
        optimum, optsite,converged = runtest(myrun, optimum, optsite, count, bcok, mctable=property_table, array = array)
        if converged==1: break
        count +=1
        if count > param['maxiter']:
            print "maxiterations is reached"
            print "optimum is: ", optimum
            break
    # ---------------------------- #
    # ------ END OF LOOPING ------ # 
    # ---------------------------- #
    print "DONE"
    return

# 2: genconf
def genconf(param):
    myrun = FrameRun(**param)
    param = myrun.__dict__
    TZmat = r.geometry(**param)
    conf = zcon.indtocon(param['startind'])
    print "in GENCONF: conf is:", conf
    param['workdir'] = os.getcwd()
    path = param['workdir']
    param["path"]=str(path)

    c = deepcopy(TZmat['core'])
    a = deepcopy(TZmat['active'])
    p = deepcopy(TZmat['passive'])
    mat = zcon.constructor2(conf,c,a,p, links=myrun.symlinks)
    if myrun.program=='gaussian':
        import gaussian as program
    elif myrun.program=='nwchem':
        import nwchem as program
    else:
        raise SystemExit('program not recognized')
    from CINDES.utils.molecule import Molecule
    mol = Molecule(index=param['startind'])
    mol.zmat = mat
    program.filewriter(mol, **param)
    return

# 3: generate
def generate_procedure(param,array):
    ''' calculate all possible structures '''

    myrun = FrameRun(**param)
    table = set_table(myrun)
    print myrun

    # get all structures
    print "len table:", len(table)

    mols = get_all_molecules(array)
    #1b check already in database
    mols_todo, mols_nodo = zcon.check_in_table(mols, table, myrun)
    #1c eventueel predictions
    mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, 0, array=array)

    # calculate them
    #print_title("COUNT: " + str(count),outline='l',signator="-")

    count=0
    if False:
        if not myrun.nosub==1:
            mols_all = submittingprocedure(mols_tocal,
                                           mols_nocal,
                                           myrun,
                                         **myrun.TZmat     ) # here call submitting procedure
        else: mols_all = skipper(mols_tocal,mols_nocal)
    elif False:
        # use batches
        batchsize=49
        def chunks(l,n):
            '''yields successive n-sized chunks of l'''
            for i in range(0, len(l), n):
                yield l[i:i+n]

        mols_all=[]
        for i,batch in enumerate(chunks(mols_tocal,batchsize)):
            print "chunk nr:", i, "with ", len(batch), "structures"
            if not myrun.nosub==1:
                batch = submittingprocedure(batch,
                                               [],
                                               myrun,
                                             **myrun.TZmat     ) # here call submitting procedure
            else: batch = skipper(batch,[])
            mols_all.extend(batch)

    elif True:
        # use job arrays. 
        # make a jobscript with the line:
        # qsub -t 1-njobs
        # $PBS_ARRAY_INDEX has to be used in the submitscript to submit each job. 
        # each number has to refer to a certain filename.
        # make file with all the filenames:
        with open('./CALC/filenames.txt','w') as fout:
            fout.write( '\n'.join([myrun.identify + mol.index for mol in mols_tocal]) )
        print "njobs:", len(mols_tocal)

        # generate all inputfiles
        from calculator import filemaker, geommaker
        geommaker(mols_tocal,myrun,**myrun.TZmat)
        filemaker(mols_tocal,myrun) #----------------------------------HERE IS THE FILEWRITER CALL


    print "mols_all:", mols_all
    table = loggings(mols_all,table,count,1,1, made_pred=made_pred )

    print "DONE"
    return

def get_all_molecules(array):
    from CINDES.utils.molecule import Molecule
    print "in get_all_molecules"
    #A = [ map(''.join,item) for item in array ]
    A = array
    C=[[]]
    for site in A:
        D=[]
        for item in C:
            for group in site:
                conf = item+[group]
                D.append(conf)
        C=D
    print "molecules:"
    for i, item in enumerate(C): print i, item
    mols = [ Molecule(conf=conf, dihedral=True) for conf in C ]

    print mols
    return mols

def generate2(core,active,passive,converter,**kwargs):
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
            mat = zcon.constructor2(confs[i],c,a,p, **kwargs)
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

# 4: genrandom
def genrandom(param,array):
    ''' generate random structures and print '''
    #set restart and startind to default:
    print
    for _ in xrange(param['nrandom']):
            conf = []
            for i in range(len(array)):
                group = random.choice(array[i])
                conf.append( group )
            print zcon.contoind(conf)
    print
    return True

# 5: testrun
def testrun(param,array):
    pass

# 5b: testpred
def testpred(param,array):
    ''' run the predictions on the tablebin file '''
    print_title("Testing Prediction procedure activated!", outline='l', signator=':')
    class Mol(object):
        def __init__(self):
            self.predictions = {}
        def __repr__(self): return "<empty molecule object>"
    myrun = FrameRun(**param)
    print myrun
    table = set_table(myrun, datacolumn=param['datacolumn'])
    sprint(10,table)
    mols_todo, mols_nodo = ([Mol(),],[Mol(),])
    mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, 99, array=array, nsite=0)
    return

# 6: steepest descent
def SteepestDescent(param,array):
    bcok=0 #TO REMOVE LATER
    param['bcok']=0

    #SET MYRUN CLASS and assign all necessary attributes
    myrun = FrameRun(**param)
    print(myrun) #this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun)
    # set optimum
    optimum = set_optimum(myrun,table)
    #set calculation properties
    startconf = get_startconf(param,array)
    #END MYRUN CLASS assignments. from now myrun should contain all the necessary information to work with during the whole program run.

    myrun.restingsites = range( myrun.nsites ) # defines which sites will be changed. only relevant for steepest2 algorithm

            # ------------------------------------- # 
            # --- HERE THE MAIN LOOP STARTS --- --- #
            # ------------------------------------- #
    count = 1 # so we start counting at 1!
    while True:
        print_title("COUNT: " + str(count),outline='l',signator="-")

        if count > 1: #define new startconfiguration if not first cycle
            # define new starting geometry
            print "optsite[0]",optsite[0]
            del startconf
            startconf = optsite.conf
            #startconf = zcon.indtocon(optsite[0])
            print "newconf: ", startconf

        # STEP 1: INDEXMAKER
        #get indices_all and the indices that still need to be calculated
        #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker_SD(startconf, array, table, myrun )
        mols_todo, mols_nodo = zcon.classmaker2_SD(startconf,array,table, myrun )

        # STEP 2: PREDICTOR
        # perform prescreaning in a predictions. 
        mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, count, array=array)
        #data_nocal,indices_tocal, predict = predictor(myrun, table, indices_todo,data_nodo, count, array=array)

        # STEP 3: SUBMITTING PART
        if not myrun.nosub==1:
            mols_all = submittingprocedure(mols_tocal,
                                           mols_nocal,
                                           myrun,
                                         **myrun.TZmat     ) # here call submitting procedure
        else: mols_all = skipper(mols_tocal,mols_nocal)

        # STEP 5: UPDATE DATABASE and LOG results of microiteration
        # logs new elements in data to table and tablebin and whole data to cyclesinfo
        table = loggings(mols_all,table,count,1,1, made_pred=made_pred )

        # STEP 6: UPDATE OPTIMUM STRUCTURE
        # decide what the optimum site is and if the bc if fullfilled
        print "BCOK:", bcok
        optsite, bcok = testmax(myrun, mols_all, bcok)

        if myrun.procedure=='steepest2':
            maxconf = zcon.indtocon(optsite[0])
            print "maxconf:", maxconf, 'while startconf:', startconf
            try:
                changedsite = [ siteM == siteS for siteM,siteS in zip(maxconf,startconf) ].index(False)
            except ValueError:
                print "no site changed"
            else:
                myrun.restingsites.remove(changedsite)
            if myrun.restingsites==[]:
                print "all sites changed once"
                break


        print("--- %s seconds ---" % (time.time() - myrun.starttime))
        print(myrun.currenttime())
        # END LOOP OVER SITES

        #get optimum and test convergence
        optimum, optsite,converged = runtest(myrun, optimum, optsite, count, bcok, mctable=table, array = array)
        if converged==1: break
        count +=1
        if count > param['maxiter']:
            print "maxiterations is reached"
            print "optimum is: ", optimum
            break
    # ---------------------------- #
    # ------ END OF LOOPING ------ # 
    # ---------------------------- #
    print "DONE"
    return

