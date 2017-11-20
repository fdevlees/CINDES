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
# from CINDES4.utils.molecule import Molecule
from CINDES4.utils.writings import log_io, print_title, sprint, dump

# initial global variables
logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
once = 0
zmatrixfile = "ZMAT"
#pp = pprint.PrettyPrinter(indent=4, width=100)

print "time for imports:", time.clock()-start


class Run(object):
    ''' This is the main object for all the parameters used during the process
    this object is initiated with a dictionary from the inputreader '''
    def __init__(self,**entries):
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

        #zmatrix reading and splitting needs: self.-ncore / -line1 / -nch3
        #self.TZmat = r.geometry(param)
        self.TZmat = r.geometry(**entries)

        self.adj = self.set_adj(self.TZmat['core'], self.TZmat['active'])
        self.corresp = self.set_corresp( self.TZmat['active'], self.TZmat['passive'])

        #sets Gaussian09 input lines
        self.set_calculation_properties()
        return

    def __str__(self):
        sb=['Run object with the following attributes:']
        for key,value in sorted(self.__dict__.items()):
            if key in ['predictions']:
                sb.append("{key:20}=".format(key=key))
                sb.append( dump( value ) )
            elif key in ['TZmat','genalg', 'adj', 'jobs']:
                sb.append("{key:20}=".format(key=key))
                sb.append( pprint.pformat(value, width=150) )
            else:
                sb.append("{key:20}='{value}'".format(key=key, value=value))
        return '\n'.join(sb)

    def __repr__(self):
        return self.__str__()

    def currenttime(self):
        return "Current time %s" % str(time.time() - self.starttime)

    @log_io()
    def setup_filesystem(self):
        param = self.__dict__
        if self.nosub==1:
            path=''
        else:
            #param['workdir'] = os.getcwd()
            self.workdir = os.getcwd()
            path = self.workdir  + '/CALC'
            param["path"]=str(path)
            logging.info("PATH:"+str(path))
            if not os.path.exists(path):
                os.makedirs(path)
            if param['program'] == 'gaussian':
                shutil.copy(os.getcwd()+'/ID_gauss',path)
            elif param['program'] == 'orca':
                shutil.copy(os.getcwd() + '/ID_orca',path)
            else:
                raise SystemExit('ERROR: No valid program specified')
        self.path = path
        return param, path

    def set_calculation_properties(self):
        if self.program in ['Gaussian','gaussian']:
            self.runspecs_gaussian()
        elif self.program in ['ORCA','orca','Orca']:
            runspecs_orca(param)
        else:
            raise SystemExit('PROGRAM NOT RECOGNIZED')
        return

    def set_adj(self, core, active):
        debug=0
        from CINDES4.utils.converter import Converter
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
            corresp[ site[0][1] ] = site[1][1]
        for site in passive:
            corresp[ site[0][1] ] = site[1][1]
        return corresp

    def runspecs_gaussian(self):
        param=self.__dict__
        if param['gaussianlines']:
            pass
        elif param['stab']==1:
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
                self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '#p geom=allcheck guess=read polar volume=tight '+param['functional']+'/'+param['basisset']+'\n'
            else:
                self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '#p geom=allcheck guess=read polar '+param['functional']+'/'+param['basisset']+'\n'
        elif param['aip']==1 or param['aea']==1:
            self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['twojob']=1
            self.gaussianline2 = '# geom=check guess=read opt scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
        elif param['ip']==1 or param['ea']==1:
            self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['twojob']=1
            self.gaussianline2 = '# geom=check guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
        else:
            if param['property']=='dipole':
                logging.warning('NO Geometry optimization will be performed!!!')
                self.gaussianline = '# ' + param['functional'] +'/'+ param['basisset'] +'\n'
            else:
                if param['basisset'] in [ None, 0, '0', 'none', 'nalse', False, 'off' , 'n', 'na' ]:
      
                    self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            if param['twojob'] == 1:
                self.gaussianline2 = '# geom=check guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            elif param['twojob'] == 2:
                self.multiplejobs = 2
                if param['semiempirical'] == 1:
                    self.gaussianline  = '# opt=(maxcycle=' + param['maxcycles'] + ') ' + 'pm6' +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '# geom=allcheck guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            elif param['twojob'] == 4:
                self.multiplejobs = 4
                if param['semiempirical'] == 1:
                    self.gaussianline  = '# opt=(maxcycle=' + param['maxcycles'] + ') ' + 'pm6' +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=' + param['maxcycles'] + ') scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '# geom=check scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                if param['solv']:
                    if True:
                        self.gaussianline_solv0 = '# geom=allcheck guess=read pm6\n'
                        self.gaussianline_solv1 = '# geom=allcheck guess=read pm6 scrf=(smd, solvent=aceticacid)\n'
                    else:
                        self.gaussianline_solv0 = '# geom=allcheck guess=read scf=xqc b3lyp/6-31G(d,p)\n'
                        self.gaussianline_solv1 = '# geom=allcheck scf=xqc scrf=(smd, solvent=aceticacid) b3lyp/6-31G(d,p)\n'

        for key in ['ip','ea','polar']:
            if param[key]==1:
                self.multiplejobs +=1
        return
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
def set_table(myrun):
    '''this function loads molecules from a given database it uses a few runattributes:
        - tablename (str)
        - restart (int)
        - props (set)

    '''
    def try_oldstyle(tablename):
        import pickle
        print tablename
        with open(tablename,'rb') as f:
            pickle_db = pickle.load(f)
        print "pickled table is loaded"
        print "pickle_db:", pickle_db
        json_db = dict()
        #tableprops=['mw','solv', 'e0_solv', 'e1_solv', 'lumo', 'solv']
        tableprops=['omega']
        for item in pickle_db:
            key=item[0]
            value={prop:prop_value for prop,prop_value in zip(tableprops,item[1:])}
            json_db[key]=value
        print "an old_style formatted tablefile was loaded with props:", tableprops
        # touch new json file
        open('{}.json'.format(tablename),'w').close()
        return json_db

    #--------
    print "------------"
    # look if tablename given in INPUT otherwise default
    try:
        tablename = myrun.tablename
    except AttributeError:
        tablename = 'table.json'
    # look if extension is used otherwise set it automatically
    if not tablename[-5:]=='.json': tablename='{}.json'.format(tablename)

    if myrun.restart>0:
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
            optsite = max(mols, key = lambda mol:mols.Pvalue)

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
    global once
    # here submitting thing knows at least the path
    fileparameters = myrun.__dict__
    if myrun.program in ['ORCA','orca','Orca']:
        import orcafunctions
        data = orcafunctions.submittingprocedure(confs,mols_tocal,mols_nocal,fileparameters,**self.TZmat)
    elif myrun.program in ['Gaussian','gaussian']:
        import gaussianfunctions as gausf
        data = gausf.procedure(myrun,mols_tocal,mols_nocal,kwargs)
    elif myrun.program == 'molpro':
        raise SystemExit('molpro not implemented')
    else:
        print "program not recognized!:", myrun.program
        raise SystemExit('no program recognized')
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
    myrun = Run(**param)
    print(myrun) #this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun)
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
            print_title("k(site)= " + str(k) + " l(nsite)= "+ str(l),outline='l',signator='=')
            if not l == 0 or count > 1: #define new startconfiguration if not first cycle
                # define new starting geometry
                print "optsite:",optsite
                del startconf
                startconf = optsite.conf

            # STEP 1: INDEXMAKER
            #get indices_all and the indices that still need to be calculated
            # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
            #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table )
            mols_todo, mols_nodo = zcon.classmaker2(startconf,array,k,table, myrun )
            if 1 in myrun.restrictions:
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
    myrun = Run(**param)
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
    zcon.filewriter2(mat,param['startind'],**param)
    return

# 3: generate
def generate_procedure(param,array):
    ''' calculate all possible structures '''

    myrun = Run(**param)
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
    if True:
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
        from gaussianfunctions import filemaker, geommaker
        geommaker(mols_tocal,myrun,**myrun.TZmat)
        filemaker(mols_tocal,myrun) #----------------------------------HERE IS THE FILEWRITER CALL


    print "mols_all:", mols_all
    table = loggings(mols_all,table,count,1,1, made_pred=made_pred )

    print "DONE"
    return

def get_all_molecules(array):
    from CINDES4.utils.molecule import Molecule
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
    mols = [ Molecule(conf=conf) for conf in C ]

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
    myrun = Run(**param)
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
    myrun = Run(**param)
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
            startconf = zcon.indtocon(optsite[0])
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

# procedure 7. Farthest Point Selection based on the diversity index. 
def database_construction(param,array):
    debug=True
    myrun = Run(**param)
    myrun.divers_discardCH=True
    print(myrun) #this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun)
    nmax = myrun.divers_nmax

    if not table:
        # run startind calculation
        table = None #make initial calculation. get table with length 1.

    count=1
    maxiter=10
    while True: #later while True
        print_title("COUNT: " + str(count),outline='l',signator="-")
        print "len table:", len(table)
        #1. get new structure(s) to calculate
        #if myrun.divindex==1:
            #1.1. run diversity on table and get occupancy per site. 
            #1.2. select for each site the least occuring group.
        mols = getdivers(array,table, myrun)
        if debug: print mols
        #1b check already in database
        mols_todo, mols_nodo = zcon.check_in_table(mols, table, myrun)
        #1c eventueel predictions
        mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, count, array=array)
        #2. run new structure
        mols_all = submittingprocedure(mols_tocal,
                                       mols_nocal,
                                       myrun,
                                     **myrun.TZmat     ) # here call submitting procedure
        #3. add structure to table
        if debug: print "after calculation:"
        if debug: print mols_all
        table = loggings(mols_all,table,count,1,1, made_pred=made_pred )
        if debug:
            for item in table: print item
        #4. stop if maxstructures is obtained. or other convergence criteria is met. 
        count += 1
        if len(table)>= nmax:
            print "desired number of samples reached!"
            break
        if count>maxiter:
            print "maxiterations reached"
            break
    return

def getdivers(array, table, myrun):
    '''give the next n most divers molecules

    This function should be replaced later to ... ?
    '''
    #0.
    index = myrun.divers_divindex
    batchsize = myrun.divers_batchsize
    discardCH = myrun.divers_discardCH

    #1.
    from CINDES4.utils.molecule import Molecule
    from CINDES4.utils.diversity1 import Diversifier
    from CINDES4.utils.table import Tablebin
    diversifier = Diversifier(index=index)
    mols = []
    #confs = [ item[0].split('_') for item in table ]
    confs = [ item.conf for item in table ]
    print "confs:", confs[:5]
    #seq = list(set([ group for conf in confs for group in conf ]))
    seq = [ 'CH','B','O','S','N','P',
                               'CNHH', 'CNOO','COH','CSH','CPh','CCHO','CSOOOH']
    print "seq:", seq
    for _ in range(batchsize):
        if index==1 or index==2:
            occupancy_sum, occupancy = diversifier.get_occupancy12(seq, confs, discardCH=discardCH)
            print "in get divers: occupancy:", occupancy
            conf = make_molecule12(occupancy, seq, array, confs, index=myrun.divers_divindex)
        elif index==3:
            occupancy, occupancy_percentages = diversifier.get_occupancy3(seq, confs)
            #print "in get divers: occupancy:", occupancy
            conf = make_molecule3(occupancy, seq, array, confs)
        if discardCH:
            from numpy.random import binomial, shuffle, seed
            seed(40)
            # adjust conf and place CH groups in it via a binomial distribution
            nsites=10
            while True:
                #nch = binomial(nsites+3,0.5) ###### here tuning factor. 
                nch = binomial(nsites,0.5) ###### here tuning factor. 
                if nch<nsites: break
            positions = range(nsites)
            shuffle(positions)
            conf = [ item if i<nch else 'CH' for item,i in zip(conf,positions)]
        fconf = zcon.indtocon( '_'.join(conf))
        mols.append(Molecule(conf=fconf))

        # now append to conf so div values can change.
        confs.append(conf)
    print "molecules:", mols
    return mols

def make_molecule12(occupancy, seq, array, confs, index=1):
    iarray = [ map("".join,item) for item in array ]
    if index==1:
        new_conf=[]
        for site_occ, site_array in zip(occupancy,iarray):
            #most_divers_group = min(zip(site_occ,seq))[1]
            # search for group that has the least occurance, but not zero because that indicates that is cannot occur on that site
            #most_divers_group = min(filter(lambda x:not x[0]==0.,zip(site_occ,seq)))[1]
            # see if in array for that site. 
            most_divers_group = min(filter(lambda x:x[1] in site_array,zip(site_occ,seq)))[1]

            new_conf.append(most_divers_group)
        if new_conf in confs:
            print "most divers mol already in table:"
            new_conf = take_nth(occupancy, seq, array, confs)
    return new_conf

@log_io()
def make_molecule3(occupancy, seq, array, confs):
    import numpy as np
    np.random.seed(40)
    #print seq
    from collections import OrderedDict
    occD = dict()
    for group1, occ1 in zip(seq,occupancy):
        #print group1, occ1
        if group1 in ['S','O','CO']: continue #these are not possible so do not include
        for group2, occ12 in zip(seq,occ1):
            key = '{}_{}'.format(group1,group2)
            occD[key]=occ12
    occD = OrderedDict( sorted( occD.iteritems(), key=lambda x:x[1] ) )
    #print occD
    new_conf = ['']*len(confs[0])
    keys_visited=[]
    for ibond, gbond in zip(( (0,7),(1,5),(2,6),(3,9) ), occD.items()[:4]):
        keys_visited.append(gbond[0])
        i3,i2 = ibond
        g3,g2 = gbond[0].split('_')
        new_conf[i3]=g3
        new_conf[i2]=g2
    #now still two places have to be filled.
    # for bond4:
    # bond4 neighbors 1 and 3 so the tert group has to be g1 or g3
    #both_bonds_done=[False,False]
    #for i,gbond in enumerate(occD.items()[4:]):
    #    if gbond[0] in keys_visited:
    #        print "already visited:", gbond[0]
    #        continue
    #    g3,g2= gbond[0].split('_')
    #    if g3==new_conf[1] or g3==new_conf[3]:
    #        new_conf[4]=g2
    #        both_bonds_done[0]=True
    #        keys_visited.append(gbond[0]) # i think this line is not necessary
    #    elif g3==new_conf[2] or g3==new_conf[3]:
    #        new_conf[8]=g2
    #        both_bonds_done[1]=True
    #        keys_visited.append(gbond[0]) # i think this line is not necessary
    #    #else:
    #    #    keys_visited.append(gbond[0]) # i think this line is not necessary
    #    #    continue
    #    print "i, both_bonds_done:", i, both_bonds_done
    #    keys_visited.append(gbond[0]) # i think this line is not necessary
    #    if all(both_bonds_done):
    #        break
    #else:
    #    raise StandardError
    for gbond in occD.items():
        if gbond[0] in keys_visited:continue
        g3,g2= gbond[0].split('_')
        if g3==new_conf[1] or g3==new_conf[3]:
            new_conf[4]=g2
            break
    else:
        raise StandardError
    #for bond8
    # bond 8 neighbors 2 and 3 so the tert group has to be g2 or g3
    for gbond in occD.items():
        if gbond[0] in keys_visited:continue
        g3,g2= gbond[0].split('_')
        if g3==new_conf[2] or g3==new_conf[3]:
            new_conf[8]=g2
            break
    else:
        raise StandardError
    print new_conf
    assert not '' in new_conf, "one group not defined!"
    # now the change that a group appears on 4 8 is different from appearing on the others? so randomly symmetry permutation:
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
    # randomly select one item from list
    new_order = Adasym[np.random.choice(range(10))]
    # permute new_conf according to new_order
    randomized = np.array(new_conf)[np.array(new_order)-1]
    return list(randomized)

def take_nth(occupancy, seq, array, confs):
    ''' it is possible that the previous function returns a molecule that already exists in the database
    so here i'll write a clever method to come up with a nth-but-most divers structure. '''
    iarray = [ map("".join,item) for item in array ]
    nsites = len(iarray)

    # the next line combines the seq and occupancy and filters only the ones that also are allowed on that site. i.e. they are in array
    sel_occseq = [ filter(lambda x:x[1] in iarray_site,zip(occ_site, seq)) for occ_site, iarray_site in zip(occupancy, iarray) ]
    # sort each element:
    sorted_occseq = [ sorted(item) for item in sel_occseq ]

    #make a list of possibly one-but-lasts. 
    nthbutbestconfs = []

    # loop over sites
    identity = [[0]*n + [1] + [0]*(nsites-n-1) for n in range(nsites) ]
    print "identity:", identity
    for row,sorted_occseq_site in zip(identity,sorted_occseq):
        divconf = [ sorted_occseq_site[i] for i in row ]
        nthbutbestconfs.append(divconf)
    print "nth but best confs:", nthbutbestconfs

    divvaluesnth = [ sum(zip(*item)[0]) for item in zip(*nthbutbestconfs) ]

    for item in sorted( zip(divvaluesnth, zip(*nthbutbestconfs))):
        conf = zip(*item[1])[1]
        if not conf in confs:
            return conf
    else:
        print("# every nth also in conf. not possible to find a good structure")
        SystemExit('stop')













