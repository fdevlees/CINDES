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
import pickle # for saving and getting the tablebin
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

        #sets Gaussian09 input lines
        self.set_calculation_properties()
        return

    def __str__(self):
        sb=['Run object with the following attributes:']
        for key,value in sorted(self.__dict__.items()):
            if key in ['predictions']:
                sb.append("{key:20}=".format(key=key))
                sb.append( dump( value ) )
            elif key in ['TZmat','genalg', 'adj']:
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
            path = self.workdir  + '/databc'
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
            print "self.sites:", self.line1
            print "active: ", active
            print "sites: ", sites
            print "sites_adj:", sites_adj
        return sites_adj

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
                self.gaussianline2 = '#p geom=allcheck guess=read polar volume=tight '+param['functional']+'/'+param['basisset']+'\n'
            else:
                self.gaussianline = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '#p geom=allcheck guess=read polar '+param['functional']+'/'+param['basisset']+'\n'
        elif param['aip']==1 or param['aea']==1:
            self.gaussianline = '# opt=(maxcycle=100) ' + param['functional'] +'/'+ param['basisset'] +'\n'
            param['twojob']=1
            #param['multiplejobs']=1
            self.gaussianline2 = '# geom=check guess=read opt ' + param['functional'] +'/'+ param['basisset'] +'\n'
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
                if param['basisset'] in [ None, 0, '0', 'none', 'nalse', False, 'off' , 'n', 'na' ]:
                    self.gaussianline = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            if param['twojob'] == 1:
                self.gaussianline2 = '# geom=check guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            elif param['twojob'] == 2:
                self.multiplejobs = 2
                if param['semiempirical'] == 1:
                    self.gaussianline  = '# opt=(maxcycle=100) ' + 'pm6' +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '# geom=allcheck guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
            elif param['twojob'] == 4:
                self.multiplejobs = 4
                if param['semiempirical'] == 1:
                    self.gaussianline  = '# opt=(maxcycle=100) ' + 'pm6' +'\n'
                else:
                    self.gaussianline = '# opt=(maxcycle=100) scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                self.gaussianline2 = '# geom=allcheck guess=read scf=xqc ' + param['functional'] +'/'+ param['basisset'] +'\n'
                if param['solv']:
                    self.gaussianline_solv0 = '# geom=allcheck guess=read scf=xqc b3lyp/6-31G(d,p)\n'
                    self.gaussianline_solv1 = '# geom=allcheck guess=read scf=xqc scrf=(smd, solvent=aceticacid) b3lyp/6-31G(d,p)\n'

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
# 6. maximum configuration within run
# 7. global maximum
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
    try:
        tablename = myrun.tablename
    except AttributeError:
        tablename = 'tablebin'
    if myrun.restart>0:
        with open(tablename,'rb') as f:
            table = pickle.load(f)
        if True:
            # if all item[1] are ones:
            if all( item[1]==1 for item in table ):
                table = [[item[0]] + item[2:] for item in table ]
            else:
                table = [[item[0]] + item[1:] for item in table ]

    else:
        table = []
        open(tablename,'wb').close()
    return table

# 5 optimum at the start of the run
def set_maximum(myrun,table):
    if myrun.restart >= 3:
        tabledict = dict( [ item[0:2] for item in table ] )
        maximum = [ myrun.startconf, tabledict[myrun.startconf] ]
        logging.info("maximum:"+ str(maximum))
    else:
        maximum = 0
    return maximum

# 6 optimum within the global iterations
def testmax(myrun, mols, bcok):
    param = myrun.__dict__

    data = [ molecule.log() for molecule in mols ]
    for item in data:
        print item

    if 'bcprop' in param:
        if param['bcoptimum'] in ['min','Min','MIN']:
        #test if BC fullfilled. 
            try:
                voldoende = [ it for it in data if float(it[3]) < float(param['bcval']) ]
            except TypeError: pass
        else:
            assert param['bcoptimum'] in ['max','Max','MAX']
            try:
                voldoende = [ it for it in data if float(it[3]) > float(param['bcval']) ]
            except TypeError: pass

        print "voldoende:\n", pprint.pformat(voldoende, width=100)
        print "the BC condition is bc<:", param['bcval']
        if voldoende==[]: #so if there is at least one fullfilling BC
            bcok=0
            print "BC not fullfilled:"
            #maxsite = min(data,key = lambda x:x[2])
            maxsite = min(data,key = lambda x:abs( float(x[3]) - float(param['bcval']) ) )
        else: #BC nog niet
            bcok=1
            print "BC fullfilled; voldoende is not empty:", pprint.pformat(voldoende, width=100)
            #maxsite = max(voldoende,key = lambda x:x[1])
            if param['optimum']== 'minimum':
                maxsite = min(voldoende,key = lambda x:x[2])
            else:
                maxsite = max(voldoende,key = lambda x:x[2])
    else:
        bcok=1 #no BC but need this variable to test later on
        # maximum of the list or MINIMUM
        if param['optimum']== 'minimum':
            if param['cutoff'] == 0:
                maxsite = min(data,key = lambda x:x[2])
            else:
                testdata = [ item for item in data if abs(item[2]) > param['cutoff'] ]
                maxsite = min(testdata,key = lambda x:x[2])
        else:
            maxsite = max(data,key = lambda x:x[2])
    logging.warning('maxisite:' + pprint.pformat(maxsite, width=100))
    return maxsite,bcok

# 7 set global optimum and define convergence and redirect to Monte Carlo component
def runtest(run, maximum, maxsite, count, bcok,mctable=[], array=[]):
    param = run.__dict__
    TZmat= run.TZmat
    converged=0
    # test if this is same as previous maximum. if so then converged and break
    print
    if (count > 1 and bcok) or param['restart']>=3: #BCOK is a test of the boundary condition is already fullfilled
        if maximum[2] == maxsite[2]:  #test the property value! not 1 anymore!
            print "maximum is the same!"
            print "converged to a maximum configuration!"
            if param['montecarlo'] == 0:
               converged = 1
            else:
               if param['ml']==0:
                   maxsite = montecarloprocedure(param, array, maximum, mctable)
               else:
                   maxsite = montecarloprocedure(param, array, maximum, mctable, **TZmat)
               print "maxsite:", pprint.pfomat( maxsite, width=100 )
        else:
            print "maximum and maxsite are not the same yet"
            print "maximum:" ,maximum
            print "maxsite:" ,maxsite
    else: #except NameError:
        print "NameError no maximum or BC not yet fullfilled."
        #pass
    maximum = maxsite[:]
    return maximum, maxsite, converged

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
        for i in replaced:
            try:
                output += string.uppercase.index(i)
            except ValueError:
                output += string.lowercase.index(i)
        propx = output
        try:
            if 'bcprop' in param:
                propy= len(item.replace('_',''))
                molecule.boundaries = [ propy ]
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
    # set maximum
    maximum = set_maximum(myrun,table)
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
                print "maxsite[0]",maxsite[0]
                del startconf
                startconf = zcon.indtocon(maxsite[0])
                print "newconf: ", startconf


            # STEP 1: INDEXMAKER
            #get indices_all and the indices that still need to be calculated
            # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
            #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table )
            mols_todo, mols_nodo = zcon.classmaker2(startconf,array,k,table, myrun )
            if 1 in myrun.restrictions:
                mols_todo, mols_nodo = restriction1(mols_todo, mols_nodo, myrun )
            print "----- END random start configurations -----"
            print "indices_todo:",mols_todo
            print "data_nodo:", mols_nodo #all item[1]==1 in data_nodo 

            # STEP 2: PREDICTOR
            # perform prescreaning in a predictions. 
            mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, count, array=array, nsite=l)



            # STEP 3: SUBMITTING PART
            if not myrun.nosub==1:
                mols_all = submittingprocedure(mols_tocal,
                                               mols_nocal,
                                               myrun,
                                             **myrun.TZmat     ) # here call submitting procedure
            else: mols_all = skipper(mols_tocal,mols_nocal)
            print "mols_all:",mols_all

            # STEP 4: SORT
            # sort data in same order as allindices:
            # not necessary anymore in molsclass
            #mols_all = sorted(mols_all, key=lambda x:x.Pvalue)

            # STEP 5: UPDATE DATABASE and LOG results of microiteration
            # logs new elements in data to table and tablebin and whole data to cyclesinfo
            table = loggings(mols_all,table,count,k,l, made_pred )

            # STEP 6: UPDATE OPTIMUM STRUCTURE
            # decide what the maximum site is and if the bc if fullfilled
            print "BCOK:", bcok
            maxsite, bcok = testmax(myrun, mols_all, bcok)

            print("--- %s seconds ---" % (time.time() - myrun.starttime))
            print(myrun.currenttime())
        # END LOOP OVER SITES

        #get maximum and test convergence
        maximum, maxsite,converged = runtest(myrun, maximum, maxsite, count, bcok, mctable=table, array = array)
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
    '''generate all structures and print in format'''
    from CINDES4.utils.converter import Converter
    from CINDES4.utils import writings
    from CINDES4.predictor import learning_skl as learning
    converter = Converter()
    #get structure
    TZmat = r.geometry(**param)
    #to get an xyz file with the data from tablebin do generate1()
    with open('tablebin','rb') as f:
        table = pickle.load(f)
    #print table[508:510]
    learning.generate1(converter=converter,table=table,**TZmat)

    #to get an xyz file with all the possible structures possible:
    #generate2(core,active,passive,converter, **param)
    return

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
    myrun = Run(**param)
    print myrun
    table = set_table(myrun)
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
    # set maximum
    maximum = set_maximum(myrun,table)
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
            print "maxsite[0]",maxsite[0]
            del startconf
            startconf = zcon.indtocon(maxsite[0])
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
        # decide what the maximum site is and if the bc if fullfilled
        print "BCOK:", bcok
        maxsite, bcok = testmax(myrun, mols_all, bcok)

        if myrun.procedure=='steepest2':
            maxconf = zcon.indtocon(maxsite[0])
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

        #get maximum and test convergence
        maximum, maxsite,converged = runtest(myrun, maximum, maxsite, count, bcok, mctable=table, array = array)
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
    return

