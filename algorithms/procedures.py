#!/bin/env python
''' This module contains some other procedures that can be useful exept for searching:
    1. genconf: to make inputfiles for a single molecular index given by startind
    2. generate all molecules given in input
    3. generate all possible molecules
    4. get n random structures
    5. test predictions from a dataset
    6. do a Steepest Descent algorithm.
'''
# this line must be at the beginning of the file!
from __future__ import division

# debug flag
debug = 1

# import python libraries
import time  # for getting time/date and time delays
start = time.clock()
import pprint  # pretty printer for printing lists
import os  # for getting window width and testing existence of files
import random  # for obtaining random geometry
import logging  # instead of the large amount of print statements not using it at the moment
from copy import deepcopy  # for keeping matrices while changing others

# import my own modules
from CINDES.evaluation import construction as zcon  # all functions needed for constructing new geometries
from CINDES.evaluation.predictions import predictor
from CINDES.algorithms.loggings import loggings
from CINDES.evaluation import calculator

# import utils
from CINDES.utils.writings import print_title
from CINDES.utils.utils import skipper
from CINDES.utils.table import set_table, get_property_table
from BFS import testmax, runtest

once = 0

print "time for imports:", time.clock() - start


def genconf(myrun):
    """ 2. Procedure to generate the inputfiles for a single index """
    from CINDES.utils.molecule import Molecule
    param = myrun.__dict__
    conf = zcon.indtocon(param['startind'])
    mol = Molecule(conf=conf)
    print "defaultgroups:", myrun.defaultgroups
    print "in GENCONF: conf is:", conf
    param['workdir'] = os.getcwd()
    path = param['workdir']
    param["path"] = str(path)

    if isinstance(myrun.zmatrixfile, list):
        print "myrun.zmatrixfile:", myrun.zmatrixfile
        for zmatfile in myrun.zmatrixfile:
            tzmat = myrun.TZmatrices[zmatfile]
            c, a, p = map(deepcopy, (tzmat['core'], tzmat['active'], tzmat['passive']))
            mtzmat = zcon.constructor2(mol.conf, c, a, p, links=myrun.symlinks, defaultgroups=myrun.defaultgroups)
            setattr(mol, zmatfile, mtzmat)
    else:
        print "zmatrixfile = single file"
        TZmat = myrun.TZmat
        c = deepcopy(TZmat['core'])
        a = deepcopy(TZmat['active'])
        p = deepcopy(TZmat['passive'])
        zmat = zcon.constructor2(conf, c, a, p, links=myrun.symlinks, defaultgroups=myrun.defaultgroups)
        setattr(mol, myrun.zmatrixfile, zmat)
        setattr(mol, 'zmat', zmat)
    if myrun.program == 'gaussian':
        from CINDES.evaluation import gaussian as program
    elif myrun.program == 'nwchem':
        from CINDES.evaluation import nwchem as program
    else:
        raise SystemExit('program not recognized')
    calcs = myrun.calcs
    for calc in calcs:
        if isinstance(calc, list):
            calc = calc
            for cal in calc:
                program.filewriter(mol, cal)
                print "\n\tprinter file with:"
                pprint.pprint(cal)
        else:
            calc = calcs[0]
            program.filewriter(mol, calc)
            print "\n\tprinter file with:"
            pprint.pprint(calc)
    return


def generate_procedure(myrun):
    ''' 3. calculate all possible structures '''

    table = set_table(myrun)
    array = myrun.array
    print myrun

    # get all structures
    print "len table:", len(table)

    if hasattr(myrun, 'ngenerate'):
        from CINDES.utils.molecule import Molecule
        mols = [ Molecule(conf=zcon.indtocon(ind)) for ind in myrun.generatemols ]
        print "mols:", mols
    else:
        mols = get_all_molecules(array)
    # 1b check already in database
    mols_todo, mols_nodo = zcon.check_in_table(mols, table, myrun.props)
    # 1c eventueel predictions
    mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo, mols_nodo, 0, array=array)

    # calculate them
    #print_title("COUNT: " + str(count),outline='l',signator="-")

    count = 0
    if myrun.batchsize is None:
        print "molecules:"
        for i, mol in enumerate(mols):
            print i, mol.index
        if not myrun.nosub == 1:
            mols_all = calculator.procedure(myrun, mols_tocal, mols_nocal)
            #mols_all = submittingprocedure(mols_tocal,
            #                               mols_nocal,
            #                               myrun,
            #                               )  # here call submitting procedure
        else:
            mols_all = skipper(mols_tocal, mols_nocal, myrun)
    elif myrun.batchsize:
        # use batches

        def chunks(l, n):
            '''yields successive n-sized chunks of l'''
            for i in range(0, len(l), n):
                print "    yielding:", i, "to:", i+n, "from total:", len(l)
                yield l[i:i + n]

        mols_all = []
        for i, batch in enumerate(chunks(mols_tocal, myrun.batchsize)):
            print "chunk nr:", i, "with ", len(batch), "structures"
            if not myrun.nosub == 1:
                batch = calculator.procedure(myrun, batch, [])
                #batch = submittingprocedure(batch,
                #                            [],
                #                            myrun,
                #                            )  # here call submitting procedure
            else:
                batch = skipper(batch, [], myrun)

            # intermediate logging:
            table = loggings(batch, table, i, 1, 1, tablename=myrun.tablename)
            mols_all.extend(batch)
    else:
        # use job arrays.
        # make a jobscript with the line:
        # qsub -t 1-njobs
        # $PBS_ARRAY_INDEX has to be used in the submitscript to submit each job.
        # each number has to refer to a certain filename.
        # make file with all the filenames:
        with open('./CALC/filenames.txt', 'w') as fout:
            fout.write('\n'.join([myrun.identify + mol.index for mol in mols_tocal]))
        print "njobs:", len(mols_tocal)

        # generate all inputfiles
        calculator.geommaker(mols_tocal, myrun)
        calculator.filemaker(mols_tocal, myrun)  # ----------------------------------HERE IS THE FILEWRITER CALL

    # final logging
    print "mols_all:", mols_all
    table = loggings(mols_all, table, count, 1, 1, made_pred=made_pred, tablename=myrun.tablename)
    print "DONE"
    return


def get_all_molecules(array):
    from CINDES.utils.molecule import Molecule
    print "in get_all_molecules"
    #A = [ map(''.join,item) for item in array ]
    A = array
    C = [[]]
    for site in A:
        D = []
        for item in C:
            for group in site:
                conf = item + [group]
                D.append(conf)
        C = D
    mols = [Molecule(conf=conf, dihedral=True) for conf in C]

    return mols


def generate2(core, active, passive, converter, **kwargs):
    from itertools import product
    confs = []
    for item in product(*array):
        confs.append(item)
    printindices = 1
    # confs to xyz:
    with open('chemspace.xyz', 'a') as fid:

        for i in range(len(confs)):
            c = deepcopy(core)
            a = deepcopy(active)
            p = deepcopy(passive)
            logging.debug("i=" + str(i))
            mat = zcon.constructor2(confs[i], c, a, p, **kwargs)
            xyz = zmatoxyz(converter, mat)
            if printindices == 1:
                fid.write('{:05d} {:4d} {:s}\n'.format(i + 1, len(xyz), zcon.contoind(confs[i])))
            else:
                fid.write('{:05d} {:4d}\n'.format(i + 1, len(xyz)))
            # fid.write('{:12.8f}\n'.format(Y[i])) #no property for chemspace
            for item in xyz:
                fid.write('{:3s} {:10.4f} {:10.4f} {:10.4f}\n'.format(item[0], item[1][0], item[1][1], item[1][2]))
            fid.write('\n')
            # now for each mat transform to xyz
    return confs

# 4: genrandom


def genrandom(param, array):
    ''' generate random structures and print '''
    print
    for _ in xrange(param['nrandom']):
        conf = []
        for i in range(len(array)):
            group = random.choice(array[i])
            conf.append(group)
        print zcon.contoind(conf)
    print
    return True

# 5: testrun


def testrun(param, array):
    pass

# 5b: testpred


def testpred(param, array):
    ''' run the predictions on the tablebin file '''
    print_title("Testing Prediction procedure activated!", outline='l', signator=':')

    class Mol(object):
        def __init__(self):
            self.predictions = {}

        def __repr__(self): return "<empty molecule object>"
    myrun = FrameRun(**param)
    print myrun
    table = set_table(myrun)
    property_table = get_property_table(table, myrun)
    mols_todo, mols_nodo = ([Mol(), ], [Mol(), ])
    mols_nocal, mols_tocal, made_pred = predictor(myrun, property_table, mols_todo, mols_nodo, 99, array=array, nsite=0)
    return

# 6: steepest descent

def SteepestDescent(param, array):
    raise NotImplementedError('this method is too old and out of date compared to the BFS procedure')
    bcok = 0  # TO REMOVE LATER
    param['bcok'] = 0

    # SET MYRUN CLASS and assign all necessary attributes
    myrun = FrameRun(**param)
    print(myrun)  # this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun)
    # set optimum
    optimum = None
    # set calculation properties
    startconf = get_startconf(param, array)
    # END MYRUN CLASS assignments. from now myrun should contain all the
    # necessary information to work with during the whole program run.

    # defines which sites will be changed. only relevant for steepest2 algorithm
    myrun.restingsites = range(myrun.nsites)

    # ------------------------------------- #
    # --- HERE THE MAIN LOOP STARTS --- --- #
    # ------------------------------------- #
    count = 1  # so we start counting at 1!
    ncalcs = 0
    while True:
        print_title("COUNT: " + str(count), outline='l', signator="-")

        if count > 1:  # define new startconfiguration if not first cycle
            # define new starting geometry
            print "optsite[0]", optsite[0]
            del startconf
            startconf = optsite.conf
            #startconf = zcon.indtocon(optsite[0])
            print "newconf: ", startconf

        # STEP 1: INDEXMAKER
        # get indices_all and the indices that still need to be calculated
        #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker_SD(startconf, array, table, myrun )
        mols = get_molecules_SD(startconf, array, myrun)

        # STEP 2&3: predict and calculate
        mols_all, nnewcalcs, made_pred = calculator.evaluate_mols(myrun, mols, table, property_table, count, nsite=0)
        ncalcs += nnewcalcs

        # STEP 4: UPDATE DATABASE and LOG results of microiteration
        # logs new elements in data to table and tablebin and whole data to cyclesinfo
        table = loggings(mols_all, table, count, 1, 1, made_pred=made_pred, tablename=myrun.tablename)

        # STEP 5: UPDATE OPTIMUM STRUCTURE
        # decide what the optimum site is and if the bc if fullfilled
        print "BCOK:", bcok
        optsite, bcok = testmax(myrun, mols_all, bcok)

        if myrun.procedure == 'steepest2':
            maxconf = zcon.indtocon(optsite[0])
            print "maxconf:", maxconf, 'while startconf:', startconf
            try:
                changedsite = [siteM == siteS for siteM, siteS in zip(maxconf, startconf)].index(False)
            except ValueError:
                print "no site changed"
            else:
                myrun.restingsites.remove(changedsite)
            if myrun.restingsites == []:
                print "all sites changed once"
                break

        print("--- %s seconds ---" % (time.time() - myrun.starttime))
        print(myrun.currenttime())
        # END LOOP OVER SITES

        # get optimum and test convergence
        optimum, optsite, converged = runtest(myrun, optimum, optsite, count, bcok, mctable=table, array=array)
        if converged == 1:
            break
        count += 1
        if count > param['maxiter']:
            print "maxiterations is reached"
            print "optimum is: ", optimum
            break
    # ---------------------------- #
    # ------ END OF LOOPING ------ #
    # ---------------------------- #
    print "DONE"
    return
