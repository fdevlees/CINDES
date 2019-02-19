'''This is the module containing the Best First Search algorithm'''

import logging
import random
import pprint
import time

from CINDES.utils.writings import print_title

from CINDES.evaluation import construction as zcon
from CINDES.evaluation import calculator
from CINDES.utils.table import set_table, get_property_table
from loggings import loggings

def get_startconf(param, array):
    logging.info("random start molecule: ")
    startconf = []
    if not param['startind'] == '':
        startconf = zcon.indtocon(param['startind'])
        logging.info("read startconf from input")
    else:
        for i in range(len(array)):
            startconf.append(random.choice(array[i]))
        logging.info("constructed random start configuration")
    logging.info("startconf:" + pprint.pformat(zcon.contoind(startconf)))
    return startconf


def testmax(myrun, mols, bcok=None):
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
    mols = [mol for mol in mols if mol.predicted == False and not mol.Pvalue is None]

    if 'bcprop' in param:
        if param['bcoptimum'] in ['min', 'Min', 'MIN']:
            # test if BC fullfilled.
            try:
                satisfactory = [mol for mol in mols if mol.boundaries[0] < float(param['bcval'])]
            except TypeError as e:
                print e
                pass
        else:
            assert param['bcoptimum'] in ['max', 'Max', 'MAX']
            try:
                satisfactory = [mol for mol in mols if mol.boundaries[0] > float(param['bcval'])]
            except TypeError as e:
                print e
                pass

        print "satisfactory:\n", pprint.pformat(satisfactory, width=100)
        print "the BC condition is bc<:", param['bcval']
        if satisfactory == []:  # so if there is at least one fullfilling BC
            bcok = 0
            print "BC not fullfilled:"
            #optsite = min(data,key = lambda x:x[2])
            optsite = min(mols, key=lambda mol: abs(mol.boundaries[0] - float(param['bcval'])))
        else:  # BC not yet
            bcok = 1
            print "BC fullfilled; satisfactory is not empty:", pprint.pformat(satisfactory, width=100)
            #optsite = max(satisfactory,key = lambda x:x[1])
            if param['optimum'] in ['minimum', 'min']:
                optsite = min(satisfactory, key=lambda x: x[2])
            else:
                optsite = max(satisfactory, key=lambda x: x[2])
    else:
        bcok = 1  # no BC but need this variable to test later on
        # optimum of the list or MINIMUM
        if param['optimum'] in ['minimum', 'min']:
            if param['cutoff'] == 0:
                optsite = min(mols, key=lambda x: x.Pvalue)
            else:
                testmols = [mol for mol in mols if abs(mol.Pvalue) > param['cutoff']]
                optsite = min(testmols, key=lambda mol: mol.Pvalue)
        else:
            optsite = max(mols, key=lambda mol: mol.Pvalue)

    optsite.opt = True
    #logging.warning('optsite:' + pprint.pformat(optsite, width=100))
    return optsite, bcok

# 7 set global optimum and define convergence and redirect to Monte Carlo component


def runtest(run, optimum, optsite, count, bcok, table=[], array=[]):
    param = run.__dict__
    converged = 0

    #raise SystemExit('optimum and optsite should be Molecule instances now')
    if (count > 1 and bcok):  # BCOK is a test of the boundary condition is already fullfilled
        if optimum == optsite:  # test the property value! not 1 anymore!
            logging.warning("optimum is the same! converged to a optimum configuration!")
            if param['montecarlo'] == 0:
                converged = 1
            else:
                from montecarlo import montecarloprocedure
                property_table = get_property_table(table, run)
                if param['ml'] == 0:
                    optsite = montecarloprocedure(run, array, optimum, property_table)
                else:
                    optsite = montecarloprocedure(run, array, optimum, property_table, **run.TZmat)
                logging.info("optimal_after_this_site:" + pprint.pformat(optsite, width=100))
        else:
            logging.info("global_iteration_optimum and optimum_after_this_site are not the same yet")
            logging.info("gi_optimum:" + pprint.pformat(optimum))
            logging.info("current optimum:" + pprint.pformat(optsite))
    else:  # except NameError:
        logging.info("NameError no optimal structure or BC not yet fullfilled.")
        # pass
    optimum = optsite.copy()
    return optimum, optsite, converged



# THERE ARE DIFFERENT GLOBAL PROGRAM FLOW PROCEDURES:
# 1: STANDARD PROCEDURE: Best First Search: BFS()
# 2: Generate 1 Configuration input file: genconf
# 3: Generate total chemical space defined by the sites and functionalisations: generate
# 4: Generate a number of random structures and print them to screen: genrandom
# 5: A testrun. Not implemented. a helper function for the test functions in ./tests/tests.py: testrun
# 6: Steepest Descent algorithm. Looks like BFS but there is no loop over sites
# 7: Generate database based on farthest point selection. (based on diversity index)

def BFS(run):
    """ 1. This is the standard BFS procedure """

    logging.info(run)  # this should print all the class elements via the __str__ function
    mybfs = BestFirstSearch(run)
    result = mybfs.evolve()
    return result




class BestFirstSearch(object):
    def __init__(self, run):
        self.run = run
        self.array = self.run.array
        self.bcok = 0  # TO REMOVE LATER
        self.startconf = get_startconf(run.__dict__, self.run.array)

        # the table with all the results of all calculated configs
        self.table = set_table(self.run, self.run.array)

        # set initial optimum
        self.optimum = None
        self.history = []
        return

    def evolve(self):
        # ------------------------------------- #
        # --- HERE THE MAIN LOOP STARTS --- --- #
        # ------------------------------------- #
        count = 1  # so we start counting at 1!
        ncalcs = 0
        while True:
            print_title("Global Iteration No.: " + str(count), outline='l', signator="-")

            # set site order in sequence INPUT: param, count
            sequence = self.get_sequence(count)

            # for each site in sequence:
            for l in range(len(sequence)):
                k = sequence[l]
                print_title("k(site)= {} l(nsite)= {} (c={})".format(k, l, count), outline='l', signator='=')
                if not l == 0 or count > 1:  # define new startconfiguration if not first cycle
                    # define new starting geometry
                    logging.info("optsite:" + pprint.pformat(optsite))
                    del self.startconf
                    self.startconf = zcon.indtocon(optsite.index)

                # STEP 1: INDEXMAKER
                # get indices_all and the indices that still need to be calculated
                # if table is correctly formatted all second element item[1]==1. meaning they are ab-initio calculated
                #indices_todo,data_nodo,configurations,indices_all = zcon.indexmaker2(startconf,array,k,table )
                mols = zcon.get_molecules(self.startconf, self.array, k, self.run)

                # STEP 2&3: predict and calculate
                mols_all, nnewcalcs, made_pred = calculator.evaluate_mols(self.run, mols, self.table, count, nsite=l)
                ncalcs += nnewcalcs


                # STEP 4: UPDATE OPTIMUM STRUCTURE
                # decide what the optimum site is and if the bc if fullfilled
                optsite, self.bcok = testmax(self.run, mols_all, self.bcok)

                # STEP 5: UPDATE DATABASE and LOG results of microiteration
                # logs new elements in data to table and tablebin and whole data to cyclesinfo
                self.table = loggings(mols_all,
                                 self.table,
                                 count,
                                 k, l,
                                 made_pred,
                                 tablename=self.run.tablename)

                logging.info("BCOK:{:d}".format(self.bcok))
                self.history.append({
                    'count':count,
                    'p':optsite.Pvalue,
                    'index':optsite.index,
                    'site':k,
                    'ncalcs':ncalcs,
                    'startind':zcon.contoind(self.startconf)})
                logging.info("--- %s seconds ---" % (time.time() - self.run.starttime))
                logging.debug(self.run.currenttime())
            # HERE ENDS LOOP OVER SITES

            # get optimum and test convergence
            self.optimum, optsite, converged = runtest(self.run,
                    self.optimum,
                    optsite,
                    count,
                    self.bcok,
                    table=self.table,
                    array=self.array)

            if converged == 1:
                break
            count += 1
            if count > self.run.maxiter:
                logging.warning("maxiterations is reached")
                logging.warning("optimum is:"+ pprint.pformat(self.optimum))
                break
        # ---------------------------- #
        # ------ END OF LOOPING ------ #
        # ---------------------------- #
        results = {
                'index':self.optimum.index,
                'gen':count,
                'history':self.history,
                'ncalcs':ncalcs,
                'p':self.optimum.Pvalue}
        logging.warning("BFS DONE")
        return results

    def get_sequence(self, count):
        # START set sequence INPUT: param
        param = self.run.__dict__
        nsites = param['nsites']
        if 'sequences' in param:
            try:
                sequence = param['sequences'][count - 1]  # accounting for the fact count starts counting at 1
            except IndexError:
                sequence = random.sample(range(nsites), nsites)
            finally:
                logging.warning("SEQUENCE: " + str(sequence))
                return sequence
        if param['norandom'] == 1:
            sequence = range(nsites)
        else:
            sequence = random.sample(range(nsites), nsites)
        # output sequence
        logging.info("SEQUENCE: " + str(sequence))
        return sequence


