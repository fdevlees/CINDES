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

def BFS(run):
    """ 1. This is the standard BFS procedure """
    logging.info(run)  # this should print all the class elements via the __str__ function
    mybfs = BestFirstSearch(run)
    result = mybfs.evolve()
    return result

def SD(run):
    """ 1. This is the SD procedure """
    logging.info(run)  # this should print all the class elements via the __str__ function
    mySD = SteepestDescent(run)
    result = mySD.evolve()
    return result

def runEmpty(run):
    logging.info(run)
    myEmpty = Empty(run)
    result = myEmpty.evolve()
    return result

class Algorithm(object):
    def __init__(self, run):

        self.run = run
        self.array = self.run.array

        self.bcok = 0  # TO REMOVE LATER
        self.count = 0
        self.ncalcs = 0

        # the table with all the results of all calculated configs
        self.table = set_table(self.run, self.run.array)

        # set initial optimum
        self.optimum = None
        self.history = []
        return

    def evaluateMols(self, mols, k=0, l=0):
        ''' k:number of site
            l:l^th site, if l==0 new ML model will be trained
        '''

        # STEP 2&3: predict and calculate
        mols_all, nnewcalcs, made_pred = calculator.evaluate_mols(self.run, mols, self.table, self.count, nsite=l)
        self.ncalcs += nnewcalcs

        # STEP 4: UPDATE OPTIMUM STRUCTURE
        # decide what the optimum site is and if the bc if fullfilled
        if not self.run.procedure=='empty':
            optsite = self.testmax(mols_all)
        else:
            optsite = mols[0]

        # STEP 5: UPDATE DATABASE and LOG results of microiteration
        # logs new elements in data to table and tablebin and whole data to cyclesinfo
        self.table = loggings(mols_all,
                         self.table,
                         self.count,
                         k, l,
                         made_pred,
                         tablename=self.run.tablename,
                         write = self.run.write)

        logging.info("BCOK:{:d}".format(self.bcok))
        self.history.append({
            'count':self.count,
            'p':optsite.Pvalue,
            'index':optsite.index,
            'site':k,
            'ncalcs':self.ncalcs })
        logging.info("--- %s seconds ---" % (time.time() - self.run.starttime))
        logging.debug(self.run.currenttime())
        return optsite

    def testmax(self, mols):
        ''' sets optsite
            NOTE: multiple boundary conditions are not yet implemented
        '''
        param = self.run.__dict__
        # select only molecules that are really calculated:
        mols = [mol for mol in mols if mol.predicted == False and not mol.Pvalue is None]

        if 'bcprop' in param:
            if param['bcoptimum'] == 'min':
                # test if BC fullfilled.
                satisfactory = [mol for mol in mols if mol.boundaries[0] < float(param['bcval'])]
            else:
                assert param['bcoptimum'] == 'max'
                satisfactory = [mol for mol in mols if mol.boundaries[0] > float(param['bcval'])]

            print "satisfactory:\n", pprint.pformat(satisfactory, width=100)
            print "the BC condition is bc<:", param['bcval']

            if satisfactory == []:  # so if there is at least one fullfilling BC
                self.bcok = 0
                print "BC not fullfilled:"
                optsite = min(mols, key=lambda mol: abs(mol.boundaries[0] - float(param['bcval'])))
            else:
                self.bcok = 1
                print "BC fullfilled; satisfactory is not empty:", pprint.pformat(satisfactory, width=100)
                if param['optimum'] == 'minimum':
                    optsite = min(satisfactory, key=lambda x: x[2])
                else:
                    optsite = max(satisfactory, key=lambda x: x[2])
        else:
            self.bcok = 1  # no BC but need this variable to test later on
            # optimum of the list or MINIMUM
            if param['optimum'] == 'minimum':
                if param['cutoff'] == 0:
                    optsite = min(mols, key=lambda x: x.Pvalue)
                else:
                    testmols = [mol for mol in mols if abs(mol.Pvalue) > param['cutoff']]
                    optsite = min(testmols, key=lambda mol: mol.Pvalue)
            else:
                optsite = max(mols, key=lambda mol: mol.Pvalue)
        optsite.opt = True
        return optsite

class Empty(Algorithm):
    def evolve(self):
        from CINDES.molecule import Molecule
        mols = [Molecule([])]
        optsite = self.evaluateMols(mols)
       

class BestFirstSearch(Algorithm):
    def __init__(self, run):
        super(BestFirstSearch, self).__init__(run)
        return

    def evolve(self):
        # ------------------------------------- #
        # --- HERE THE MAIN LOOP STARTS --- --- #
        # ------------------------------------- #
        self.count = 1  # so we start counting at 1!
        startconf = None

        while True:
            print_title("Global Iteration No.: " + str(self.count), outline='l', signator="-")

            optsite, startconf = self.runCycle(startconf)

            # get optimum and test convergence
            optsite, converged = self.runtest(optsite)

            if converged == 1:
                break
            self.count += 1
            if self.count > self.run.maxiter:
                logging.warning("maxiterations is reached")
                logging.warning("optimum is:"+ pprint.pformat(self.optimum))
                break
        # ---------------------------- #
        # ------ END OF LOOPING ------ #
        # ---------------------------- #
        results = {
                'index':self.optimum.index,
                'gen':self.count,
                'history':self.history,
                'ncalcs':self.ncalcs,
                'p':self.optimum.Pvalue}
        logging.warning("ALGORITHM DONE")
        return results

    def runCycle(self, startconf, *args):
        # set site order in sequence
        sequence = self.get_sequence()

        # for each site in sequence:
        for l in range(len(sequence)):
            k = sequence[l]
            print_title("k(site)= {} l(nsite)= {} (c={})".format(k, l, self.count), outline='l', signator='=')
            if not startconf:
                startconf = self.get_startconf()

            optsite = self.optimizeSite(k, l, startconf)
            logging.info("optsite:" + pprint.pformat(optsite))

            startconf = zcon.indtocon(optsite.index)

        return optsite, startconf

    def optimizeSite(self, k, l, startconf):
        '''optimizes the k^th site'''
        mols = self.getMols(startconf, k)
        optsite = self.evaluateMols(mols, k, l)
        return optsite

    def getMols(self, startconf, k):
        mols = zcon.get_molecules(startconf, self.array, k, self.run)
        return mols

    def get_startconf(self):
        logging.info("random start molecule: ")
        startconf = []
        startind = self.run.startind
        if not startind == '':
            startconf = zcon.indtocon(startind)
            logging.info("startconf read from input")
        else:
            for i in range(len(self.run.array)):
                startconf.append(random.choice(self.run.array[i]))
            logging.info("constructed random start configuration")
        logging.info("startconf:" + pprint.pformat(zcon.contoind(startconf)))
        return startconf

    def get_sequence(self):
        # START set sequence INPUT: param
        param = self.run.__dict__
        nsites = param['nsites']
        if 'sequences' in param:
            try:
                sequence = param['sequences'][self.count - 1]  # accounting for the fact count starts counting at 1
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

    # 7 set global optimum and define convergence and redirect to Monte Carlo component
    def runtest(self, optsite):
        converged = 0

        if (self.count > 1 and self.bcok):  # bcok is a test for the boundary condition
            if self.optimum == optsite:
                logging.warning("optimum is the same as in previous iteration: converged to a optimum configuration!")
                if self.run.montecarlo == 0:
                    converged = 1
                else:
                    from montecarlo import montecarloprocedure
                    property_table = get_property_table(self.table, self.run)
                    if self.run.ml == 0:
                        optsite = montecarloprocedure(self.run, self.array, self.optimum, property_table)
                    else:
                        optsite = montecarloprocedure(self.run, self.array, self.optimum, property_table, **self.run.TZmat)
                    logging.info("optimum of current iteration:" + pprint.pformat(optsite, width=100))
            else:
                logging.info("optimum last global iteration:" + pprint.pformat(self.optimum))
                logging.info("optimum  of current iteration:" + pprint.pformat(optsite))
        elif not self.bcok:
            logging.info("boundary conditions are not yet fullfilled.")
        self.optimum = optsite.copy()
        return optsite, converged


class SteepestDescent(BestFirstSearch):
    def __init__(self, run):
        super(SteepestDescent, self).__init__(run)

        # defines which sites will be changed. only relevant for reduced steepest algorithm
        self.restingsites = range(self.run.nsites)
        return

    def runtest(self, optsite):
        # for RSD algorithm:
        if self.restingsites == []:
            return optsite, True
        else:
            return super(SteepestDescent, self).runtest(optsite)

    def runCycle(self, startconf):
        # define new starting geometry
        print "startconf:", startconf
        if not startconf:
            startconf = self.get_startconf()
        print "start configuration: ", startconf

        optsite = self.optimize(startconf)

        if self.run.procedure == 'rsd':
            maxconf = zcon.indtocon(optsite.index)
            print "maxconf:", maxconf, 'while startconf:', startconf
            try:
                changedsite = [siteM == siteS for siteM, siteS in zip(maxconf, startconf)].index(False)
            except ValueError:
                print "no site changed"
            else:
                self.restingsites.remove(changedsite)
            if self.restingsites == []:
                print "all sites changed once"
            else:
                print "changed site:", changedsite, "sites still to be optimized:", self.restingsites

        startconf = optsite
        return optsite, startconf

    def optimize(self, startconf):
        mols = self.getMols(startconf)
        optsite = self.evaluateMols(mols, startconf)
        return optsite

    def getMols(self, startconf):
        mols = zcon.get_molecules_SD(startconf, self.array, self.run, self.restingsites)
        return mols
