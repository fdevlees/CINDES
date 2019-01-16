#!/bin/env python
# this line must be at the beginning of the file!
from __future__ import division

# debug flag
debug = 1

# import python libraries
import logging  # instead of the large amount of print statements not using it at the moment
import numpy as np

# import my own modules
import construction as zcon  # all functions needed for constructing new geometries
from predictions import predictor
from loggings import loggings
import calculator
from run import FrameRun

# import utils
from CINDES.utils.writings import print_title
from CINDES.utils.table import set_table, get_property_table
from CINDES.utils.molecule import Molecule


class Algorithm(object):
    def logpopulation(self, mols_todo, mols_nodo):
        logpop = []
        p=logpop.append
        p("|      NEW POPULATION CONSTRUCTED:")
        p("|   mols_todo:")
        if mols_todo:
            for mol in mols_todo:
                p("|      {}".format(mol))
        else:
            p("|      -")
        p("|   mols_nodo:")
        if mols_nodo:
            for mol in mols_nodo:
                p("|      {}".format(mol))
        else:
            p("|      -")
        logging.info('\n'.join(logpop))
        return


def GP(param, array=None):
    """ This will be a Bayes Optimization: Gaussian Process """
    if array is None:
        self.array = param['array']
    myrun = FrameRun(**param)
    logging.info(myrun)
    my_gp = GaussianProcess(myrun)
    result = my_gp.evolve()
    return result


class GaussianProcess(Algorithm):
    def __init__(self, run):
        self.run = run
        self.array = self.run.array
        # make an array with CNOO instead of ['C','N','O','O'] etc
        self.array_joined = [ tuple("".join(item) for item in site ) for site in self.array ]
        self.bcok = 0  # TO REMOVE LATER
        self.ncalcs = 0
        self.skip = self.run.nosub==1

        # the table with all the results of all calculated configs
        self.table = set_table(self.run, self.run.array)
        self.property_table = get_property_table(self.table, self.run)

        if self.run.optimum=='maximum':
            raise NotImplementedError('only minimization implemented')

        # set initial optimum
        self.optimum = None
        self.history = []
        return

    def set_space(self):
        print "array:", self.array
        from skopt.utils import normalize_dimensions
        space = normalize_dimensions(self.array_joined)
        print "space:", space
        self.space = space
        return

    def get_Y(self, X, gen):

        # 1. get set of unique confs
        indices = tuple( '_'.join(item) for item in X )
        unique_confs = [ zcon.indtocon(item) for item in set(indices)]

        # 2. make molecules
        individuals = [Molecule(conf=conf) for conf in unique_confs]

        # 3. check if already in table
        mols_todo, mols_nodo = zcon.check_in_table(individuals, self.table, self.run.props)
        self.logpopulation(mols_todo, mols_nodo)

        # irrelevant for now, just renaming:
        if self.run.predictions:
            property_table = get_property_table(self.table, self.run)
            mols_nocal, mols_tocal, made_pred = predictor(
                self.run,
                property_table,
                mols_todo, mols_nodo,
                gen,
                array=self.array,
                nsite=0
            )
        else:
            mols_nocal, mols_tocal = mols_nodo, mols_todo

        # 4. calculate configurations
        self.ncalcs += len(mols_tocal)
        if self.skip:
            mols_all = skipper(mols_tocal, mols_nocal, self.run)
        else:
            mols_all = calculator.procedure(self.run, mols_tocal, mols_nocal)
            #mols_all = submittingprocedure(mols_tocal,
            #                               mols_nocal,
            #                               self.run)
        self.optimum, _ = testmax(self.run, mols_all)
        print "mols_all:", mols_all

        # 5. log new results
        self.table = loggings(mols_all,
                              self.table,
                              gen,
                              1, 1,
                              made_pred=False,
                              tablename=self.run.tablename,
                              write=self.run.write)
        self.history.append({
                    'count':gen,
                    'p':self.optimum.Pvalue,
                    'index':self.optimum.index,
                    'ncalcs':self.ncalcs
                    })

        # 6. set Y
        y_dict = {mol.index: mol.Pvalue for mol in mols_all}
        #print "y_dict:", y_dict, "indices:", indices
        Y=[]
        for index in indices:
            Y.append(y_dict[index])
        return Y

    def evolve(self):
        batchsize=10
        self.set_space()

        from skopt import Optimizer
        optimizer = Optimizer(
            dimensions=self.space,
            base_estimator="GP",
            n_initial_points=batchsize,
            random_state=1,
            acq_optimizer='sampling'
            )


        print "optimizer model:", optimizer.base_estimator_
        print "eta: {}, acq-function: {}, acq-optimizer: {}".format(optimizer.eta, optimizer.acq_func, optimizer.acq_optimizer)


        for gen in range(self.run.maxiter):
            print_title("BATCH-NO: " + str(gen), outline='l', signator="-")

            X = optimizer.ask(n_points=batchsize)
            Y = self.get_Y(X, gen=gen)
            X, Y = self.only_finite(X, Y)
            optimizer.tell(X, Y)

            # save model:
            with open('my-optimizer.pkl', 'wb') as f:
                pickle.dump(optimizer, f)

            logging.info("--- %s seconds ---" % (time.time() - self.run.starttime))

        print "\n\tOptimum:", min(zip(optimizer.yi, optimizer.Xi))

        return

    @staticmethod
    def only_finite(X, Y):
        print "X, Y before:", X, Y
        X, Y = zip(*[ (x,y) for x,y in zip(X,Y) if np.isfinite(y) ])
        print "X, Y after:", X, Y
        return X, Y

