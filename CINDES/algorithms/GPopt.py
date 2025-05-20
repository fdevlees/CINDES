#!/bin/env python
# this line must be at the beginning of the file!


# debug flag
debug = 1

# import python libraries
import logging  # instead of the large amount of print statements not using it at the moment
import numpy as np
import pickle
import time

# import my own modules
from CINDES.evaluation import construction as zcon  # all functions needed for constructing new geometries
from CINDES.evaluation.predictions import predictor
from .loggings import loggings
from .BFS import testmax
from CINDES.evaluation.calculator import evaluate_mols

# import utils
from CINDES.utils.writings import print_title
from CINDES.utils.table import set_table, get_property_table
from CINDES.molecule import Molecule
from CINDES.utils.utils import skipper

from skopt.utils import normalize_dimensions
from skopt import Optimizer

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


def GP(run):
    """ This will be a Bayes Optimization: Gaussian Process """
    my_gp = GaussianProcess(run)
    result = my_gp.evolve()
    return result


class GaussianProcess(Algorithm):
    def __init__(self, run):
        self.run = run
        self.param = self.run.gprf
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
        print("array:", self.array)
        space = normalize_dimensions(self.array_joined)
        print("space:", space)
        self.space = space
        return

    def get_Y(self, X, gen):

        # 1. get set of unique confs
        indices = tuple( '_'.join(item) for item in X )
        unique_confs = [ zcon.indtocon(item) for item in set(indices)]

        # 2. make molecules
        mols = [Molecule(conf=conf) for conf in unique_confs]

        # 3. check if already in table
        mols, nnewcalcs, made_pred = evaluate_mols(self.run, mols, self.table, gen, nsite=0)
        self.ncalcs += nnewcalcs
        self.optimum, _ = testmax(self.run, mols)
        print("mols_all:", mols)

        # 5. log new results
        self.table = loggings(mols,
                              self.table,
                              gen,
                              1, 1,
                              made_pred=False,
                              tablename=self.run.tablename,
                              write=self.run.write)
        # 6. set Y
        y_dict = {mol.index: mol.Pvalue for mol in mols}
        Y=[]
        for index in indices:
            Y.append(y_dict[index])
        return Y

    def evolve(self):
        best = None
        self.set_space()

        # algorithms: GP(default), RF, ET, GBRT
        # acq_func: gp_hedge(default), LCB, EI, PI
        # acq_optimizer: sampling(for categorical) lbfgs
        print("self.param:", self.param)

        optimizer = Optimizer(
            dimensions=self.space,
            base_estimator=self.param['algorithm'],
            n_initial_points=self.param['n_initial_points'],
            acq_func=self.param['acq_func'],
            random_state=1287294368,
            acq_optimizer=self.param['acq_optimizer']
            )


        print("optimizer model:", optimizer.base_estimator_)
        print("eta: {}, acq-function: {}, acq-optimizer: {}".format(optimizer.eta, optimizer.acq_func, optimizer.acq_optimizer))

        for gen in range(self.run.maxiter):
            print_title("BATCH-NO: " + str(gen), outline='l', signator="-")

            X = optimizer.ask(n_points=self.param['batchsize'])
            Y = self.get_Y(X, gen=gen)
            X, Y = self.only_finite(X, Y)
            optimizer.tell(X, Y)

            for x, y in zip(X, Y):
                if best is None or best[1] > y:
                    best = ('_'.join(x), y)

            # save model:
            with open('my-optimizer.pkl', 'wb') as f:
                pickle.dump(optimizer, f)

            self.history.append({
                    'count':gen,
                    'p':best[1],
                    'index':best[0],
                    'ncalcs':self.ncalcs
                    })


            logging.info("--- %s seconds ---" % (time.time() - self.run.starttime))

        print("\n\tOptimum:", min(list(zip(optimizer.yi, optimizer.Xi))))
        result = {'p':best[1], 'index':best[0], 'history':self.history, 'ncalcs':self.ncalcs, 'count':self.run.maxiter}
        return result

    @staticmethod
    def only_finite(X, Y):
        print("X, Y before:", X, Y)
        # Y  is converted to np.float array so None values become np.nan values since
        # isfinite cannot handle None's
        X, Y = list(zip(*[ (x,y) for x,y in zip(X,np.array(Y, dtype=np.float)) if np.isfinite(y) ]))
        print("X, Y after:", X, Y)
        return X, Y

