import pickle
import time

import numpy as np
import pandas as pd
pd.set_option('display.width',150)

from sklearn.model_selection import KFold

from IO import get_XY, get_X

debug=True


class Experiment(object):

    def __init__(self,
                 run,
                 name,
                 table=[],
                 array=[],
                 n_folds=5,
                 retrain=True,
                 descriptor='BoB',
                 **kwargs):
        """
        Initialize experiment by reading/constructing data.

        params:
            - setting: 'IP' or 'HLG'
            - n_folds: number of folds for splitting training and test data
        """
        np.random.seed(run.seed)
        self.reoptimize = True
        self.getR = True
        self.retrain = retrain
        self.run = run
        self.descriptor = descriptor
        self.name = name
        self.n_folds = n_folds
        self.array = array # only for 'int' descriptor
        self.table = table

        if self.retrain:
            #self.X, self.y = read_BoB_data(setting, '../data')
            self.X, self.y = get_XY(table=self.table, descriptor=self.descriptor, identify = self.run.identify, array=self.array, TZmat=run.TZmat, **kwargs)

    def train(self, X, y, **kwargs):
        """
        Interface for training.

        Params:
            - X: training inputs [(n, D) array]
            - y: training outputs [(n,) array]
        Returns:
            - model, log
        """
        pass

    def test(self, X, y, model, **kwargs):
        """
        Interface for training.

        Params:
            - X: py inputs [(n, D) array]
            - y: py outputs [(n,) array]
            - model
        Returns:
            - predicted outputs
        """
        pass

    def load_model(self, count):
        """
        Interface for loading the model

        Returns: 
            - model
        """
        pass

    def save_model(self, model):
        """
        Interface for saving the model
        
        Params:
            - model
        """
        pass

    def get_model(self, count=0, nsite=0, *args, **kwargs):
        """
        Interface to get the model to use for prediction.
        either:
            - load old model
            - fit new model
        """
        if not self.retrain:
            try:
                self.model = self.load_model(count)
                return
            except IOError as e:
                print "tried to load model but not found:", e
                print "going to train model:"
                self.X, self.y = get_XY(table=self.table, descriptor=self.descriptor, identify = self.run.identify, array = self.array, TZmat=self.run.TZmat, **kwargs)

        if self.reoptimize:
            self.get_best_hyperparams()
        elif self.getR:
            self.cross_val()
        
        # and always do a refit on total database:
        print "Training for final model..."
        stime = time.time()
        self.model  = self.train(verbose=True, **kwargs)
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # save model:
        self.save_model(count)
        return

    def predict(self, molecules):
        ''' get molecules list 

        NB: for the LinRegOneExperiment this function is overwritten because it uses get_X_1D
        '''
        indices = [ mol.index for mol in molecules ]
        X_pred = get_X(indices, array=self.array, descriptor=self.descriptor, identify=self.run.identify, **self.run.TZmat )
        y_pred = self.test( X_pred)
        print "y_pred:", y_pred
        for molecule, y in zip(molecules, y_pred):
            molecule.predictions[self.name] = y
            print molecule, y
        return y_pred

    def cross_val(self, write_log=False, **kwargs):
        # save R**2, MAE and percentiles of each fold to a row in a dataframe.
        stats_df_train = pd.DataFrame(columns=('r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))
        stats_df_test =  pd.DataFrame(columns=('r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))

        # Randomly split training and py data
        for fold, (X_train, y_train, X_test, y_test) in enumerate(self.get_fold()):
            print "FOLD ", fold
            
            # this 5 lines are also in get_model without crosval
            print "Training ...",
            stime = time.time()
            model = self.train(X_train, y_train, **kwargs)
            time_to_fit = time.time() - stime
            print "\tTime to fit: ", time_to_fit, ' s'
            
            print "Testing ...",
            stime = time.time()
            y_train_pred = self.test(X_train, model)
            y_test_pred = self.test(X_test, model)
            time_to_fit = time.time() - stime
            print "\tTime to predict: ", time_to_fit, ' s'
     
            if debug:
                print "y_train:", y_train
                print "y_train_pred:", y_train_pred
                print "y_test:", y_test
                print "y_test_pred:", y_test_pred


            #print "Training accuracy:"
            #print "Testing accuracy:"
            stats_df_train.loc[fold] = print_stats(y_train, y_train_pred)
            stats_df_test.loc[fold] = print_stats(y_test, y_test_pred)

            #if write_log:
            #    log = {"seed": seed, "fold": fold, "model": log_model, "x_train": X_train, "y_train": y_train, "x_test": X_test, "y_test": y_test, "y_pred": y_test_pred}
            #    with open(''.join(["log_", self.setting, "_", str(seed) , "_", str(fold), ".pic"]), 'wb') as f:
            #        pickleA.dump(log, f)

        stats_df_train.loc['means']= stats_df_train.mean()
        stats_df_test.loc['means'] = stats_df_test.mean()
        print "train statistics:\n", stats_df_train
        print " test statistics:\n", stats_df_test

        return stats_df_train.loc['means'], stats_df_test.loc['means']

    def get_best_hyperparams(self):
        # save R**2, MAE and percentiles of each fold to a row in a dataframe.
        stats_df_train = pd.DataFrame(columns=('C', 'r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))
        stats_df_test =  pd.DataFrame(columns=('C', 'r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))

        C_s = np.logspace(-10, 10, 21)
        gamma_s = np.logspace( -10,10,21)
        #for i, C in enumerate(C_s):
        for i, gamma in enumerate(gamma_s):
            #df_train, df_test = self.cross_val(C=C)
            df_train, df_test = self.cross_val(gamma=gamma)

            stats_df_train.loc[i] = df_train
            stats_df_test.loc[i] = df_test
            stats_df_train['C'][i] = gamma
            stats_df_test['C'][i] = gamma

        print "end of hyper opt:"
        print "stats_df_train:\n", stats_df_train
        print "stats_df_test:\n", stats_df_test


        raise SystemExit('stop')
        return

    def get_fold(self):
        for train_ind, test_ind in KFold(n_splits=self.n_folds, shuffle=True).split(self.X):
            print "shape self.X:", self.X.shape, "     shape self.y:", self.y.shape,
            yield (self.X[train_ind,:], self.y[train_ind], self.X[test_ind,:], self.y[test_ind])



# move later to a util module:

def print_stats(y, pred):
    from scipy.stats import pearsonr
    sq_err = (y - pred)**2

    r = pearsonr(y, pred)
    mean_err = (np.mean(sq_err), np.var(sq_err))
    perc = tuple(np.percentile(sq_err, quantile) for quantile in [25, 50, 75])

    #print "\tPearson's R: ", r
    #print "\tMean error: ", mean_err
    #print "\tPercentiles: ", perc

    ret = [ r[0], r[1], mean_err[0] ]
    ret.extend(perc)
    return ret



