import pickle
import time

import numpy as np
import pandas as pd
pd.set_option('display.width',150)

from sklearn.model_selection import KFold
from scipy.stats import pearsonr

from IO import get_XY, get_X
from CINDES4.utils.utils import processify
from CINDES4.utils.statistics import print_stats

debug=False


class Experiment(object):

    def __init__(self,
                 run,
                 name,
                 table=[],
                 array=[],
                 n_folds=5,
                 retrain=False,
                 getR = False,
                 reoptimize = True,
                 descriptor='BoB',
                 **kwargs):
        """
        Initialize experiment by reading/constructing data.

        params:
            - setting: 'IP' or 'HLG'
            - n_folds: number of folds for splitting training and test data
        """
        np.random.seed(run.seed)
        self.reoptimize = reoptimize
        self.getR = getR
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

    def get_best_hyperparams(self):
        """
        Interface for optimization of hyperparameters of the Experiment

        This function will change self.hparams to best performing hparams

        NB! this function is used for the sklearn based methods. 
         - Gaussian Processes has its own implementation of this function! see gp.py
        
        hyperparameter search with use of the sklearn GridSearchCV function
        """
        from sklearn.model_selection import GridSearchCV
        import time

        # 1. set hyperparamter search
        clf_gs = GridSearchCV( self.get_estimator(),
                               cv= self.n_folds,
                               n_jobs=8,
                               param_grid = self.hparam_grid )

        # 2. set hyper_param set and validation set
        # if X.shape[0] > 400: do 300 for training. else do 75% for training
        ind = np.arange(self.X.shape[0])
        np.random.shuffle(ind)
        if self.X.shape[0] > 400:
            n_opt = 300
        else:
            n_opt = int( 0.75 * self.X.shape[0] )
        X_hyp = self.X[ind[:n_opt],:]
        y_hyp = self.y[ind[:n_opt]]
        X_test= self.X[ind[n_opt:],:]
        y_test= self.y[ind[n_opt:]]
        print "\tn X:", self.X.shape, "n X_hyp:", X_hyp.shape, "n X_test:", X_test.shape

        # 3. Fit the GridSearch
        stime = time.time()
        clf_gs.fit(X_hyp, y_hyp)

        # 4. determine R on validation set:
        self.R = pearsonr( y_test, clf_gs.predict(X_test))[0]

        # 5. print time consumed
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # 6. print results
        print "clf_gs:", clf_gs
        print "best_params_:", clf_gs.best_params_
        print "best_score_:", clf_gs.best_score_
        # print "R**2 on validation set:", self.R
        # print "cv_results_", clf_gs.cv_results_ # too verbose

        # 7. update model.hparams to best ones. 
        self.hparam.update(clf_gs.best_params_)

        return


    def get_model(self, count=0, nsite=0, *args, **kwargs):
        """
        Interface to get the model to use for prediction.
        either:
            - load old model
            - fit new model
        """

        # 1. If self.retrain=False: try to load model. but if not found do nevertheless a training with hparam opt.
        if not self.retrain:
            try:
                self.model = self.load_model(count)
                # it gets the R**2 value from the moment where the model was created.
                print "self.R:", self.R
                return
            except IOError as e:
                print "tried to load model but not found:", e
                print "going to train model:"
                self.X, self.y = get_XY(table=self.table,
                                        descriptor=self.descriptor,
                                        identify = self.run.identify,
                                        array = self.array,
                                        TZmat=self.run.TZmat,
                                        **kwargs )
                self.reoptimize = True

        if self.reoptimize: # sets self.hparam
            self.get_best_hyperparams()
        elif self.getR:
            self.cross_val()
        print "R**2 value is:", self.R
        
        # and always do a refit on total database:
        print "Training for final model..."
        stime = time.time()
        self.model  = self.train(verbose=True, **kwargs)
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # save model:
        self.save_model(count)
        return

    def predict(self, molecules, **kwargs):
        ''' get molecules list 

        NB: for the LinRegOneExperiment this function is overwritten because it uses get_X_1D
        '''
        indices = [ mol.index for mol in molecules ]
        X_pred = get_X(indices, array=self.array, descriptor=self.descriptor, identify=self.run.identify, **self.run.TZmat )
        y_pred = self.test( X_pred, **kwargs)
        if debug: print "y_pred:", y_pred
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
        self.R = stats_df_test.loc['means']['r']

        return stats_df_train.loc['means'], stats_df_test.loc['means']


    def get_fold(self):
        for train_ind, test_ind in KFold(n_splits=self.n_folds, shuffle=True).split(self.X):
            print "shape self.X:", self.X.shape, "     shape self.y:", self.y.shape,
            yield (self.X[train_ind,:], self.y[train_ind], self.X[test_ind,:], self.y[test_ind])



# move later to a util module:

#def print_stats(y, pred):
#    from scipy.stats import pearsonr
#    sq_err = (y - pred)**2
#
#    r = pearsonr(y, pred)
#    mean_err = (np.mean(sq_err), np.var(sq_err))
#    perc = tuple(np.percentile(sq_err, quantile) for quantile in [25, 50, 75])
#
#    #print "\tPearson's R: ", r
#    #print "\tMean error: ", mean_err
#    #print "\tPercentiles: ", perc
#
#    ret = [ r[0], r[1], mean_err[0] ]
#    ret.extend(perc)
#    return ret



