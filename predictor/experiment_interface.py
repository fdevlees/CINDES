import pickle
import time

import numpy as np
from sklearn.model_selection import KFold

from util.io import read_BoB_data, print_stats

class Experiment(object):
    
    def __init__(self, setting, n_folds):
        """
        Initialize experiment by reading BoB data.
        
        params:
            - setting: 'IP' or 'HLG'
            - n_folds: number of folds for splitting training and test data
        """
        self.setting = setting
        self.X, self.y = read_BoB_data(setting, '../data')
        self.n_folds = n_folds

    def train(self, X, y):
        """
        Interface for training.
        
        Params:
            - X: training inputs [(n, D) array]
            - y: training outputs [(n,) array]
        Returns:
            - model, log
        """
        pass
    
    def test(self, X, y, model):
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
    
    def run(self, seed, write_log=False):
        np.random.seed(seed)
        
        # Randomly split training and py data
        for fold, (X_train, y_train, X_test, y_test) in enumerate(self.get_fold()):
            print "FOLD ", fold
        
            print "Training ..."
            stime = time.time()
            model, log_model = self.train(X_train, y_train)
            time_to_fit = time.time() - stime
            print "\tTime to fit: ", time_to_fit, ' s'
            
            print "Testing ..."
            stime = time.time()
            y_train_pred = self.test(X_train, model)
            y_test_pred = self.test(X_test, model)
            time_to_fit = time.time() - stime
            print "\tTime to predict: ", time_to_fit, ' s'
            
            print "Training accuracy:"
            print_stats(y_train, y_train_pred)
            print "Testing accuracy:"
            print_stats(y_test, y_test_pred)

            if write_log:
                log = {"seed": seed, "fold": fold, "model": log_model, "x_train": X_train, "y_train": y_train, "x_test": X_test, "y_test": y_test, "y_pred": y_test_pred}
                with open(''.join(["log_", self.setting, "_", str(seed) , "_", str(fold), ".pic"]), 'wb') as f:
                    pickle.dump(log, f)


    
    def get_fold(self):
        for train_ind, test_ind in KFold(n_splits=self.n_folds, shuffle=True).split(self.X):
            print self.X.shape
            print self.y.shape
            yield (self.X[train_ind,:], self.y[train_ind], self.X[test_ind,:], self.y[test_ind])

