
import pickle
import time

import numpy as np
import pandas as pd
pd.set_option('display.width',150)

from sklearn.model_selection import KFold, train_test_split
from scipy.stats import pearsonr

from IO import get_XY, get_X
from CINDES.utils.utils import processify
from CINDES.utils.statistics import print_stats
from CINDES.utils.writings import sprint, log_io, dump

debug=False


class Experiment(object):

    def __init__(self,
                 run=None,
                 name='default_estimator',
                 table=[],
                 array=[],
                 n_folds=5,
                 retrain=False,
                 getR = False,
                 reoptimize = True,
                 multiple = False,
                 descriptor='BoB',
                 weights= False,
                 plots=[],
                 **kwargs):
        """
        Initialize experiment by reading/constructing data.

        params:
            - setting: 'IP' or 'HLG'
            - n_folds: number of folds for splitting training and test data
        """
        if run: np.random.seed(run.seed)
        self.reoptimize = reoptimize
        self.getR = getR
        self.retrain = retrain
        self.run = run
        self.descriptor = descriptor
        self.name = name
        self.n_folds = n_folds
        self.array = array # only for 'int' descriptor
        self.table = table
        self.plots = plots
        self.retrained = False
        self.multiple = multiple
        self.R = None
        self.weights = weights
        if debug: print "self.multiple:", self.multiple

        if self.retrain:
            #self.X, self.y = read_BoB_data(setting, '../data')
            self.X, self.y = get_XY(table=self.table, descriptor=self.descriptor, identify = self.run.identify, array=self.array, TZmat=run.TZmat, **kwargs)

    def get_XY(self, **kwargs):
        X, y = get_XY(table=self.table,
                                descriptor=self.descriptor,
                                identify = self.run.identify,
                                array = self.array,
                                TZmat=self.run.TZmat,
                                **kwargs )
        return X, y

    def get_X(self, indices, **kwargs):
        X_pred = get_X(indices, array=self.array, descriptor=self.descriptor, identify=self.run.identify, **self.run.TZmat )
        return X_pred

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

        # 1. If self.retrain=False: try to load model. but if not found do nevertheless a training with hparam opt.
        if not self.retrain and not self.multiple:
            try:
                self.model = self.load_model(count=count)
                print "     LOAD succesful!"
                print "     self.model:", self.model
                # it gets the R**2 value from the moment where the model was created.
                print "self.R:", self.R
                return
            except IOError as e:
                print "tried to load model but not found:", e
                print "going to train model:"
                self.X, self.y = self.get_XY(**kwargs)
                self.reoptimize = True

        if self.multiple:
            self.do_multiple(**kwargs)
        elif self.reoptimize: # sets self.hparam
            print "self.hparam before h_opt:", self.hparam
            clf_gs = self.get_best_hyperparams( split=True)
            print "self.hparam after h_opt:", self.hparam
        elif self.getR:
            self.cross_val()
        print "R**2 value is:", self.R

        # and always do a refit on total database:
        if not self.multiple:
            if True:
                print "Training for final model..."
                stime = time.time()
                self.model  = self.train(verbose=True, **kwargs)
                time_to_fit = time.time() - stime
                print "\tTime to fit: ", time_to_fit, ' s'
            else:
                self.model = clf_gs

            # save model:
            self.save_model(count=count)

        return

    def get_best_hyperparams(self, split=True):
        """
        Interface for optimization of hyperparameters of the Experiment

        This function will change self.hparams to best performing hparams

        NB! this function is used for the sklearn based methods.
         - Gaussian Processes has its own implementation of this function! see gp.py 
        hyperparameter search with use of the sklearn GridSearchCV function
        """
        from sklearn.model_selection import GridSearchCV
        import time
        stime = time.time()

        # 1. set hyperparamter search
        clf_gs = GridSearchCV( self.get_estimator(),
                               cv= self.n_folds,
                               n_jobs=8,
                               param_grid = self.hparam_grid )

        # 2. set hyper_param set and validation set
        # if X.shape[0] > 400: do 300 for training. else do 75% for training

        if split: # use a validation set. after the hyper_opt. This is used in normal run to decide based on R**2 value
            ind = np.arange(self.X.shape[0])
            np.random.shuffle(ind)
            if self.X.shape[0] > 600: # use only 300 or if nX<400 only 75% of items. 
                n_opt = 300
            else:
                n_opt = int( 0.75 * self.X.shape[0] )
            X_hyp, X_test, y_hyp, y_test = train_test_split( self.X, self.y, train_size=n_opt, random_state = self.run.seed )
            #X_hyp = self.X[ind[:n_opt],:]
            #y_hyp = self.y[ind[:n_opt]]
            #X_test= self.X[ind[n_opt:],:]
            #y_test= self.y[ind[n_opt:]]
            print "\tn X:", self.X.shape, "n X_hyp:", X_hyp.shape, "n X_test:", X_test.shape

            # 3. Fit the GridSearch
            gs_results = clf_gs.fit(X_hyp, y_hyp)

            # 4. determine R on validation set:
            y_pred = clf_gs.predict(X_test)
            y_pred_hyp = clf_gs.predict(X_hyp)
            self.R = pearsonr( y_test, y_pred)[0]
            self.plot1 = np.c_[ y_test, y_pred ]
            self.plot2 = np.c_[ y_hyp, y_pred_hyp ]
            for i in range(10):
                try:print y_test[i], y_pred[i]
                except IndexError: pass
            print "8"*30

        else:   # no split. because in do_multiple there is already another validation set specified. 
            clf_gs.fit(self.X, self.y)

        # 5. print time consumed
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # 6. print results
        print "clf_gs:", clf_gs
        print "best_params_:", clf_gs.best_params_
        print "best_score_:", clf_gs.best_score_
        self.retrained = True
        # print "R**2 on validation set:", self.R
        # print "cv_results_", clf_gs.cv_results_ # too verbose

        # 7. update model.hparams to best ones. 
        self.hparam.update(clf_gs.best_params_)

        return clf_gs

    def do_multiple(self, n_max_train=None, **kwargs):
        results = dict()

        # 1. get X and y
        X, y = get_XY(table=self.table,
                      descriptor=self.descriptor,
                      identify = self.run.identify,
                      array = self.array,
                      TZmat=self.run.TZmat,
                      **kwargs )
        n = y.shape[0]

        # determine the numbers based on the datasize or based just on powers of 2 from 8 to 1024
        if n<300:
            ns = [16,32,64] 
        elif True: # maksimum n_totaal - 200 samples
            ns = np.logspace( 4, np.log2(n-200), base=2, num=10, dtype=int)
        else: # maksimum 1024 samples
            ns = np.logspace(3,10,base=2,num=8,dtype=int)
        if debug: print "ns:", ns

        # 2. split in data_train and data_test
        for N in ns: # for every n in Ns
            N_bu = N

            # take a maximum number of training samples
            if n_max_train:
                if N > n_max_train:
                    print "number of training data set from {:d} to {:d}".format(N, n_max_train)
                    N=n_max_train

            # split in train / test
            print "\n    N:", N
            if True:
                result = []
                #for fold, (train_ind, test_ind) in enumerate( KFold(n_splits=self.n_folds, shuffle=True).split(X) ):
                #for fold, (X_train, X_test, y_train, y_test) 
                for fold in range(3):
                    print "\n        FOLD:", fold

                    X_train, X_test, y_train, y_test = train_test_split( X, y, train_size=N, random_state = self.run.seed + fold )

                    #X_train, X_test = X[train_ind], X[test_ind]
                    #y_train, y_test = y[train_ind], y[test_ind]
                    self.X = X_train
                    self.y = y_train

                    # 3. hparam opt on data_train
                    best_estimator = self.get_best_hyperparams(split=False)
                    print "best_estimator:", best_estimator

                    # 4. validate model on data_test
                    y_test_pred = best_estimator.predict(X_test)
                    y_train_pred = best_estimator.predict(X_train)

                    # 5. return R and y_train_pred / y_test_pred
                    R_test = pearsonr( y_test, y_test_pred)[0]
                    R_train = pearsonr( y_train, y_train_pred)[0]
                    print "    R_test :", R_test
                    print "    R_train:", R_train

                    result.append({ 'R_test' :R_test  , 'R_train'     :R_train,
                                    'y_test' :y_test  , 'y_test_pred' :y_test_pred,
                                    'y_train':y_train , 'y_train_pred':y_train_pred })
                results[N_bu]=result
                if False: # if i want to let it write intermediate data
                    self.save_json(results, extra_identifier=str(N_bu))


            else:
                ind = np.arange(n)
                np.random.shuffle(ind)
                X_train = X[ind[:N],:]
                y_train = y[ind[:N]]
                X_test= X[ind[N:],:]
                y_test= y[ind[N:]]
                self.X = X_train
                self.y = y_train
             
                # 3. hparam opt on data_train
                best_estimator = self.get_best_hyperparams(split=False)
                print "best_estimator:", best_estimator
             
                # 4. validate model on data_test
                y_test_pred = best_estimator.predict(X_test)
                y_train_pred = best_estimator.predict(X_train)
             
                # 5. return R and y_train_pred / y_test_pred
                R_test = pearsonr( y_test, y_test_pred)[0]
                R_train = pearsonr( y_train, y_train_pred)[0]
                print "    R_test :", R_test
                print "    R_train:", R_train
             
                # 6. save in a data format Not anymore: Note that an integer is used as dict key!
                results[N] = [{ 'R_test': R_test , 'R_train':R_train,
                               'y_test' :y_test , 'y_test_pred' :y_test_pred,
                               'y_train':y_train, 'y_train_pred':y_train_pred }]

        print "results:", dump(results, n=10)
        self.plot3 = results

        # jump results as json formatted file
        self.save_json(results)

        return

    def do_predict(self, molecules, MC=False, **kwargs):
        ''' get molecules list 

        NB: for the LinRegOneExperiment this function is overwritten because it uses get_X_1D
        '''
        indices = [ mol.index for mol in molecules ]
        X_pred = self.get_X(indices)
        #sprint(10, X_pred)
        y_pred = self.test( X_pred, **kwargs)

        if not MC:
            # test if numbers are different enough.
            if len(set(y_pred.round(2)))<3 and len(y_pred)>10:
                self.R = -0.5
                print "TOO MUCH SIMILAR PREDICTIONS> OVERFITTING?! > self.R set to -0.5"

        if debug and not MC: print "y_pred:", y_pred
        for molecule, y in zip(molecules, y_pred):
            molecule.predictions[self.name] = y
            #print molecule, y
        return y_pred

    def cross_val(self, write_log=False, plot=True, **kwargs):
        # save R**2, MAE and percentiles of each fold to a row in a dataframe.
        stats_df_train = pd.DataFrame(columns=('r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))
        stats_df_test =  pd.DataFrame(columns=('r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))

        # for plotting:
        r_test = []
        r_test_pred = []
        r_train = []
        r_train_pred = []

        # Randomly split training and py data
        for fold, (X_train, y_train, X_test, y_test) in enumerate(self.get_fold()):
            print "FOLD ", fold
            
            # this 5 lines are also in get_model without crosval
            print "Training ...",
            r_test.extend(y_test)
            r_train.extend(y_train)
            stime = time.time()
            model = self.train(X_train, y_train, **kwargs)
            time_to_fit = time.time() - stime
            print "\tTime to fit: ", time_to_fit, ' s'
            
            print "Testing ...",
            stime = time.time()
            y_train_pred = self.test(X_train, model)
            y_test_pred = self.test(X_test, model)
            r_test_pred.extend(y_test_pred)
            r_train_pred.extend(y_train_pred)
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

        if plot:
            self.plot1 = np.c_[ r_test, r_test_pred ]
            self.plot2 = np.c_[ r_train, r_train_pred ]
            self.retrained = True

        return stats_df_train.loc['means'], stats_df_test.loc['means']


    def get_fold(self):
        for train_ind, test_ind in KFold(n_splits=self.n_folds, shuffle=True).split(self.X):
            print "shape self.X:", self.X.shape, "     shape self.y:", self.y.shape,
            yield (self.X[train_ind,:], self.y[train_ind], self.X[test_ind,:], self.y[test_ind])

    def save_json(self, results, extra_identifier=None):
        import json


        json_results = jsonify(results)
        if extra_identifier:
            name_json_out = 'results_{}_{}'.format(self.name, extra_identifier)
        else:
            name_json_out = 'results_{}.json'.format(self.name)

        with open(name_json_out,'w') as f:
            try:
                json.dump(json_results, f, sort_keys=True, separators=(',', ':'), indent=0)
            except TypeError:
                print "json error"
                raise
        return

    def get_weights(self,n):
        n_higher = int( 0.20 * n )
        higher = np.full(n_higher, 20.)
        lower =  np.ones( n-n_higher)
        sample_weights = np.concatenate( (higher, lower), axis=0)
        print "    sample weights:\n", sample_weights

        return sample_weights

def jsonify(data):
    json_data = dict()
    for key, value in data.iteritems():
        if isinstance(value, list):
            value = [ jsonify(item) if isinstance(item, dict) else item for item in value ]
        if isinstance(value, dict):
            value = jsonify(value)
        if isinstance(key, int):
            key = str(key)
        if type(value).__module__=='numpy':
            value = value.tolist()
        json_data[key] = value
    return json_data


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



