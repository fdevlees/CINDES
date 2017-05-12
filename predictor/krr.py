from sklearn.decomposition.pca import PCA
from sklearn.externals import joblib
from sklearn.kernel_ridge import KernelRidge
from sklearn.utils import resample
import pandas as pd
import numpy as np

from experiment_interface import Experiment

# krr with rbf kernel. C controls simplisity or decision surface. High C will
# try to fit all data and select more support vector. Low C will give a more
# smooth surface. gamma parameter defines how far the influence of a single
# training example reaches, with low values meaning far and high values
# meaning close. low gamma value high bias, high values high variance.

class KernelRidgeExperiment(Experiment):

    def __init__(self, **kwargs):
        """
        In **kwargs:
            - run=run
            - n_folds=n_folds
            - retrain
            - table
        optional also:
            - n_principal_components
        """
        super(KernelRidgeExperiment, self).__init__(**kwargs)

        # set default hyperparam ( super sets self.hparams to dict() ) so beware of order.
        self.hparam = { 'kernel':'rbf',
                        'alpha': 1.e-3,
                        'gamma':1.e-6 }  # gamma parameter is specific for rbf/laplacian kernel

        self.hparam_grid = {'alpha': np.logspace(-10,5,10),
                            'gamma': np.logspace(-10,5,10),
                            'kernel': ['rbf','laplacian'] }
        
        # see if new defaults are given via input
        for key in self.hparam:
            if key in kwargs:
                self.hparam[key] = kwargs[key]
                print "new default hyperparameter:", key, kwargs[key]

        return

    def train(self, X=None, y=None, verbose=True, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        krr_rbf = self.get_estimator()
        print "Fitting...",
        krr_rbf.fit(X, y)

        if verbose: print "\tLearned model: ", krr_rbf

        return krr_rbf

    def get_estimator(self):
        krr_rbf = KernelRidge(kernel='rbf',
                      alpha = self.hparam['alpha'],
                      gamma = self.hparam['gamma'],
                      )
        return krr_rbf
                       

    def test(self, X, model=None, **kwargs ):
        if model is None: model=self.model
        return model.predict(X).flatten()

    def save_model(self, count, model=None):
        if model is None: model=self.model
        model.R = self.R

        modelname = '{}_{}.pkl'.format(self.name, count)
        joblib.dump(model,modelname)
        return

    def load_model(self,count):
        modelname = '{}_{}.pkl'.format(self.name, count)
        model = joblib.load(modelname)
        self.R = model.R
        return model



    def get_best_hyperparams_old(self):
        ''' hyperparameter search '''
        # save R**2, MAE and percentiles of each fold to a row in a dataframe.
        stats_df_train = pd.DataFrame(columns=('C', 'r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))
        stats_df_test =  pd.DataFrame(columns=('C', 'r','p-value','mae','perc_25', 'perc_50', 'perc_75' ))

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

        # get hyperparameter with highest R**2
        C = stats_df_test.loc[stats_df_test['r'].idxmax()]['C']
        print "self.C:", C
        self.hparam['C'] = C

        #raise SystemExit('stop')
        return

class KernelRidgeWithPCAExperiment(KernelRidgeExperiment):

    def __init__(self, n_principal_components, **kwargs):
        super(KernelRidgeWithPCAExperiment, self).__init__(**kwargs)
        self.n_principal_components = n_principal_components

    def get_XY(self, **kwargs):
        X,y = super(KernelRidgeWithPCAExperiment, self).get_XY(**kwargs)
        return self.do_PCA(X), y

    def get_X(self, *args, **kwargs):
        X = super(KernelRidgeWithPCAExperiment, self).get_X(*args, **kwargs)
        return self.do_PCA(X, fit=False)

    def do_PCA(self, X, fit=True):
        if fit:
            F = PCA(self.n_principal_components)
            F.fit(X)
            print "\tLeast explained variance:", F.explained_variance_[-1]
            self.F = F
        X_F = self.F.transform(X)
        print "\tDimensionality reduction: ", X_F.shape
        return X_F

#    def train(self, X=None, y=None, **kwargs):
#        if X is None: X=self.X
#        if y is None: y=self.y
#        # Dimensionality reduction
#        F = PCA(self.n_principal_components)
#        F.fit(X)
#        X_F = F.transform(X)
#
#        print "\tLeast explained variance:", F.explained_variance_[-1]
#        print "\tDimensionality reduction: ", X_F.shape
#
#        # Nearest neighbor
#        krr = super(KernelRidgeWithPCAExperiment, self).train(X_F, y, **kwargs)
#        self.F = F
#
#        return krr
#
#    def test(self, X, model=None, **kwargs):
#        if model is None: model=self.model
#        krr = model
#        X_F = self.F.transform(X)
#        return super(KernelRidgeWithPCAExperiment, self).test(X_F, krr, **kwargs)

    def save_model(self, count, model=None):
        if model is None: model=self.model
        model.R = self.R

        modelname = '{}_pca_{}.pkl'.format(self.name,count)
        joblib.dump( (model, self.F) ,modelname)
        return

    def load_model(self,count):
        modelname = '{}_pca_{}.pkl'.format(self.name, count)
        (model, self.F) = joblib.load(modelname)
        self.R = model.R
        return model



if __name__=="__main__":
    import pickle
    print "SVR:"
    class Run:
        pass
    retrain = True
    n_folds = 3
    run = Run()
    run.identify = 'ada_'
    table = pickle.load(open('tablebin','rb'))
    regressor = KernelRidgeExperiment(table=table, retrain=True, n_folds=n_folds, run=run, identify='ada_',
                                        descriptor= '1DL' )

    print "regressor:", regressor



