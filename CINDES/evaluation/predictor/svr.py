from sklearn.decomposition.pca import PCA
from sklearn.externals import joblib
from sklearn.svm import SVR
from sklearn.utils import resample
import pandas as pd
import numpy as np

from .experiment_interface import Experiment

# svr with rbf kernel. C controls simplisity or decision surface. High C will
# try to fit all data and select more support vector. Low C will give a more
# smooth surface. gamma parameter defines how far the influence of a single
# training example reaches, with low values meaning far and high values
# meaning close. low gamma value high bias, high values high variance.

class SupportVectorExperiment(Experiment):

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
        super(SupportVectorExperiment, self).__init__(**kwargs)

        # set default hyperparam ( super sets self.hparams to dict() ) so beware of order.
        self.hparam = { 'C':1.e4,
                        'gamma':1.e-6,
                        'tol':0.001  }

        self.hparam_grid = {'C': np.logspace(-5,5,3),
                            'gamma': np.logspace(-2,2,3) }
        
        #C_s = np.logspace(-10, 10, 3)
        #gamma_s = np.logspace( -10,10,3)

        # see if new defaults are given via input
        for key in self.hparam:
            if key in kwargs:
                self.hparam[key] = kwargs[key]
                print("new default hyperparameter:", key, kwargs[key])

        return

    def train(self, X=None, y=None, verbose=True, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        svr_rbf = self.get_estimator()
        print("Fitting...", end=' ')
        svr_rbf.fit(X, y)

        if verbose: print("\tLearned model: ", svr_rbf)

        return svr_rbf

    def get_estimator(self):
        svr_rbf = SVR(kernel='rbf',
                      C= self.hparam['C'],
                      gamma= self.hparam['gamma'],
                      tol = self.hparam['tol'],
                      cache_size=200 # number of megabytes
                      )
        return svr_rbf
                       

    def test(self, X, model=None, **kwargs):
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
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        from sklearn.model_selection import GridSearchCV
        import time

        # 1. set hyperparamter search
        svr = GridSearchCV( self.get_estimator(),
                            cv= self.n_folds,
                            n_jobs=8,
                            param_grid = self.hparam_grid )

        # 2. do search on dataset
        if True:
            n_train = 300
            #X = self.X[:n_train]
            #y = self.y[:n_train]
            X, y = resample(self.X, self.y, n_samples=n_train)
            print("restricted hparamopt to only {} samples".format(n_train))
        else:
            X = self.X
            y = self.y
        stime = time.time()
        svr.fit(X, y)
        time_to_fit = time.time() - stime
        print("\tTime to fit: ", time_to_fit, ' s')

        # 3. print results
        print("svr:", svr)
        print("n svr.best_estimator_.support_", len(svr.best_estimator_.support_))
        print("best_params_:", svr.best_params_)
        print("best_score_:", svr.best_score_)
        self.R = svr.best_score_
        #print "cv_results_", svr.cv_results_ # too verbose

        # 4. update model.hparams to best ones. 
        self.hparam.update(svr.best_params_)

        return

    def get_best_hyperparams_older(self):
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

        print("end of hyper opt:")
        print("stats_df_train:\n", stats_df_train)
        print("stats_df_test:\n", stats_df_test)

        # get hyperparameter with highest R**2
        C = stats_df_test.loc[stats_df_test['r'].idxmax()]['C']
        print("self.C:", C)
        self.hparam['C'] = C

        #raise SystemExit('stop')
        return


class SupportVectorWithPCAExperiment(SupportVectorExperiment):

    def __init__(self, n_principal_components, **kwargs):
        super(SupportVectorWithPCAExperiment, self).__init__(**kwargs)
        self.n_principal_components = n_principal_components

    def train(self, X=None, y=None, **kwargs):
        if X is None: X=self.X
        if y is None: y=self.y
        # Dimensionality reduction
        F = PCA(self.n_principal_components)
        F.fit(X)
        X_F = F.transform(X)

        print("\tLeast explained variance:", F.explained_variance_[-1])
        print("\tDimensionality reduction: ", X_F.shape)

        # Nearest neighbor
        svr = super(SupportVectorWithPCAExperiment, self).train(X_F, y, **kwargs)
        self.F = F

        return svr

    def test(self, X, model=None, **kwargs):
        if model is None: model=self.model
        svr = model
        X_F = self.F.transform(X)
        return super(SupportVectorWithPCAExperiment, self).test(X_F, svr, **kwargs)

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
    print("SVR:")
    class Run:
        pass
    retrain = True
    n_folds = 3
    run = Run()
    run.identify = 'ada_'
    table = pickle.load(open('tablebin','rb'))
    regressor = SupportVectorExperiment(table=table, retrain=True, n_folds=n_folds, run=run, identify='ada_',
                                        descriptor= '1DL' )

    print("regressor:", regressor)



