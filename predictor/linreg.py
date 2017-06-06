from sklearn.decomposition.pca import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.externals import joblib
from sklearn import linear_model
from sklearn.utils import resample

import numpy as np

from experiment_interface import Experiment
from descriptor import get_X_1D

class LinRegOneExperiment(Experiment):

    def __init__(self, subtype='ridge', **kwargs):
        """
        In **kwargs:
            - run=run
            - n_folds=n_folds
            - retrain
            - table
        optional also:
            - n_principal_components
        """
        super(LinRegOneExperiment, self).__init__(**kwargs)

        self.subtype = subtype

        self.hparam = { 'alpha':1e4,
                        'tol': 0.001,
                        'intercept':True }
        self.hparam_grid = {'alpha': np.logspace(-5,5,10) }

        for key in self.hparam:
            if key in kwargs:
                self.hparam[key] = kwargs[key]
                print "new default hyperparameter:", key, kwargs[key]
        return


    def get_estimator(self, **kwargs):
        alpha = self.hparam['alpha']
        intercept = self.hparam['intercept']
        tol = self.hparam['tol']
        subtype = self.subtype

        if subtype in ['linreg','ols']:
            clf = linear_model.LinearRegression(fit_intercept=intercept,
                                                copy_X=True)
        elif subtype in ['ridge']:
            clf = linear_model.Ridge(alpha=alpha,
                                     fit_intercept=intercept,
                                     tol=tol,
                                     solver='auto',
                                     copy_X=True)
        elif subtype in ['ridgecv']:
            clf = linear_model.RidgeCV(alphas=alpha,
                                       fit_intercept=intercept,
                                       store_cv_values=True)
        elif subtype in ['lasso']:
            clf = linear_model.Lasso(alpha=alpha,
                                     fit_intercept=intercept,
                                     tol=tol)
        elif subtype in ['ElasticNet']:
            l1_ratio = 0.1 #default 0.5
            clf = linear_model.ElasticNet(alpha=alpha,
                                          l1_ratio=l1_ratio,
                                          fit_intercept=intercept,
                                          tol=tol)
        return clf

    def train(self, X=None, y=None, verbose=False, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        clf = self.get_estimator(**kwargs)

        clf.fit(X,y)

        #print clf.coef_

        if verbose: print "\tLearned model: ", clf

        return clf

    def test(self, X, model=None, **kwargs):
        if model is None: model=self.model
        y_pred = model.predict(X)
        #print y_pred, model, model.coef_
        return y_pred.flatten()

    def predict_old(self, molecules, **kwargs):
        ''' get molecules list '''
        # same
        indices = [ mol.index for mol in molecules ]

        # different
        #X_pred = get_X_1D(indices=indices, descriptor=self.descriptor, identify=self.run.identify)
        X_pred = get_X(indices, array=self.array, descriptor=self.descriptor, identify=self.run.identify, **self.run.TZmat )

        # same
        y_pred = self.test( X_pred, **kwargs)
        if debug: print "y_pred:", y_pred
        for molecule, y in zip(molecules, y_pred):
            molecule.predictions[self.name] = y
            print molecule, y
        return y_pred

    def save_model(self, count, model=None):
        if model is None: model=self.model
        model.R = self.R

        try:
            model.plot1 = self.plot1
            model.plot2 = self.plot2
        except AttributeError:
            pass

        modelname = '{}_{}.pkl'.format(self.name, count)
        joblib.dump(model,modelname)
        return

    def load_model(self, count):
        modelname = '{}_{}.pkl'.format(self.name, count)
        model = joblib.load(modelname)
        self.R = model.R

        try:
            self.plot1 = model.plot1
            self.plot2 = model.plot2
        except AttributeError:
            pass

        return model

    def get_best_hyperparams_old(self):
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        from sklearn.model_selection import GridSearchCV
        import time

        # 1. set hyperparamter search
        clf = GridSearchCV( self.get_estimator(),
                            cv= self.n_folds,
                            n_jobs=8,
                            param_grid = self.hparam_grid )

        # 2. do search on dataset
        if True:
            n_train = 300
            #X = self.X[:n_train]
            #y = self.y[:n_train]
            X, y = resample(self.X, self.y, n_samples=n_train)
            print "restricted hparamopt to only {} samples".format(n_train)
        else:
            X = self.X
            y = self.y
        stime = time.time()
        clf.fit(X, y)
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # 3. print results
        print "clf:", clf
        #print "n clf.best_estimator_.support_", len(clf.best_estimator_.support_)
        print "best_params_:", clf.best_params_
        print "best_score_:", clf.best_score_
        self.R = clf.best_score_
        # print "cv_results_", clf.cv_results_ # too verbose

        # 4. update model.hparams to best ones. 
        self.hparam.update(clf.best_params_)

        return



class LinRegOneWithPCAExperiment(LinRegOneExperiment):

    def __init__(self, n_principal_components=50, **kwargs):
        self.F = None
        super(LinRegOneWithPCAExperiment, self).__init__(**kwargs)
        self.n_principal_components = n_principal_components
        #self.hparam_grid = {'linreg__alpha': np.logspace(-10,5,15) }

    def get_XY(self, **kwargs):
        X,y = super(LinRegOneWithPCAExperiment, self).get_XY(**kwargs)
        return self.do_PCA(X), y

    def get_X(self, *args, **kwargs):
        X = super(LinRegOneWithPCAExperiment, self).get_X(*args, **kwargs)
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
#
#        # Nearest neighbor
#        krr = super(LinRegOneWithPCAExperiment, self).train(X_F, y, **kwargs)
#        self.F = F
#        return krr

#    def test(self, X, model=None, **kwargs):
#        if model is None: model=self.model
#        X_F = self.F.transform(X)
#        return super(LinRegOneWithPCAExperiment, self).test(X_F, model)

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

#    def get_estimator(self, **kwargs):
#        clf = super(LinRegOneWithPCAExperiment, self).get_estimator()
#        pca = PCA(self.n_principal_components)
#        from sklearn.pipeline import Pipeline
#        pipe = Pipeline( steps=[ ('pca', pca), ('linreg', clf)])
#        return pipe

