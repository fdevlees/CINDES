from sklearn.decomposition.pca import PCA
from sklearn.neighbors import NearestNeighbors, KNeighborsRegressor
from sklearn.externals import joblib
from sklearn.utils import resample
import numpy as np

from experiment_interface import Experiment

supervised = True

class NearestNeighborExperiment(Experiment):

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
        super(NearestNeighborExperiment, self).__init__(**kwargs)

        self.hparam = { "n_neighbors":1,
                        "metric": "minkowski",
                        "weights": "uniform"}

        self.hparam_grid = { "n_neighbors": [ 1, 2, 3, 4, 5, 7, 9 ],
                             "metric": ['minkowski', 'euclidean'],
                             "weights": ['uniform', 'distance'] }

        # see if new defaults are given via input
        for key in self.hparam:
            if key in kwargs:
                self.hparam[key] = kwargs[key]
                print "new default hyperparameter:", key, kwargs[key]

        return

    def get_estimator(self,**kwargs):
        if supervised:
            estimator = KNeighborsRegressor(n_neighbors=self.hparam['n_neighbors'],
                                            metric=self.hparam['metric'])
        else:
            estimator = NearestNeighbors(n_neighbors=self.hparam['n_neighbors'],
                                         metric=self.hparam['metric'])
        return estimator

    def train(self, X=None, y=None, verbose=False, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y
        
        if supervised:
            NN = self.get_estimator(**kwargs).fit(X,y)
        else:
            NN = self.get_estimator(**kwargs).fit(X)

        if verbose: print "\tLearned model: ", NN

        return (NN,y)

    def test(self, X, model=None ):
        if model is None: model=self.model
        NN, y_train = model
        if supervised:
            y_test = NN.predict(X).flatten()
        else:
            _, ind = NN.kneighbors(X)
            y_test = y_train[ind[:,0]]
        return y_test

    def save_model(self, count, model=None):
        if model is None: model=self.model

        modelname = '{}_{}.pkl'.format(self.name,count)
        joblib.dump(model,modelname)
        return

    def load_model(self,count):
        modelname = '{}_{}.pkl'.format(self.name, count)
        model = joblib.load(modelname)
        return model

    def get_best_hyperparams(self):
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        from sklearn.model_selection import GridSearchCV
        import time

        if not supervised: return

        # 1. set hyperparamter search
        knn = GridSearchCV( self.get_estimator(),
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
        knn.fit(X, y)
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # 3. print results
        print "knn:", knn
        #print "n knn.best_estimator_.support_", len(knn.best_estimator_.support_)
        print "best_params_:", knn.best_params_
        print "best_score_:", knn.best_score_
        # print "cv_results_", knn.cv_results_ # too verbose

        # 4. update model.hparams to best ones. 
        self.hparam.update(knn.best_params_)

        return


class NearestNeighborWithPCAExperiment(NearestNeighborExperiment):

    def __init__(self, n_principal_components, **kwargs):
        super(NearestNeighborWithPCAExperiment, self).__init__(**kwargs)
        self.n_principal_components = n_principal_components

    def train(self, X=None, y=None, **kwargs):
        if X is None: X=self.X
        if y is None: y=self.y
        # Dimensionality reduction
        F = PCA(self.n_principal_components)
        F.fit(X)
        X_F = F.transform(X)

        print "\tLeast explained variance:", F.explained_variance_[-1]
        print "\tDimensionality reduction: ", X_F.shape

        # Nearest neighbor
        (NN, _) = super(NearestNeighborWithPCAExperiment, self).train(X_F, y, **kwargs)
        self.F = F

        return (NN, y)

    def test(self, X, model=None):
        if model is None: model=self.model
        NN, y_train = model
        X_F = self.F.transform(X)
        return super(NearestNeighborWithPCAExperiment, self).test(X_F, (NN, y_train))

    def save_model(self, count, model=None):
        if model is None: model=self.model

        modelname = '{}_pca_{}.pkl'.format(self.name,count)
        joblib.dump( (model, self.F) ,modelname)
        return

    def load_model(self,count):
        modelname = '{}_pca_{}.pkl'.format(self.name, count)
        (model, self.F) = joblib.load(modelname)
        return model
