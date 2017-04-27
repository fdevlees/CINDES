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
                        "weights": "uniform",
                        "supervised": True}

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
        if self.hparam['supervised']:
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
        
        if self.hparam['supervised']:
            NN = self.get_estimator(**kwargs).fit(X,y)
        else:
            NN = self.get_estimator(**kwargs).fit(X)

        if verbose: print "\tLearned model: ", NN

        return (NN,y)

    def test(self, X, model=None, **kwargs ):
        if model is None: model=self.model
        NN, y_train = model
        if self.hparam['supervised']:
            y_test = NN.predict(X).flatten()
        else:
            _, ind = NN.kneighbors(X)
            y_test = y_train[ind[:,0]]
        return y_test

    def save_model(self, count, model=None):
        if model is None: model=self.model
        R = self.R

        modelname = '{}_{}.pkl'.format(self.name,count)
        joblib.dump((model,R),modelname)
        return

    def load_model(self,count):
        modelname = '{}_{}.pkl'.format(self.name, count)
        model, self.R = joblib.load(modelname)
        return model

    def get_best_hyperparams(self, **kwargs):
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        if not self.hparam['supervised']:
            self.cross_val(plot=True)
        else:
            cls_gs = super(NearestNeighborExperiment, self).get_best_hyperparams(**kwargs)
        return cls_gs


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

    def test(self, X, model=None, **kwargs):
        if model is None: model=self.model
        NN, y_train = model
        X_F = self.F.transform(X)
        return super(NearestNeighborWithPCAExperiment, self).test(X_F, (NN, y_train))

    def save_model(self, count, model=None):
        if model is None: model=self.model
        R = self.R

        modelname = '{}_pca_{}.pkl'.format(self.name,count)
        joblib.dump( (model, self.F, R) ,modelname)
        return

    def load_model(self,count):
        modelname = '{}_pca_{}.pkl'.format(self.name, count)
        (model, self.F, self.R) = joblib.load(modelname)
        return model
