from sklearn.decomposition.pca import PCA
from sklearn.externals import joblib
from sklearn.utils import resample
from sklearn.base import BaseEstimator
import numpy as np
from qml.kernels import gaussian_kernel
from qml.math import cho_solve

from experiment_interface import Experiment

# krr with rbf kernel. C controls simplisity or decision surface. High C will
# try to fit all data and select more support vector. Low C will give a more
# smooth surface. gamma parameter defines how far the influence of a single
# training example reaches, with low values meaning far and high values
# meaning close. low gamma value high bias, high values high variance.

class KernelRidgeExperiment(Experiment, BaseEstimator):

    def __init__(self, sigma=1000, labda=1e-8, **kwargs):
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
        self.sigma = sigma
        self.labda = labda
        #self.kernel= 'arad'
        self.kernel= 'gaussian'

        # set default hyperparam ( super sets self.hparams to dict() ) so beware of order.
        self.hparam = { 'labda':1e-8,
                        'sigma':1000 }  # gamma parameter is specific for rbf/laplacian kernel

        self.hparam_grid = {'sigma': np.logspace(-10,5,10),
                            'labda': np.logspace(-10,5,10) }
                           # 'kernel': ['rbf','laplacian'] }
        # see if new defaults are given via input
        #for key in self.hparam:
        #    if key in kwargs:
        #        self.hparam[key] = kwargs[key]
        #        print "new default hyperparameter:", key, kwargs[key]

        return

    def train(self, X=None, y=None, verbose=True, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y
        if self.kernel=='gaussian':
            K = gaussian_kernel(X, X, self.sigma)
            K[np.diag_indices_from(K)] += self.labda
        elif self.kernel=='arad':
            #from qml.arad_kernels import get_atomic_kernels_arad
            from qml.arad import get_local_kernels_arad
            sigmas=[50.0,100.0,200.0]
            K = get_local_kernels_arad(X,X,sigmas)
            K[np.diag_indices_from(K)] += self.labda
        self.alpha = cho_solve(K, y)
        model = (X,self.alpha)

        # extra;
        Yss = np.dot(K, self.alpha)
        mae_training_error = np.mean(np.abs(Yss - y))
        rmse_training_error = np.sqrt(np.mean(np.square(Yss - y)))
        #if verbose: print "\tLearned model: ", krr_rbf

        return model

    def test(self, X, model=None, **kwargs ):
        if model is None: model=self.model
        X_training, alpha = model
        if self.kernel=='gaussian':
            Ks = gaussian_kernel(X, X_training, self.sigma)
        elif self.kernel=='arad':
            #from qml.arad_kernels import get_atomic_kernels_arad
            from qml.arad import get_local_kernels_arad
            sigmas=[1.0,100.0,10000.0]
            Ks = get_local_kernels_arad(X,X_training ,sigmas)
        #Ks = gaussian_kernel(X, X_training, self.sigma)
        Y_predicted = np.dot(Ks, alpha)
        return Y_predicted

    def load_model(self, count, path='.', load_kernel=False, **kwargs):
        filename = "{}/{}_{}_".format(path, self.name, count)
        self.X = np.load(filename + "X.npy")
        self.y = np.load(filename + "Y.npy")
        self.alpha = np.load(filename + "alpha.npy")
        #if load_kernel:
        #    self.K = np.load(path + "/K.npy")
        return


    def save_model(self, count, path='.', save_kernel=False, **kwargs):
        filename = "{}/{}_{}_".format(path, self.name, count)
        np.save(filename + "X.npy", self.X)
        np.save(filename + "Y.npy", self.y)
        np.save(filename + "alpha.npy", self.alpha)
        #if save_kernel:
        #    np.save(path + "/K.npy")
        return

#    def save_model(self, count, model=None):
#        if model is None: model=self.model
#        model.R = self.R
#
#        modelname = '{}_{}.pkl'.format(self.name, count)
#        joblib.dump(model,modelname)
#        return
#
#    def load_model(self,count):
#        modelname = '{}_{}.pkl'.format(self.name, count)
#        model = joblib.load(modelname)
#        self.R = model.R
#        return model

    def fit(self, X, y=None):
        K = gaussian_kernel(X, X, self.sigma)
        K[np.diag_indices_from(K)] += self.labda
        alpha = cho_solve(K, y)
        self.model=(X,alpha)
        return self

    def predict(self, X, y=None, **kwargs):
        X_training, alpha = self.model
        Ks = gaussian_kernel(X, X_training, self.sigma)
        return np.dot(Ks, alpha)

    def score(self, X, y=None):
        from scipy.stats import pearsonr
        self.R = pearsonr( self.predict(X), y)[0]
        return self.R

    def get_estimator(self):
        return self


