import numpy as np

from GPy.kern.src.rbf import RBF
from GPy.models.gp_regression import GPRegression
from sklearn.decomposition.pca import PCA

from experiment_interface import Experiment

class GaussianProcessExperiment(Experiment):

    def __init__(self, white_noise=1e-1, **kwargs):
        """
        white_noise is a constant parameter in the GP to model observational noise. [Can be used as a regularizer] 
        """
        super(GaussianProcessExperiment, self).__init__(**kwargs)
        self.white_noise = white_noise

    def train(self, X=None, y=None, verbose=False, **kwargs):

        if X is None: X=self.X
        if y is None: y=self.y

        # Get subset for hyperparameter optimization
        ind = np.arange(X.shape[0])
        np.random.shuffle(ind)
        n_opt = min((300, X.shape[0]))
        X_hyp = X[ind[:n_opt],:]
        y_hyp = y[ind[:n_opt]]

        # Define kernel
        kernel = RBF(X.shape[1], ARD=True)

        # Optimize regularized GP hyperparameters
        gpr = GPRegression(X_hyp, y_hyp[:,None], kernel=kernel)
        gpr.Gaussian_noise.variance.constrain_fixed(self.white_noise) # White noise as regularizer
        gpr.optimize('scg', max_iters=500)

        # Fit GP through entire training data set
        gpr.set_XY(X, y[:,None])

        if verbose: print "\tLearned model: ", gpr

        return gpr

    def test(self, X, model=None):
        if model is None: model = self.model
        return model.predict(X)[0].flatten()

    def save_model(self, count, model=None):
        if model is None: model = self.model

        # other option:
        #modelname = 'gp_{}.npy'.format(1)
        # model = gpr = GPy.models.GPRegression
        #np.save(modelname, model.param_array)

        modelname = '{}_{}.npz'.format(self.name, count)
        np.savez(modelname, X=self.X, y=self.y, param_array=model.param_array )
        return

    def load_model(self, count):
        modelname = '{}_{}.npz'.format(self.name,count)

        # other option
        # m = GPy.models(GPRegression(X,Y, initialize=False)
        #model = GPRegression(self.X, self.y, initialize=False)
        #model[:] = np.load(modelname)

        import os.path
        print "exist:", os.path.exists(modelname)

        npzfile = np.load(modelname)
        model = GPRegression( npzfile['X'], npzfile['y'], initialize=False)
        model[:] = npzfile('param_array')
        print "loaded model:", model
        return model


class GaussianProcessWithPCAExperiment(GaussianProcessExperiment):

    def __init__(self, setting, n_folds, white_noise, n_principal_components):
        super(GaussianProcessWithPCAExperiment, self).__init__(setting, n_folds, white_noise)
        self.n_principal_components = n_principal_components

    def train(self, X, y):
        # Dimensionality reduction
        F = PCA(self.n_principal_components)
        F.fit(X)
        X_F = F.transform(X)

        print "\tLeast explained variance:", F.explained_variance_[-1]
        print "\tDimensionality reduction: ", X_F.shape

        gp, log = super(GaussianProcessWithPCAExperiment, self).train(X_F, y)
        log["pca"] = F

        return (F, gp), log

    def test(self, X, model):
        F, gp = model
        X_F = F.transform(X)
        return super(GaussianProcessWithPCAExperiment, self).test(X_F, gp)
