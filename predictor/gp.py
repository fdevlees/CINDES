import numpy as np
import pickle

from GPy.kern.src.rbf import RBF
from GPy.models.gp_regression import GPRegression
from sklearn.decomposition.pca import PCA

from experiment_interface import Experiment

from CINDES4.utils.statistics import print_stats

class GaussianProcessExperiment(Experiment):

    def __init__(self, white_noise=1e-1, max_iters=500, **kwargs):
        """
        white_noise is a constant parameter in the GP to model observational noise. [Can be used as a regularizer] 
        """
        super(GaussianProcessExperiment, self).__init__(**kwargs)
        self.white_noise = white_noise
        self.max_iters = max_iters

    def train(self, X=None, y=None, verbose=True, **kwargs):

        if X is None: X=self.X
        if y is None: y=self.y
         
        if verbose:
            print "an example input vector:", X[3]

        # Get subset for hyperparameter optimization
        ind = np.arange(X.shape[0])
        np.random.shuffle(ind)
        if X.shape[0] > 400:
            n_opt = 300
        else:
            n_opt = int( 0.75 * X.shape[0] )
        X_hyp = X[ind[:n_opt],:]
        y_hyp = y[ind[:n_opt]]
        X_test= X[ind[n_opt:],:]
        y_test= y[ind[n_opt:]]

        # Define kernel
        kernel = RBF(X.shape[1], ARD=True)

        # Optimize regularized GP hyperparameters
        gpr = GPRegression(X_hyp, y_hyp[:,None], kernel=kernel)
        gpr.Gaussian_noise.variance.constrain_fixed(self.white_noise) # White noise as regularizer
        gpr.optimize('scg', max_iters=self.max_iters)

        # get an R**2 value:
        if True:
            y_pred = gpr.predict(X_test)[0].flatten()
            stats =  print_stats(y_test, y_pred)
            self.R = stats[0]

        # Fit GP through entire training data set
        gpr.set_XY(X, y[:,None])

        if verbose: print "\tLearned model: ", gpr

        return gpr

    def test(self, X, model=None, return_std=False, **kwargs):
        if model is None: model = self.model
        if return_std:
            return model.predict(X)
        else:
            return model.predict(X)[0].flatten()

    def save_model(self, count, model=None):
        if model is None: model = self.model

        # try to attach R to the model
        model.R = self.R

        # other option:
        #modelname = 'gp_{}.npy'.format(1)
        # model = gpr = GPy.models.GPRegression
        #np.save(modelname, model.param_array)
        modelname2 = '{}_2_{}.npz'.format(self.name, count)
        with open(modelname2,'wb') as f:
            pickle.dump(model,f)

        #modelname = '{}_{}.npz'.format(self.name, count)
        #np.savez(modelname, X=self.X, y=self.y, param_array=model.param_array )
        return

    def load_model(self, count):
        #modelname = '{}_{}.npz'.format(self.name,count)
        modelname2 = '{}_2_{}.npz'.format(self.name,count)

        # other option
        # m = GPy.models(GPRegression(X,Y, initialize=False)
        #model = GPRegression(self.X, self.y, initialize=False)
        #model[:] = np.load(modelname)

        #import os.path
        #print "exist:", os.path.exists(modelname)

        #npzfile = np.load(modelname)
        #y = npzfile['y'][:,None]
        #model = GPRegression( npzfile['X'], y, initialize=False)
        #model.update_model(False)
        #model.initialize_parameter()
        #array = tuple(npzfile['param_array'])
        #print "array:", array
        #model[:] = array
        #model.update_model(True)

        model = pickle.load(open(modelname2,'rb'))
        print "loaded model:", model

        # try to take R from the model:
        try:
            self.R = model.R
        except AttributeError:
            self.R = 0.0
        return model


class GaussianProcessWithPCAExperiment(GaussianProcessExperiment):

    def __init__(self, white_noise=1e-1, n_principal_components=50, **kwargs):
        super(GaussianProcessWithPCAExperiment, self).__init__(white_noise, **kwargs)
        self.n_principal_components = n_principal_components

    def train(self, X=None, y=None, **kwargs):
        if X is None: X=self.X
        if y is None: y=self.y
        # Dimensionality reduction
        F = PCA(self.n_principal_components)
        F.fit(X)
        X_F = F.transform(X)
        self.F = F

        print "\tLeast explained variance:", F.explained_variance_[-1]
        print "\tDimensionality reduction: ", X_F.shape

        gp = super(GaussianProcessWithPCAExperiment, self).train(X_F, y)

        return gp

    def test(self, X, model=None, **kwargs):
        if model is None: model = self.model
        gp = model
        X_F = self.F.transform(X)
        return super(GaussianProcessWithPCAExperiment, self).test(X_F, gp, **kwargs)

    def save_model(self, count, model=None):
        if model is None: model = self.model
        model.R = self.R

        modelname2 = '{}_pca_{}.npz'.format(self.name, count)
        with open(modelname2,'wb') as f:
            pickle.dump((model,self.F),f)
        return

    def load_model(self, count):
        modelname2 = '{}_pca_{}.npz'.format(self.name,count)

        model, self.F  = pickle.load(open(modelname2,'rb'))

        print "loaded model:", model
        try:
            self.R = model.R
        except AttributeError:
            self.R = 0.0
        return model




class GaussianProcessExperiment_skl(Experiment):

    def __init__(self, white_noise=1e-1, **kwargs):
        """
        white_noise is a constant parameter in the GP to model observational noise. [Can be used as a regularizer] 
        """
        super(GaussianProcessExperiment_skl, self).__init__(**kwargs)
        self.white_noise = white_noise

    def train(self, X=None, y=None, verbose=True, **kwargs):

        if X is None: X=self.X
        if y is None: y=self.y

        # Get subset for hyperparameter optimization
        ind = np.arange(X.shape[0])
        np.random.shuffle(ind)
        n_opt = min((300, X.shape[0]))
        X_hyp = X[ind[:n_opt],:]
        y_hyp = y[ind[:n_opt]]

        # Define kernel # RBF ( alias Squared Exponential (SE) )
        from sklearn.gaussian_process.kernels import RBF, WhiteKernel
        #kernel = RBF( length_scale=1.0, length_scale_bounds=(1.e-5, 100000.0)
        kernel = 1.0 * RBF( length_scale=1.0, length_scale_bounds=(1.e-5, 100000.0) )\
                + WhiteKernel( noise_level=1.e-5, noise_level_bounds=(1e-10, 1e+1 ) )
        #kernel = RBF(X.shape[1], ARD=True)

        # set estimator
        from sklearn.gaussian_process import GaussianProcessRegressor
        gpr = GaussianProcessRegressor(kernel=kernel, alpha= 0.0, n_restarts_optimizer=9)

        # Optimize regularized GP hyperparameters
        #gpr = GPRegression(X_hyp, y_hyp[:,None], kernel=kernel)
        #gpr.Gaussian_noise.variance.constrain_fixed(self.white_noise) # White noise as regularizer
        #gpr.optimize('scg', max_iters=500)

        # Fit GP through entire training data set
        #gpr.set_XY(X, y[:,None])
        gpr.fit(X,y)

        if verbose: print "\tLearned model: ", gpr
        print "score:", gpr.score(X,y)

        return gpr

    def test(self, X, model=None, **kwargs):
        if model is None: model = self.model
        #return model.predict(X)[0].flatten()
        return model.predict(X, return_std=False )

    def save_model(self, count, model=None):
        if model is None: model = self.model
        model.R = self.R
        
        # other option:
        #modelname = 'gp_{}.npy'.format(1)
        # model = gpr = GPy.models.GPRegression
        #np.save(modelname, model.param_array)
        modelname2 = '{}_2_{}.npz'.format(self.name, count)
        with open(modelname2,'wb') as f:
            pickle.dump(model,f)

        #modelname = '{}_{}.npz'.format(self.name, count)
        #np.savez(modelname, X=self.X, y=self.y, param_array=model.param_array )
        return

    def load_model(self, count):
        #modelname = '{}_{}.npz'.format(self.name,count)
        modelname2 = '{}_2_{}.npz'.format(self.name,count)

        # other option
        # m = GPy.models(GPRegression(X,Y, initialize=False)
        #model = GPRegression(self.X, self.y, initialize=False)
        #model[:] = np.load(modelname)

        #import os.path
        #print "exist:", os.path.exists(modelname)

        #npzfile = np.load(modelname)
        #y = npzfile['y'][:,None]
        #model = GPRegression( npzfile['X'], y, initialize=False)
        #model.update_model(False)
        #model.initialize_parameter()
        #array = tuple(npzfile['param_array'])
        #print "array:", array
        #model[:] = array
        #model.update_model(True)

        model = pickle.load(open(modelname2,'rb'))
        print "loaded model:", model
        try:
            self.R = model.R
        except AttributeError:
            self.R = 0.0
        return model
