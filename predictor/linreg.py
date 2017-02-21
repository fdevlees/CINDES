from sklearn.decomposition.pca import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.externals import joblib

from experiment_interface import Experiment
from descriptor import get_X_1D

class LinRegOneExperiment(Experiment):

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
        super(LinRegOneExperiment, self).__init__(**kwargs)

    def train(self, X=None, y=None, verbose=False, subtype='ridge', intercept=False, twosite=False, alpha=0,  **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        if subtype == 'LinearRegression':
            clf = linear_model.LinearRegression(fit_intercept=intercept)
        elif subtype in ['Ridge']:
            clf = linear_model.Ridge(alpha=alpha,fit_intercept=intercept,tol=0.001,solver='auto')
        elif subtype in ['RidgeCV','ridgecv']:
            clf = linear_model.RidgeCV(alphas=alpha, fit_intercept=intercept, store_cv_values=True)
        elif subtype in ['Lasso']:
            clf = linear_model.Lasso(alpha=alpha,fit_intercept=intercept,tol=0.001)
        elif subtype in ['ElasticNet']:
            l1_ratio = 0.1 #default 0.5
            alpha = 1e-2
            clf = linear_model.ElasticNet(alpha=alpha,l1_ratio=l1_ratio, fit_intercept=intercept,tol=0.001)

        clf.fit(X,y)

        if verbose: print "\tLearned model: ", clf

        return clf

    def test(self, X, model=None ):
        if model is None: model=self.model
        y_pred = model.predict(X)
        return y_pred.flatten()

    def predict(self, molecules):
        ''' get molecules list '''
        indices = [ mol.index for mol in molecules ]
        X_pred = get_X_1D(indices,self.descriptor, self.run.identify)
        y_pred = self.test( X_pred)
        print "y_pred:", y_pred
        return y_pred

    def save_model(self, model=None):
        if model is None: model=self.model
        modelname = 'linreg1_{}.pkl'.format(1)
        joblib.dump(model,modelname)
        return

    def load_model(self):
        modelname = 'linreg1_{}.pkl'.format(1)
        model = joblib.load(modelname)
        return model

