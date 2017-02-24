from sklearn.decomposition.pca import PCA
from sklearn.externals import joblib
from sklearn.svm import SVR

from experiment_interface import Experiment

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

    def train(self, X=None, y=None, verbose=False, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        svr_rbf = SVR(kernel='rbf', C=1e4, gamma=1e-6)
        print "Fitting...",
        svr_rbf.fit(X, y)

        if verbose: print "\tLearned model: ", svr_rbf

        return svr_rbf

    def test(self, X, model=None ):
        if model is None: model=self.model
        return model.predict(X).flatten()

    def save_model(self, count, model=None):
        if model is None: model=self.model

        modelname = '{}_{}.pkl'.format(self.name, count)
        joblib.dump(model,modelname)
        return

    def load_model(self,count):
        modelname = '{}_{}.pkl'.format(self.name, count)
        model = joblib.load(modelname)
        return model


