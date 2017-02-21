from sklearn.decomposition.pca import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.externals import joblib

from experiment_interface import Experiment


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

    def train(self, X=None, y=None, verbose=False, **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        NN = NearestNeighbors(n_neighbors=1).fit(X)

        if verbose: print "\tLearned model: ", NN

        return (NN,y)

    def test(self, X, model=None ):
        if model is None: model=self.model
        NN, y_train = model
        _, ind = NN.kneighbors(X)
        return y_train[ind[:,0]]

    def save_model(self, model=None):
        if model is None: model=self.model

        modelname = 'knn_{}.pkl'.format(1)
        joblib.dump(model,modelname)
        return

    def load_model(self):
        modelname = 'knn_{}.pkl'.format(1)
        model = joblib.load(modelname)
        return model


class NearestNeighborWithPCAExperiment(NearestNeighborExperiment):

    def __init__(self, setting, n_folds, n_principal_components):
        super(NearestNeighborWithPCAExperiment, self).__init__(setting, n_folds)
        self.n_principal_components = n_principal_components

    def train(self, X, y):
        # Dimensionality reduction
        F = PCA(self.n_principal_components)
        F.fit(X)
        X_F = F.transform(X)

        print "\tLeast explained variance:", F.explained_variance_[-1]
        print "\tDimensionality reduction: ", X_F.shape

        # Nearest neighbor
        (NN, _), log = super(NearestNeighborWithPCAExperiment, self).train(X_F, y)
        log["pca"] = F

        return (F, NN, y), log

    def test(self, X, model):
        F, NN, y_train = model
        X_F = F.transform(X)
        return super(NearestNeighborWithPCAExperiment, self).test(X_F, (NN, y_train))
