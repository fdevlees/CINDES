from sklearn.decomposition.pca import PCA
from sklearn.neighbors import NearestNeighbors

from experiment_interface import Experiment
        
class NearestNeighborExperiment(Experiment):
    
    def __init__(self, setting, n_folds):
        super(NearestNeighborExperiment, self).__init__(setting, n_folds)

    def train(self, X, y):
        NN = NearestNeighbors(n_neighbors=1).fit(X)

        print "\tLearned model: ", NN
        
        return (NN, y), {"nn": NN, "y_train": y}
    
    def test(self, X, model):
        NN, y_train = model
        _, ind = NN.kneighbors(X)
        return y_train[ind[:,0]]
    
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
