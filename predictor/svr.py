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

    def train(self, X=None, y=None, verbose=True,
              C = 1e4,
              gamma=1e-6,
              **kwargs):
        """ train the KNN with parameters:
            - n_neighbors: 1
        """
        if X is None: X=self.X
        if y is None: y=self.y

        svr_rbf = SVR(kernel='rbf',
                      C=C,
                      gamma=gamma,
                      cache_size=200, # number of megabytes
                      tol = 0.001
                      )
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


class SupportVectorWithPCAExperiment(SupportVectorExperiment):

    def __init__(self, n_principal_components, **kwargs):
        super(SupportVectorWithPCAExperiment, self).__init__(**kwargs)
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
        svr = super(SupportVectorWithPCAExperiment, self).train(X_F, y, **kwargs)
        self.F = F

        return svr

    def test(self, X, model=None):
        if model is None: model=self.model
        svr = model
        X_F = self.F.transform(X)
        return super(SupportVectorWithPCAExperiment, self).test(X_F, svr)

    def save_model(self, count, model=None):
        if model is None: model=self.model

        modelname = '{}_pca_{}.pkl'.format(self.name,count)
        joblib.dump( (model, self.F) ,modelname)
        return

    def load_model(self,count):
        modelname = '{}_pca_{}.pkl'.format(self.name, count)
        (model, self.F) = joblib.load(modelname)
        return model



if __name__=="__main__":
    import pickle
    print "SVR:"
    class Run:
        pass
    retrain = True
    n_folds = 3
    run = Run()
    run.identify = 'ada_'
    table = pickle.load(open('tablebin','rb'))
    regressor = SupportVectorExperiment(table=table, retrain=True, n_folds=n_folds, run=run, identify='ada_',
                                        descriptor= '1DL' )

    print "regressor:", regressor



