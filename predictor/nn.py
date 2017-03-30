
import sys
import numpy as np
import cPickle as pickle
import os.path
import keras as nn
from sklearn.model_selection import train_test_split
from experiment_interface import Experiment

class NeuralNetworkExperiment(Experiment):

    def __init__(self, **kwargs):
        super(NeuralNetworkExperiment,self).__init__(**kwargs)
        return

    def train(self, X=None, y=None, verbose=True, **kwargs):
        if X is None: X = self.X
        if y is None: y = self.y

        print "X.shape:", X.shape
        input_dim = X.shape[1]
        if verbose:
            print "an example input vector:", X[3]

        model = get_model(input_dim = input_dim, **kwargs)

        model.compile(loss='mse', optimizer='adam')

        model.fit(X, y,
                  nb_epoch=100,
                  batch_size=256,
                  callbacks=[nn.callbacks.EarlyStopping(monitor='val_loss',verbose=1,patience=10)],
                  validation_split=.1
                  )
        print ('scoring model...')
        train_error = model.evaluate(X, y)
        #test_error = model.evaluate(X_test, y_test)
        print('\n')
        print('training error: '+str(train_error))
        #print('test error: '+str(test_error))

        return model

    def test(self, X, model=None):
        if model is None: model=self.model
        y_test = model.predict(X).flatten()
        return y_test

    def load_model(self, count):
        modelname = '{}_{}.hdf5'.format(self.name, count)
        model = nn.models.load_model(modelname)
        self.R = model.R

        # other way:
        #model = get_model()
        #model.compile(loss='mse', optimizer='adam')
        #model.load_weights(modelname)
        return model

    def save_model(self, count, model=None):
        if model is None: model=self.model
        model.R = self.R
        modelname = '{}_{}.hdf5'.format(self.name, count)
        print('saving model')
        model.save(modelname)
        return

    def get_best_hyperparams(self):
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        from sklearn.model_selection import GridSearchCV
        import time

        # 1. set hyperparamter search
        grid_nn = GridSearchCV( self.get_estimator(),
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
        grid_nn.fit(X,y)
        time_to_fit = time.time() - stime
        print "\tTime to fit: ", time_to_fit, ' s'

        # 3. print results
        print "nn_grid:", nn_grid
        print "n nn_grid.best_estimator_.support_", len(nn_grid.best_estimator_.support_)
        print "best_params_:", nn_grid.best_params_
        print "best_score_:", nn_grid.best_score_
        self.R = nn_grid.best_score_
        #print "cv_results_", nn_grid.cv_results_ # too verbose

        # 4. update model.hparams to best ones. 
        self.hparam.update(nn_grid.best_params_)

        return

    def grid_search(self):
        ''' perform a grid search on:
            - dropout_rate
            - weight_constraint
            - init_mode
            - number of neurons ( 'neurons' )
            - activation function
            - if optimizer SGD: SGD_learn_rate and SGD_momentum SGD(lr=learn_rate, momentum=momentum)
                with: learn_rate = [ 0.001, 0.01, 0.1, 0.2, 0.3 ]
                with: momentum  = [ 0.0, 0.2, 0.4, 0.6, 0.8, 0.9 ]
            '''
        #dropout_rate_grand = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ]
        #dropout_rate_small = [0.0, 0.15, 0.3 ]
        #weight_constraint = [ 1, 2, 3, 4, 5 ]
        #init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
        #activation= ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear']
        #optimizers= ['rmsprop', 'adam' ]

        kfold = KFold(n_splits=3, random_state=seed )

        estimator = KerasRegressor( build_fn=get_model, nb_epoch=10, batch_size=30, verbose=1)

        #param_grid = { 'nb_epoch' : [10,20], 'batch_size':[10, 20 ] }
        #param_grid = dict(dropout_rate=dropout_rate, weight_constraint=weight_constraint, optimizer=optimizers)
        #param_grid = dict(optimizer=optimizers)

        #param_grid = dict(
        #                 #  nb_epoch         = [10, 20 ],
        #                 #  batch_size       = [10, 20 ],
        #                 #  dropout_rate     = dropout_rate_small,
        #                 #  weight_constraint= weight_constraint  )
        #                    activation       = activation   )
        #                 #  optimizer        = optimizers   )
        #                 #  init_mode        = init_mode   )

        nn_grid = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=kfold)
        print "before fit"
        grid_result = grid.fit(X,y)

        # summarize results:
        print("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_))
        means = grid_result.cv_results_['mean_test_score']
        stds = grid_result.cv_results_['std_test_score']
        params = grid_result.cv_results_['params']
        for mean, stdev, param in zip(means, stds, params):
            print("%f (%f) with: %r" % (mean, stdev, param))

        print "grid_result:", grid_result
        return grid_result

def get_model(input_dim=9316, **kwargs):
    print ('training model...')
    dropout = .01
    model = nn.models.Sequential()
    model.add(nn.layers.Dense(128,input_dim=input_dim))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation("relu"))
    model.add(nn.layers.Dropout(dropout))
    model.add(nn.layers.Dense(16))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation("relu"))
    model.add(nn.layers.Dropout(dropout))
    model.add(nn.layers.Dense(output_dim=1))

    return model



