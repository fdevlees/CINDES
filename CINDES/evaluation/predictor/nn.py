
import sys
import numpy as np
import pickle as pickle
import os.path
import keras as nn
from sklearn.model_selection import train_test_split
from .experiment_interface import Experiment
from keras.wrappers.scikit_learn import KerasRegressor

input_dim = None

class NeuralNetworkExperiment(Experiment):

    def __init__(self, **kwargs):
        super(NeuralNetworkExperiment,self).__init__(**kwargs)
        self.hparam = dict(
                           nb_epoch         = 100,
                           batch_size       = 256,
                           dropout_rate     = 0.3,
                           weight_constraint= 3,
                           activation       = 'relu',
                           optimizer        = 'sgd',
                           init_mode        = 'uniform',
                           loss             = 'mean_squared_error'
                           )
        self.hparam_grid = get_grid()
        return

    def train(self, X=None, y=None, verbose=True, **kwargs):
        if X is None: X = self.X
        if y is None: y = self.y

        print("X[0].shape:", X[0].shape)
        print("X.shape:", X.shape)
        input_dim = X[0].shape[0]
        if verbose:
            print("an example input vector:", X[3])

        print("self.kwargs:", kwargs)
        print("self.hparam:", self.hparam)
        model = set_model(input_dim = input_dim, **self.hparam)

        #model.compile(loss='mse', optimizer='adam')

        model.fit(X, y,
                  nb_epoch=self.hparam['nb_epoch'],
                  batch_size= self.hparam['batch_size'],
                  #callbacks=[nn.callbacks.EarlyStopping(monitor='val_loss',verbose=1,patience=10)],
                  #validation_split=.1
                  )
        print ('scoring model...')
        train_error = model.evaluate(X, y)
        #test_error = model.evaluate(X_test, y_test)
        print('\n')
        print(('training error: '+str(train_error)))
        #print('test error: '+str(test_error))

        return model

    def test(self, X, model=None):
        if model is None: model=self.model
        y_test = model.predict(X).flatten()
        return y_test

    def load_model(self, count):
        modelname = '{}_{}.hdf5'.format(self.name, count)
        model = nn.models.load_model(modelname)
        try:
            self.R = model.R
        except AttributeError:
            print("model has no R value!")
            self.R = None

        # other way:
        #model = set_model()
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

    def get_best_hyperparams_old(self):
        ''' hyperparameter search with use of the sklearn GridSearchCV function '''
        from sklearn.model_selection import GridSearchCV
        import time

        # 1. set hyperparamter search
        grid_nn = GridSearchCV( set_model(),
                            cv= self.n_folds,
                            n_jobs=8,
                            param_grid = self.hparam_grid )

        # 2. do search on dataset
        if True:
            n_train = 300
            #X = self.X[:n_train]
            #y = self.y[:n_train]
            X, y = resample(self.X, self.y, n_samples=n_train)
            print("restricted hparamopt to only {} samples".format(n_train))
        else:
            X = self.X
            y = self.y
        stime = time.time()
        grid_nn.fit(X,y)
        time_to_fit = time.time() - stime
        print("\tTime to fit: ", time_to_fit, ' s')

        # 3. print results
        print("nn_grid:", nn_grid)
        print("n nn_grid.best_estimator_.support_", len(nn_grid.best_estimator_.support_))
        print("best_params_:", nn_grid.best_params_)
        print("best_score_:", nn_grid.best_score_)
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

        kfold = KFold(n_splits=3, random_state=seed )

        estimator = KerasRegressor( build_fn=set_model, nb_epoch=100, batch_size=30, verbose=0)

        #param_grid = { 'nb_epoch' : [10,20], 'batch_size':[10, 20 ] }
        #param_grid = dict(dropout_rate=dropout_rate, weight_constraint=weight_constraint, optimizer=optimizers)
        #param_grid = dict(optimizer=optimizers)


        nn_grid = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=kfold)
        print("before fit")
        grid_result = grid.fit(X,y)

        # summarize results:
        print(("Best: %f using %s" % (grid_result.best_score_, grid_result.best_params_)))
        means = grid_result.cv_results_['mean_test_score']
        stds = grid_result.cv_results_['std_test_score']
        params = grid_result.cv_results_['params']
        for mean, stdev, param in zip(means, stds, params):
            print(("%f (%f) with: %r" % (mean, stdev, param)))

        print("grid_result:", grid_result)
        return grid_result

    def get_estimator(self):
        from functools import partial, update_wrapper
        print(self.X.shape)
        print(self.X[0].shape)
        global input_dim
        input_dim = self.X[0].shape

        #model = partial( set_model, input_dim=self.X[0].shape )
        #model = update_wrapper( model, set_model)
        #model = my_partial( set_model, self.X[0].shape)

        model = my_partial( self.X[0].shape[0] )

        estimator = KerasRegressor(
                                   build_fn= model,
                                   nb_epoch=100,
                                   batch_size=256,
                                   verbose=0,
        #                           callbacks=[nn.callbacks.EarlyStopping(monitor='val_loss',verbose=1,patience=10)],
        #                           validation_split=.1
                                   )
        return estimator

def set_model(
              input_dim=None,
              dropout_rate=0.1,
              weight_constraint=3,
              activation='relu',
              optimizer='sgd',
              loss='binary_crossentropy',
              init_mode='uniform',
              **kwargs
              ):
    #print ('training model...')
    model = nn.models.Sequential()

    # layer 1
    model.add(nn.layers.Dense(128,input_dim=input_dim ))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation(activation))
    model.add(nn.layers.Dropout(dropout_rate))

    # layer 2
    #model.add(nn.layers.Dense(16, kernel_initializer=init_mode))
    model.add(nn.layers.Dense(50))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation(activation))
    model.add(nn.layers.Dropout(dropout_rate))

    # layer 3 output layer
    model.add(nn.layers.Dense(output_dim=1))

    # compilation
    model.compile(loss=loss, optimizer=optimizer)

    return model

def get_grid():
    dropout_rate_grand = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ]
    dropout_rate_small = [ 0.3, 0.45, 0.6, 0.8, 0.9 ]
    weight_constraint = [ 1, 5 ]
    init_mode  = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
    activation = ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear']
    optimizers = ['rmsprop', 'adam' ]
    nb_epochs  = [ 50, 200 ]
    batch_sizes= [ 50, 200 ]

    grid = dict(
                nb_epoch         = nb_epochs,
                batch_size       = batch_sizes,
                dropout_rate     = dropout_rate_small,
                weight_constraint= weight_constraint
                #activation       = activation
                #optimizer        = optimizers,
                #init_mode        = init_mode
                )
    return grid

def set_model_gridsearch(
              dropout_rate=0.1,
              weight_constraint=3,
              activation='relu',
              optimizer='adam',
              loss='mse',
              init_mode='uniform'
              ):
    #print ('training model...')
    model = nn.models.Sequential()

    # layer 1
    model.add(nn.layers.Dense(128,input_dim=input_dim ))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation(activation))
    model.add(nn.layers.Dropout(dropout_rate))

    # layer 2
    #model.add(nn.layers.Dense(16, kernel_initializer=init_mode))
    model.add(nn.layers.Dense(16))
    model.add(nn.layers.normalization.BatchNormalization())
    model.add(nn.layers.Activation(activation))
    model.add(nn.layers.Dropout(dropout_rate))

    # layer 3 output layer
    model.add(nn.layers.Dense(output_dim=1))

    # compilation
    model.compile(loss=loss, optimizer=optimizer)
    return model
 
def my_partial( input_dim_local, *args, **kwargs):
    global input_dim
    input_dim = input_dim_local
    return set_model_gridsearch
