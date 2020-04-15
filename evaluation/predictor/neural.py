debug=True

from CINDES.utils.writings import log_io, sprint

from keras.models import Sequential
from keras.layers import Dense, Activation, Dropout
from keras.optimizers import SGD
from keras.regularizers import l2, activity_l2
from keras.constraints import maxnorm
import numpy as np
import matplotlib.pyplot as plt

from keras.wrappers.scikit_learn import KerasRegressor, KerasClassifier
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

seed = np.random.seed(7)


def run_example(X, Y):
    # tweaked example from http://machinelearningmastery.com/regression-tutorial-keras-deep-learning-library-python/

    def baseline_model():
        # create model
        model = Sequential()
        model.add( Dense( 12, input_dim=3136, init='normal', activation='relu') )
        model.add( Dense(  1, init='normal' , activation = 'sigmoid') )
        # compile model
        model.compile(loss='mean_squared_error', optimizer='adam')
        return model

    kfold = KFold(n_splits=5, random_state=seed )
    if False:
        estimator = KerasRegressor( build_fn = baseline_model, nb_epoch=50, batch_size=10, verbose=1)
        results = cross_val_score(estimator, X, Y, cv=kfold)
        print(("Results: %.2f (%.2f) MSE" % (results.mean(), results.std())))
    else:
        estimators = []
        estimators.append(( 'standardize', StandardScaler() ) )
        estimators.append(( 'mlp'        , KerasRegressor( build_fn = baseline_model, nb_epoch=50, batch_size=5, verbose=1)))
        pipeline = Pipeline(estimators)
        results = cross_val_score(pipeline, X, Y, cv=kfold)
        print("Standardized: %.2f (%.2f) MSE" % (results.mean(), results.std() ))

    return


def plot_error(data1,data2, *args, **kwargs):
    plt.plot(data1,data2, *args, **kwargs)
    return

def RMSE(data1,data2):
    assert len(data1)==len(data2)
    npoints = len(data1)
    rmse = 0
    for i in range(npoints):
        rmse += ( data1[i] - data2[i] )**2
    return np.sqrt( rmse / float(npoints) )

def get_model(input_dim=3136, 
              dropout_rate=0.3, 
              weight_constraint=3, 
              optimizer='rmsprop', 
              activation='relu',
              init = 'normal',
              nn_lay1 = 128,
              nn_lay2 =  64
              ):
    # 1. input layer
    model = Sequential()
    #init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
    # init = 'normal'
    #activation = 'relu' # also try relu / sigmoid
    activation_final = 'relu' # relu / sigmoid / softmax #softmax performs badly
    #dropout_final=0.02
    w_l2 = 1e-20
    a_l2 = 1e-20
    #w_l2 = 0.
    #a_l2 = 0.
    if True: #regularization
        # option 2 # init normal / uniform
        model.add(Dense(output_dim=nn_lay1, input_dim=input_dim, init=init,
                             W_regularizer=l2(w_l2),
                             activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
        model.add(Activation(activation))
        model.add(Dropout(dropout_rate))
        # 2. Hidden layers with weights. 
        model.add(Dense(output_dim=nn_lay2, input_dim=input_dim, init=init,
                             W_regularizer=l2(w_l2),
                             W_constraint=maxnorm(weight_constraint),
                             activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
        model.add(Activation(activation))
        model.add(Dropout(dropout_rate))
        # 3. Output layer
        model.add(Dense(output_dim=1, input_dim=input_dim, init=init,
                             W_regularizer=l2(w_l2),
                             activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
        model.add(Activation(activation_final))
        model.add(Dropout(dropout_rate))
    model.compile(loss='mean_squared_error', optimizer=optimizer )
    return model


def test_model(dropout_rate=0.0, weight_constraint=0):
    # 1. input layer
    input_dim=3136
    init='uniform'
    activation = 'relu' # also try relu / sigmoid / linear
    model = Sequential()
    model.add(Dense(output_dim=128, input_dim=input_dim, init=init ))
    model.add(Activation(activation))
    model.add(Dropout(dropout_rate))
    # 2. Hidden layers with weights. 
    model.add(Dense(output_dim=64, input_dim=input_dim, init=init,
                         W_constraint=maxnorm(weight_constraint)))
    model.add(Activation(activation))
    model.add(Dropout(dropout_rate))
    # 3. Output layer
    model.add(Dense(output_dim=1, input_dim=input_dim, init=init ))
    model.add(Activation(activation))
    model.add(Dropout(dropout_rate))

    model.compile(loss='mean_squared_error', optimizer='adam' )
    return model



class neural(object):

    def __init__(self, X, y, fraction=False, scaling=True, grid_search=True, *kwargs):
        self.scaling = scaling
        if self.scaling:
            self.X = self.scaler(X)
            #newX_scaled = scaler.transform(newX)
        else:
            self.X = X
        self.y = y
        assert len(self.X)==len(self.y)
        self.ndim = len(self.X)
        self.xdim = len(self.X[0])
        print("ndim,xdim", self.ndim, self.xdim)
        return

    def scaler(self,X):
        from sklearn import preprocessing
        self.scaler = preprocessing.StandardScaler().fit(X)
        X_scaled = self.scaler.transform(X)
        return X_scaled

    def run_procedure_testing(self):
        #self.setup_works()
        if False:
            #self.gridsearch(self.X, self.y) # does not work seems to come in infinite loop or so. is not exiting with ^C
            self.CV(self.X, self.y)
        if fraction:
            from sklearn.model_selection import train_test_split
            X_train, X_test, y_train, y_test = train_test_split(self.X,self.y, train_size = fraction)
        self.compile()
        self.fit(X_train, y_train, X_test, y_test)
        scores = self.test(X_test, y_test)
        if True:
            pred_train = self.test_new( X_train)
            pred_test  = self.test_new( X_test)
            sprint(10, pred_train, y_train)
            sprint(10, pred_test, y_test)
            print('train error:', RMSE( pred_train, y_train))
            print(' test error:', RMSE( pred_test, y_test))
            plot_error( pred_train, y_train, 'ro')
            plot_error( pred_test,  y_test, 'bo', alpha=0.5)
            plt.show()
        print(scores)
        self.scores = scores
        return

    def run_procedure_INDES(self, X_new):

        # optimize parameters
        if self.grid_search:
            hyperparameters = self.grid_search(self.X, self.y)
            print("hyperparameters:", hyperparameters)
        else:
            self.model = get_model(input_dim        = self.xdim,
                                   dropout_rate     = 0.3,
                                   weight_constraint= 3,
                                   optimizer        = 'adam',
                                   activation       = 'relu',
                                   init             = 'normal',
                                   nn_lay1          = 128,
                                   nn_lay2          =  64)
            self.CV(self.X, self.y)

        # train on whole set
        self.fit(self.X, self.y)

        # predict:
        if self.scaling:
            X_new = self.scaler(X_new)
        y_new  = self.test_new( X_new )
        return y_new

    def CV(self, X, y):
        estimator = KerasRegressor(build_fn=test_model, nb_epoch=100, batch_size=5, verbose=1)
        kfold = KFold(n_splits=3, random_state=seed)
        results = cross_val_score( estimator, X,y, cv=kfold)
        print(("Results: %.2f (%.2f) MSE" % (results.mean(), results.std())))
        print("results:", results)
        return

    def grid_search(self, X, y):
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
        #GS_model = KerasClassifier( build_fn=self.model, nb_epoch=10, dropout_rate=0.1 )
        dropout_rate_grand = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7 ]
        dropout_rate_small = [0.0, 0.15, 0.3 ]
        weight_constraint = [ 1, 2, 3, 4, 5 ]
        init_mode = ['uniform', 'lecun_uniform', 'normal', 'zero', 'glorot_normal', 'glorot_uniform', 'he_normal', 'he_uniform']
        activation= ['softmax', 'softplus', 'softsign', 'relu', 'tanh', 'sigmoid', 'hard_sigmoid', 'linear']
        optimizers= ['rmsprop', 'adam' ]

        kfold = KFold(n_splits=3, random_state=seed )

        estimator = KerasRegressor( build_fn=get_model, nb_epoch=10, batch_size=30, verbose=1)

        #param_grid = { 'nb_epoch' : [10,20], 'batch_size':[10, 20 ] }
        #param_grid = dict(dropout_rate=dropout_rate, weight_constraint=weight_constraint, optimizer=optimizers)
        #param_grid = dict(optimizer=optimizers)

        param_grid = dict(
                         #  nb_epoch         = [10, 20 ],
                         #  batch_size       = [10, 20 ],
                         #  dropout_rate     = dropout_rate_small,
                         #  weight_constraint= weight_constraint  )
                            activation       = activation   )
                         #  optimizer        = optimizers   )
                         #  init_mode        = init_mode   )

        grid = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=kfold)
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

    def setup_works(self):
        # 1
        #self.model.add(Dense(output_dim=25, input_dim=self.xdim, init='uniform')) #100 works. 150 not
        #self.model.add(Activation('relu'))
        #self.model.add(Dense(output_dim=12, init='uniform')) # if output_dim 100 output preds become 0.0000
        #self.model.add(Activation('relu'))
        #self.model.add(Dense(output_dim=1 ))
        #self.model.add(Activation("sigmoid")) # this sigmoid is important. 'relu' doesn't work neither does 'softmax'
        # 2
        #self.model.add(Dense(output_dim=64, input_dim=self.xdim, init='normal')) #100 works. 150 not
        #self.model.add(Activation('relu'))
        #self.model.add(Dense(output_dim=16, init='normal')) # if output_dim 100 output preds become 0.0000
        #self.model.add(Activation('relu'))
        #self.model.add(Dense(output_dim=1 , init='normal'))
        #self.model.add(Activation("sigmoid"))
        # 3 # first one that worked on the RUN42 dataset
        init = 'normal'
        activation = 'sigmoid' # also try ReLu
        w_l2 = 1e-14
        a_l2 = 1e-14
        if True: #regularization
            from keras.regularizers import l2, activity_l2
            self.model.add(Dense(output_dim=32, input_dim=self.xdim, init=init,
                                 W_regularizer=l2(w_l2),
                                 activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
            self.model.add(Activation(activation))
            self.model.add(Dropout(self.dropout))
            self.model.add(Dense(output_dim=16, input_dim=self.xdim, init=init,
                                 W_regularizer=l2(w_l2),
                                 activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
            self.model.add(Activation(activation))
            self.model.add(Dropout(self.dropout))
            self.model.add(Dense(output_dim=1, input_dim=self.xdim, init=init,
                                 W_regularizer=l2(w_l2),
                                 activity_regularizer=activity_l2(a_l2))) #100 works. 150 not
            self.model.add(Activation(activation))
            self.model.add(Dropout(self.dropout))
        return

    def compile(self):
        # different optimizers:
        # SGD = stochastic gradient descent 
        # loss: categorical_crossentropy / binary_crossentropy / crossentropy
        # for a mean squared error regression problem
        #self.model.compile(optimizer='rmsprop',
        #                        loss='mse'    , metrics=['accuracy'] )
        self.model.compile(loss='mean_squared_error', optimizer='adam', metrics=['accuracy'])
        #self.model.compile(loss='mean_squared_error', optimizer='adam')
        #self.model.compile(loss='categorical_crossentropy', optimizer='sgd', metrics=['accuracy'])
        #self.model.compile(loss='categorical_crossentropy', optimizer=SGD(lr=0.01, momentum=0.9, nesterov=True))
        return

    def fit(self, X_train, y_train, X_test=None, y_test=None):
        #model_out = self.model.fit(X_train, y_train, nb_epoch=400, batch_size=100) # each individual is seen n epoch times. 
        model_out = self.model.fit(X_train, y_train,
                                   nb_epoch=200,
                                   batch_size=32,
                                   validation_data=(X_test,y_test) ) # each individual is seen n epoch times. 
        if debug: print("model_out:", model_out)
        return

    def test(self, X_test, y_test):
        loss_and_metrics = self.model.evaluate(X_test, y_test, batch_size=32)
        if debug: print("loss_and_metrics:", loss_and_metrics)
        return loss_and_metrics

    def test_new( self, X_test_new):
        pred_y = self.model.predict(X_test_new)
        #if debug: print " in test_new; pred_y:", pred_y
        return pred_y

#### Call CINDES for predictions: ####
def main(X,y, X_new=[], fraction=0.0):
    #run_example(X,y)
    NN = neural(X,y, fraction)
    NN.run_procedure_INDES(X_new)
    return NN.scores




