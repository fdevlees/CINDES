
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

    def train(self, X=None, y=None, **kwargs):
        if X is None: X = self.X
        if y is None: y = self.y

        print "X.shape:", X.shape
        input_dim = X.shape[1]

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

        # other way:
        #model = get_model()
        #model.compile(loss='mse', optimizer='adam')
        #model.load_weights(modelname)
        return model

    def save_model(self, model, count):
        modelname = '{}_{}.hdf5'.format(self.name, count)
        print('saving model')
        model.save(modelname)
        return

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



