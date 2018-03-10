#!/bin/env python
'''python learn.py '''

import sys
import numpy as np

from CINDES.predictor import learning_int
from CINDES.predictor import learning_skl as learning
from sklearn.kernel_ridge import KernelRidge

import matplotlib.pyplot as plt

def main_single(show=True,
                validation_fraction=0.2 ):
    ''' this function is written for the perspective of having a file like table.xyz and you want to do
        1. a grid_search for optimizing the hyperparameters of the estimator
        2. test the best estimator on the validation set
        3. analyse the results on the validation set
    '''
    from CINDES.predictor import learning_xyz

    # my_ML sets the data into the ML class. self.X / self.y and possibly self.X_val / self.y_val
    my_ML = learning_xyz.MachineLearning( inputfile = args.file,
                                          descriptor= 'norm4',
                                          validation_fraction = validation_fraction )

    # going to set the estimator: 
    gamma = 1. / (2* args.sigma**2 )
    print "sigma:", args.sigma
    print "gamma:", gamma
    clf = KernelRidge( alpha = args.labda, kernel= 'rbf' , gamma = gamma)
    my_ML.set_estimator(clf)
    print my_ML


    if True: # make True for a gridsearch
        param_grid = {'alpha' : [1e4, 1e2, 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5,1e-6, 1e-8 ],
                      'gamma' : np.logspace(-8, 8, 10) }
        best_params = my_ML.GridSearch( param_grid = param_grid)
    else:
        best_params = { 'labda': args.labda,
                        'sigma': args.sigma
                      }

    if True: # 
        # test on validation set
        y_val_pred = my_ML.predict()
        print "y_val_pred:", y_val_pred


    if True: # analysis of results
        rmse = learning_xyz.RMSE(     y_val_pred, my_ML.y_val)
        mae  = learning_xyz.MAE(      y_val_pred, my_ML.y_val)
        pearsonr =  learning_xyz.pearsonr( y_val_pred, my_ML.y_val)
        print "RMSE     :", rmse
        print "MAE      :", mae
        print "pearsonr :", pearsonr
        scores = [ rmse, mae, pearsonr ]
    else: scores = None

    if False: # plot predicted values vs real values
        learning_xyz.plot_error( my_ML.y_val, y_val_pred, 'bo', alpha=0.5 )
        if True: # plot also the training set in red
            y_pred = my_ML.predict( my_ML.X )
            learning_xyz.plot_error( my_ML.y, y_pred, 'ro', alpha=0.5 )
        if show: plt.show()

    if True: # plot KDE
        logdens = learning_xyz.KDE( my_ML.y_val, y_val_pred )
        if show: plt.show()
        if True: # save KDE
            outname = 'logdens' + str(args.validation_fraction*10) + str(args.sigma) + str(args.labda) + '.b'
            import pickle
            with open(outname,'wb') as fid:
                pickle.dump(logdens,fid)
    else: logdens=None

    return logdens, scores


def main_multiple():
    '''this function will do machine learning in multiple fractions combine the results and plot final results '''
    validation_fractions = [ 0.1, 0.3, 0.5, 0.7, 0.9 ]

    kdes = []
    for fraction in validation_fractions:
        logdens, scores = main_single(show=False, validation_fraction = fraction )
        kdes.append(logdens)
    print "kdes:", kdes

    fig, ax = plt.subplots()
    xlim=20
    X_plot = np.linspace(-xlim, xlim, 1000)[:, np.newaxis]
    for logdens in kdes:
        y = np.exp(logdens)
        ax.plot(X_plot, y, '-')
    plt.show()
    return



if __name__=='__main__':
    class Unbuffered(object):
        def __init__(self,stream):
            self.stream = stream
        def write(self,data):
            self.stream.write(data)
            self.stream.flush()
        def __getattr__(self,attr):
            return getattr(self.stream, attr)
    sys.stdout = Unbuffered(sys.stdout)
    import argparse
    parser = argparse.ArgumentParser(description="reads cycles data stored in cyclesinfo")
    parser.add_argument('file', nargs='?', default='table.xyz',help='filename default is table.xyz')
    parser.add_argument("-N","--neural",action="store_true",help="going to use Neural Networks")
    parser.add_argument("-S","--use_sklearn",action="store_true",help="analyze and make a property vs property plot of the training data")
    parser.add_argument("-M","--multiple",action="store_true",help="do multiple fractions and check if learning curve")
    parser.add_argument("-T","--timer",action="store_true",help="perform some time analyses")
    parser.add_argument("-k","--kernel",action="store",type = str,default='gaussian',help="which kernel to use: (laplacian, gaussian)")
    parser.add_argument("-d","--descriptor",action="store",type = str,default='norm4',help="which descriptor to use: (BoB, Coulomb(norm1/norm3/norm4))")
    parser.add_argument("-r","--random",action="store_true",help="use a randomly selected test set and training set")
    parser.add_argument("-s","--sigma",action="store",nargs='?',type=float,default=1.e2,const=1e7,help="do a sigma default 1e2 KRR")
    parser.add_argument("-B","--bandwidth",action="store",nargs='?',type=float,default=0.5,const=1e7,help="do a sigma default 1e2 KRR")
    parser.add_argument("-c","--cutoff",nargs=2, type = float, help="cutoff values min max")
    parser.add_argument("-l","--labda",action="store",nargs='?',type=float,default=1.e-5,const=1e-5,help="do a labda default 1e-5 KRR")
    parser.add_argument("-f","--fraction",action="store",nargs='?',type=float,default=1,const=1,help="between 0-1 use this fraction as training set")
    parser.add_argument("-V","--validation_fraction",action="store",nargs='?',type=float,default=0.0,const=0.0,help="between 0-1 use this fraction as a validation set")
    args=parser.parse_args()
    # get a test c,a,p

    #if args.timer:
    if False:
        print "use sklearn:", args.use_sklearn
        from CINDES.utils.timer import Timer
        with Timer() as t:
            learning.Amachinelearning2( sigma = args.sigma,
                               labda = args.labda,
                               fraction = args.fraction,
                               descriptor = args.descriptor,
                               kernel = args.kernel ,
                               args= args)
        print "=> elapsed learning5: %s s" % t.secs
    elif args.multiple:
        main_multiple()
    else:
        main_single(validation_fraction = args.validation_fraction)
    print "DONE LEARNING.PY"
