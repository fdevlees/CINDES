import sys
import numpy as np
import pickle
sys.path.append('../data')
from read_table import MachineLearning, pearsonr
import os.path
from sklearn.model_selection import train_test_split
from sklearn.svm import SVR

data_file = '../data/BoB.pkl'

# cache data for loading speed
if  os.path.exists(data_file):
    print('loading data...')
    with open(data_file,'rb') as f:
        X,y = pickle.load(f)

else:
    print('reading data from table...')
    # read data
    ml = MachineLearning(type="BoB", data_file='../data/table_HLG.xyz')
    X,y = (ml.X,ml.y)
    # cache loaded data as pickle file
    with open(data_file,'wb') as f:
        pickle.dump((X,y),f,-1)

print(X.shape)
print(y.shape)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# svr with rbf kernel. C controls simplisity or decision surface. High C will
# try to fit all data and select more support vector. Low C will give a more
# smooth surface. gamma parameter defines how far the influence of a single
# training example reaches, with low values meaning far and high values
# meaning close. low gamma value high bias, high values high variance.
svr_rbf = SVR(kernel='rbf', C=1e4, gamma=1e-6)
print "Fitting..."
svr_rbf.fit(X_train, y_train)
#print svr_rbf.score(X_train,y_train)
#print svr_rbf.score(X_test, y_test)
print('training R: ' +str(pearsonr(y_train,svr_rbf.predict(X_train).flatten())))
print('test R: ' +str(pearsonr(y_test,svr_rbf.predict(X_test).flatten())))
