import json
from numpy import asarray as A
from CINDES.INDES.predictions import plot_predictions

# read file
with open('predictions.json') as f:
    data = json.load(f)

# process data
preds = list(data.values())
for pred in preds:
    pred['plot1'] = A( pred['plot1'] )
for pred in preds:
    pred['plot2'] = A( pred['plot2'] )

# plot
plot_predictions( preds )

# get MAEs and RMSEs
for pred in preds:
    print()
    print(pred['name'])

    #test errors
    tests = pred['plot1']
    MAE_test = sum(abs( tests[:,0] - tests[:,1] ) ) / len(tests)
    RMSE_test = sum( ( tests[:,0] - tests[:,1] )**2 ) / len(tests)
    print("MAE_test:", MAE_test)
    print("RMSE_test:", RMSE_test)

    #training errors
    trains = pred['plot2']
    MAE_train = sum(abs( trains[:,0] - trains[:,1] ) ) / len(trains)
    RMSE_train = sum( ( trains[:,0] - trains[:,1] )**2 ) / len(trains)
    print("MAE_train:", MAE_train)
    print("RMSE_train:", RMSE_train)

