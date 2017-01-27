debug = 0

# python modules
import pprint

# my own modules
from CINDES4.predictor import learning_skl as learning
from CINDES4.predictor import learning_int as ml_int
from CINDES4.predictor import tfitter

from CINDES4.utils.writings import log_io
#import learning_skl as learning
#import learning_int as ml_int
#from writings import log_io
#import tfitter

def do_ml(indices, database, **TZmat):
    # via sklearn
    preds_ml = learning.machinelearning2(indices,database,**TZmat)
    # not via sklearn
    #preds_ml = learning.normal_machinelearning(indices=indices, table=database, **TZmat)
    return preds_ml

@log_io()
def predictor(run_object,table,mols_todo,mols_nodo,count, array=[]):
    '''makes the predictions using KRR(ML) / RR(LS) / DIF(MC)
       run_object = myrun with all param elements
    '''
    print "mols_todo:", mols_todo
    TZmat = run_object.TZmat
    predict = []
    enoughdata = not table==[] and not mols_todo==[] and count > 1

    if enoughdata:
      for prediction in run_object.predictions:
        # each prediction element is a dictionary with a 'type' key. 
        ptype = prediction['type']
        pred = prediction.copy()
        if ptype in ['ML','ml']:
            preds = do_ml(mols_todo, table, **TZmat)
        elif ptype=='1D':
            run_object.printlevel=1
            param = run_object.__dict__
            preds = tfitter.regression(table,mols_todo,**param)
        elif ptype=='2D':
            param = run_object.__dict__
            preds = tfitter.twodim_regression(table,mols_todo,**param)
        elif ptype=='MC':
            instance = tfitter.get_instance()
            preds = tfitter.dif_predict(mols_todo,maximum,instance)
        elif ptype=='iML':
            preds = ml_int.learn_int_skl_procedure(table,mols_todo, array, **TZmat)
        elif ptype=='NN':
            preds = learning.ANN(mols_todo, table, **TZmat)
        pred['results'] = dict( zip( mols_todo, preds ) )
        predict.append(pred)

    ##########

    # Splitting part
    #if run_object.ml==2 and not mols_todo==[]: #prescrean calculate only the best 50 %
    #    preds_data = zip(mols_todo, preds_ml) #get indices and predictions in same list
    #    preds_data = map(list,preds_data)
    #    preds_ml_sorted = sorted(preds_data,key=lambda x:x[1] ) #sort them based on prediction
    #    if run_object.optimum=='maximum': preds_ml_sorted = preds_ml_sorted[::-1] #when not optimum minimum reverse the list
    #    for item in preds_ml_sorted:
    #        item.insert(1,0)
    #    if debug#:
    #        print "preds_ml_sorted:"
    #        for item in preds_ml_sorted: print item
    #    half = len(preds_ml_sorted) / 2
    #    print "half:", half
    #    preds_tocal = preds_ml_sorted[:half]
    #    preds_nocal = preds_ml_sorted[half:]
    #    indices_tocal = [ item[0] for item in preds_tocal ] # the first entries are the indices
    #    data_nocal = data_nodo + preds_nocal # merge the already known and the predictions
    #elif run_object.regression==2:
    #    raise SystemExit('prescreaning not implemented')
    #    pass
    #else: #indentate off
    mols_nocal = mols_nodo
    mols_tocal = mols_todo
    if debug:
        print "data_nocal",mols_nocal
        try:
            #print "data_tocal",data_tocal
            print "indices_tocal",mols_tocal
        except NameError:
            print "NameError!"
    return mols_nocal, mols_tocal, predict





