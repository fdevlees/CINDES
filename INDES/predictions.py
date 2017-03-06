debug = 0

# python modules
import pprint

# my own modules
from CINDES4.predictor import learning_skl as learning
from CINDES4.predictor import learning_int as ml_int
from CINDES4.predictor import tfitter

from CINDES4.utils.writings import log_io, print_title, dump
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


def get_experiment(prediction, table, run, retrain=True, array=[]):
    # each prediction element is a dictionary with a 'type' key. 
    ptype = prediction['type']



    kwargs = prediction # prediction is a dictonary with options specific for that experiment type
    kwargs['tableindex']=1
    
    if ptype == 'ml':
        #preds = do_ml(mols_todo, table, **TZmat)
        pass
    elif ptype=='1d':
        print "kwargs:", kwargs

        from CINDES4.predictor.linreg import LinRegOneExperiment
        regressor = LinRegOneExperiment(    table=table,
                                            retrain=retrain,
                                            array=array,
                                            run=run,
                                            **kwargs
                                            )
    elif ptype=='2d':
        regressor = LinRegOneExperiment(    table=table,
                                            retrain=retrain,
                                            array=array,
                                            run=run,
                                            **kwargs
                                            )
    elif ptype=='mc':
        #instance = tfitter.get_instance()
        #preds = tfitter.dif_predict(mols_todo,maximum,instance)
        pass
    elif ptype=='iml':
        #preds = ml_int.learn_int_skl_procedure(table,mols_todo, array, **TZmat)
        pass
    elif ptype=='nn':
        from CINDES4.predictor.nn import NeuralNetworkExperiment

        regressor = NeuralNetworkExperiment(   table=table,
                                               retrain=retrain,
                                               array=array,
                                               run=run,
                                               **kwargs
                                               )
    elif ptype=='gp':
        from CINDES4.predictor.gp import GaussianProcessExperiment, GaussianProcessWithPCAExperiment
        if prediction['pca']:
            regressor = GaussianProcessWithPCAExperiment(table=table,
                                                         retrain = retrain,
                                                         array=array,
                                                         run = run,
                                                         **kwargs )
        else:
            regressor = GaussianProcessExperiment(     table=table,
                                                       retrain = retrain,
                                                       array=array,
                                                       run = run,
                                                       **kwargs )
    elif ptype=='knn':
        #from CINDES4.predictor.knn import NearestNeighborWithPCAExperiment
        from CINDES4.predictor.knn import NearestNeighborExperiment, NearestNeighborWithPCAExperiment
        if prediction['pca']:
            print "PCA!"
            regressor = NearestNeighborWithPCAExperiment( table = table,
                                                          n_principal_components=100,
                                                          retrain=retrain,
                                                          array=array,
                                                          run = run,
                                                          **kwargs   #run=run
                                                        )
        else:
            regressor = NearestNeighborExperiment(    table = table,
                                                      retrain=retrain,
                                                      run = run,
                                                      array=array,
                                                      **kwargs   #run=run
                                                      )
    elif ptype=='svr':
        from CINDES4.predictor.svr import SupportVectorExperiment, SupportVectorWithPCAExperiment
        regressor = SupportVectorExperiment(          table=table,
                                                      n_principal_components=100,
                                                      retrain=retrain,
                                                      array=array,
                                                      run = run,
                                                      **kwargs   #run=run
                                                      )
    return regressor


@log_io()
def predictor(run,table,mols_todo,mols_nodo,count, nsite=0, array=[]):
    '''makes the predictions using KRR(ML) / RR(LS) / DIF(MC)
       run_object = myrun with all param elements
    '''
    made_pred=False
    print "mols_todo:", mols_todo
    TZmat = run.TZmat
    retrain = False
    enoughdata = not table==[] and not mols_todo==[] and count > 1

    if enoughdata:
      for prediction in run.predictions:
          # 0. log prediction:
          made_pred=True
          print_title(prediction['type'], outline='l')
          dump(prediction)


          # 1. initiate prediction experiment
          regressor = get_experiment(prediction = prediction,
                                     table= table,
                                     retrain= retrain,
                                     array=array,
                                     run=run)

          # 2. train or reload the model
          regressor.get_model(
                              write_log=True,
                              reopt_hyps=False,
                              count = count,
                              nsite = nsite,
                              )

          # 3. use model to predict
          regressor.predict(mols_todo)

      for molecule in mols_todo:
        print molecule.predictions

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
    #return mols_nocal, mols_tocal, predict
    return mols_nocal, mols_tocal, made_pred





