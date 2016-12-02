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
def predictor(run_object,table,indices_todo,data_nodo,count, array=[]):
    '''makes the predictions using KRR(ML) / RR(LS) / DIF(MC) '''
    print "indices_todo:", indices_todo
    TZmat = run_object.TZmat
    predict = {}

    # PREDICT VIA MACHINE LEARNING
    if run_object.ml==1 and not table==[] and not indices_todo==[]:
        preds_ml = do_ml(indices_todo, table, **TZmat)
        predict['ML'] = dict( zip( indices_todo, preds_ml ) )
    elif run_object.ml==2 and not table==[] and not indices_todo==[]:
        preds_ml = do_ml(indices_todo, table, **TZmat)
        predict['ML'] = dict( zip( indices_todo, preds_ml ) )
        if debug: print "preds_ml:", preds_ml
    else: preds_ml=[]

    # PREDICT VIA 1D REGRESSION
    if run_object.regression==1 and not table==[] and not indices_todo==[]:
        run_object.printlevel=1
        param = run_object.__dict__
        preds_ls = tfitter.regression(table,indices_todo,**param)
        predict['1D'] = dict( zip( indices_todo, preds_ls ) )
        # LOGGING 1
        print "Indices_todo & PREDICTIONS:"
        for index, prediction in zip(indices_todo, preds_ls):
            print index, prediction
    else: preds_ls=[]

    if run_object.tdregression==1 and not table==[] and not indices_todo==[]:
        param = run_object.__dict__
        preds_td = tfitter.twodim_regression(table,indices_todo,**param)
        predict['2D'] = dict( zip( indices_todo, preds_td) )

    # PREDICT VIA DIFMODEL
    if run_object.difmodel==1 and not table==[] and not indices_todo==[] and (run_object.restart>2 or count>1):
        instance = tfitter.get_instance()
        preds_dif = tfitter.dif_predict(indices_todo,maximum,instance)
        predict['MC'] = preds_dif

        # LOGGING 2
        print "Indices_todo & PREDICTIONS DIFMODEL:"
        for index, prediction in zip(indices_todo, preds_dif):
            print index, prediction

    else: preds_dif = []

    # OTHER PREDICTION SCHEME:
    if run_object.nosub==3 and count>2 and not indices_todo==[]:
        param = run_object.__dict__
        if False:
            preds_iml, preds_ml = ml_int.learn_int_procedure(table,indices_todo, array, **TZmat)
        else:
            preds_ml = do_ml(indices_todo, table, **TZmat)
            #preds_iml = ml_int.learn_int_procedure(table,indices_todo, array, **TZmat)
            preds_iml = ml_int.learn_int_skl_procedure(table,indices_todo, array, **TZmat)
        preds_ls = tfitter.regression(table,indices_todo,**param)
        #preds_ls = []
        preds_dif = tfitter.twodim_regression(table,indices_todo,**param)
        predict['1D'] = dict( zip( indices_todo, preds_ls ) )
        predict['2D'] = dict( zip( indices_todo, preds_dif) )
        predict['ML'] = dict( zip( indices_todo, preds_ml ) )
        predict['iML']= dict( zip( indices_todo, preds_iml) )

    # LOGGING 3
    if preds_ml==[] and preds_ls==[]:
        print "indices_todo: "
        pprint.pprint(indices_todo)

    if preds_ml==[] and preds_dif == [] and preds_ls == [] :
        predictions = []
    else:
        predictions = [indices_todo, preds_ml, preds_ls, preds_dif ]
    ##########

    # Splitting part
    if run_object.ml==2 and not indices_todo==[]: #prescrean calculate only the best 50 %
        preds_data = zip(indices_todo, preds_ml) #get indices and predictions in same list
        preds_data = map(list,preds_data)
        preds_ml_sorted = sorted(preds_data,key=lambda x:x[1] ) #sort them based on prediction
        if run_object.optimum=='maximum': preds_ml_sorted = preds_ml_sorted[::-1] #when not optimum minimum reverse the list
        for item in preds_ml_sorted:
            item.insert(1,0)
        if debug:
            print "preds_ml_sorted:"
            for item in preds_ml_sorted: print item
        half = len(preds_ml_sorted) / 2
        print "half:", half
        preds_tocal = preds_ml_sorted[:half]
        preds_nocal = preds_ml_sorted[half:]
        indices_tocal = [ item[0] for item in preds_tocal ] # the first entries are the indices
        data_nocal = data_nodo + preds_nocal # merge the already known and the predictions
    elif run_object.regression==2:
        raise SystemExit('prescreaning not implemented')
        pass
    else:
        data_nocal = data_nodo
        indices_tocal = indices_todo
    if debug:
        print "data_nocal",data_nocal
        try:
            #print "data_tocal",data_tocal
            print "indices_tocal",indices_tocal
        except NameError:
            print "NameError!"
    return data_nocal, indices_tocal, predict





