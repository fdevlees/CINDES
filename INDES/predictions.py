# python modules
import pprint
import logging

# my own modules
from CINDES.predictor import learning_skl as learning
from CINDES.predictor import learning_int as ml_int
from CINDES.utils.writings import log_io, print_title, dump
from CINDES.utils.utils import processify

debug = 0


def do_ml(indices, database, **TZmat):
    # via sklearn
    preds_ml = learning.machinelearning2(indices, database, **TZmat)
    # not via sklearn
    #preds_ml = learning.normal_machinelearning(indices=indices, table=database, **TZmat)
    return preds_ml


def get_experiment(prediction, table, run, retrain=True, array=[]):
    # each prediction element is a dictionary with a 'type' key.
    ptype = prediction['type']
    table = list(table.iteritems())
    #print "table:", table
    # prediction is a dictonary with options specific for that experiment type
    kwargs = prediction
    kwargs['tableindex'] = 1
    if ptype == 'ml':
        #preds = do_ml(mols_todo, table, **TZmat)
        pass
    elif ptype == '1d':
        #print "kwargs:", kwargs
        from CINDES.predictor.linreg import LinRegOneExperiment, LinRegOneWithPCAExperiment
        if prediction['pca']:
            regressor = LinRegOneWithPCAExperiment(table=table,
                                                   retrain=retrain,
                                                   array=array,
                                                   run=run,
                                                   **kwargs
                                                   )

        else:
            regressor = LinRegOneExperiment(table=table,
                                            retrain=retrain,
                                            array=array,
                                            run=run,
                                            **kwargs
                                            )
    elif ptype == '2d':
        regressor = LinRegOneExperiment(table=table,
                                        retrain=retrain,
                                        array=array,
                                        run=run,
                                        **kwargs
                                        )
    elif ptype == 'mc':
        #instance = tfitter.get_instance()
        #preds = tfitter.dif_predict(mols_todo,maximum,instance)
        pass
    elif ptype == 'iml':
        #preds = ml_int.learn_int_skl_procedure(table,mols_todo, array, **TZmat)
        pass
    elif ptype == 'nn':
        from CINDES.predictor.nn import NeuralNetworkExperiment

        regressor = NeuralNetworkExperiment(table=table,
                                            retrain=retrain,
                                            array=array,
                                            run=run,
                                            **kwargs
                                            )
    elif ptype == 'gp':
        from CINDES.predictor.gp import GaussianProcessExperiment, GaussianProcessWithPCAExperiment
        from CINDES.predictor.gp import GaussianProcessExperiment_skl
        if prediction['pca']:
            regressor = GaussianProcessWithPCAExperiment(table=table,
                                                         retrain=retrain,
                                                         array=array,
                                                         run=run,
                                                         **kwargs)
        else:
            regressor = GaussianProcessExperiment_skl(table=table,
                                                      retrain=retrain,
                                                      array=array,
                                                      run=run,
                                                      **kwargs)
    elif ptype == 'knn':
        #from CINDES.predictor.knn import NearestNeighborWithPCAExperiment
        from CINDES.predictor.knn import NearestNeighborExperiment, NearestNeighborWithPCAExperiment
        if prediction['pca']:
            print "PCA!"
            regressor = NearestNeighborWithPCAExperiment(table=table,
                                                         retrain=retrain,
                                                         array=array,
                                                         run=run,
                                                         **kwargs  # run=run
                                                         )
        else:
            regressor = NearestNeighborExperiment(table=table,
                                                  retrain=retrain,
                                                  run=run,
                                                  array=array,
                                                  **kwargs  # run=run
                                                  )
    elif ptype == 'svr':
        from CINDES.predictor.svr import SupportVectorExperiment, SupportVectorWithPCAExperiment
        regressor = SupportVectorExperiment(table=table,
                                            retrain=retrain,
                                            array=array,
                                            run=run,
                                            **kwargs  # run=run
                                            )
    elif ptype == 'krr':
        from CINDES.predictor.krr import KernelRidgeExperiment, KernelRidgeWithPCAExperiment
        if prediction['pca']:
            regressor = KernelRidgeWithPCAExperiment(table=table,
                                                     retrain=retrain,
                                                     array=array,
                                                     run=run,
                                                     **kwargs  # run=run
                                                     )
        else:
            regressor = KernelRidgeExperiment(table=table,
                                              retrain=retrain,
                                              array=array,
                                              run=run,
                                              **kwargs  # run=run
                                              )
    elif ptype == 'qml':
        from CINDES.predictor.krr_qml import KernelRidgeExperiment
        regressor = KernelRidgeExperiment(table=table,
                                          retrain=retrain,
                                          array=array,
                                          run=run,
                                          **kwargs  # run=run
                                          )
    return regressor


# @processify
@log_io(print_time=True)
def do_prediction(
        prediction,
        table,
        retrain,
        array,
        count,
        nsite,
        run,
        mols_todo):
    # prediction in: prediction, table, mols_todo, retrain, run, array, count,
    # nsite

    #raise SystemExit('predictions are not jet JSON implemented')
    # table has changed format is dict and has not automatically the properties that
    # must be predicted (for example function values are not included in table
    # so they must be calculated somewhere from the other items in table

    # 0. log prediction:
    made_pred = True
    print_title(prediction['name'], outline='l')
    dump(prediction)

    # 1. initiate prediction experiment
    regressor = get_experiment(prediction=prediction,
                               table=table,
                               retrain=retrain,
                               array=array,
                               run=run)

    # 2. train or reload the model
    regressor.get_model(
        write_log=True,
        reopt_hyps=False,
        count=count,
        nsite=nsite,
        n_max_train=None
    )

    # 3. use model to predict
    if not run.procedure == 'testpred':
        regressor.do_predict(mols_todo, rstd=True)

    # 4. add R to prediction and plots?
    prediction['R'] = regressor.R
    prediction['retrained'] = regressor.retrained
    try:
        prediction['plot1'] = regressor.plot1
    except AttributeError:
        print "prediction", prediction['name'], "has no plot1 attribute"
    try:
        prediction['plot2'] = regressor.plot2
    except AttributeError:
        print "prediction", prediction['name'], "has no plot2 attribute"
    try:
        prediction['plot3'] = regressor.plot3
    except AttributeError:
        print "prediction", prediction['name'], "has no plot3 attribute"
    return prediction


def do_prediction_process(*args, **kwargs):
    import concurrent.futures
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        result = executor.submit(do_prediction, *args, **kwargs).result()
    return result


def do_prediction_process2(*args, **kwargs):
    import multiprocessing
    with concurrent.futures.ProcessPoolExecutor(max_workers=1) as executor:
        result = executor.submit(do_prediction, *args, **kwargs).result()
    return result


def json_predictions(predictions):
    import json
    from CINDES.predictor.experiment_interface import jsonify

    #print "predictions:", predictions

    # try to load old predictions file:
    try:
        with open('predictions.json') as f:
            pred_dict = json.load(f)
    except IOError:
        pred_dict = {}

    for prediction in predictions:
        pred_dict[prediction['name']] = jsonify(prediction)

    with open('predictions.json', 'w') as fout:
        json.dump(pred_dict, fout, indent=4, sort_keys=True)

    return


def plot_predictions(predictions, prop=""):
    import matplotlib.pyplot as plt
    from numpy import asarray as A
    # if predictions[0]['plots']==[]:return
    plotted = True
    n = len(predictions)
    nx, ny = set_nxy(len(predictions))
    print nx, ny
    f, axs = plt.subplots(ny, nx)
    try:
        axs2d = [item for sublist in axs for item in sublist]
    except TypeError:
        try:
            axs2d = [item for item in axs]
        except TypeError:
            axs2d = [axs]
    print "axs2d:", axs2d
    for i, prediction in enumerate(predictions):
        if prediction['retrained'] or True:
            j = 0
            a = axs2d[i]
            print prediction['name']

            # p#lot test data
            try:
                #print prediction['plot1']
                x1 = A(prediction['plot1'])[:, 0]
                y1 = A(prediction['plot1'])[:, 1]
            except KeyError:
                print "no plot1 Key"
                j += 1
            else:
                a.scatter(
                    x1,
                    y1,
                    label=prediction['name'] +
                    " test",
                    s=6,
                    alpha=.8)
                plotted = True

            # plot train data
            try:
                #print prediction['plot2']
                x2 = A(prediction['plot2'])[:, 0]
                y2 = A(prediction['plot2'])[:, 1]
            except KeyError:
                print "no plot2 Key"
                j += 1
            else:
                a.scatter(
                    x2,
                    y2,
                    label=prediction['name'] +
                    " train",
                    s=6,
                    alpha=.5)
                plotted = True

            if j == 2:
                continue

            a.legend()
    if plotted:
        # pass
        plt.show()
    else:
        del f, axs, axs2d
        # plt.close()
        import gc
        gc.collect()
    return


def set_nxy(n):
    xys = {'1': (1, 1), '2': (2, 1), '3': (3, 1),
           '4': (2, 2), '5': (3, 2), '6': (3, 2),
           '7': (4, 2), '8': (4, 2), '9': (3, 3),
           '10': (4, 3), '11': (4, 3), '12': (4, 3)
           }
    return xys[str(n)]


@log_io()
def predictor(
        run,
        table,
        mols_todo,
        mols_nodo=[],
        count=0,
        nsite=0,
        array=[],
        optimum=None):
    '''makes the predictions using KRR(ML) / RR(LS) / DIF(MC)
       run_object = myrun with all param elements
    '''

    made_pred = False
    #print "mols_todo:", mols_todo
    retrain = nsite == 0
    retrain = False
    GA = False

    # there is enough data to do predictive analytics when:
    #    - there are at least 50 samples
    #    - each site is at least visited once: COUNT>1 OR run.restart is larger than zero
    #    - there should at least one prediction be made: not mols_todo==[] unless we are not interested in predictions: procedure testpred
    enoughdata = (
        (
            len(table) > 30 and (
                count > 1 or run.restart > 0
            )
        ) and
        not mols_todo == []
    ) or run.procedure == 'testpred'

    # make every item in run uncallable to be able to be pickled by the
    # subprocess.Queue
    store_function = run.function
    run.function = 'function'

    if run.predictions and enoughdata:
        # 1. do predictions
        for prediction in run.predictions:
            #do_prediction_process(prediction, table, retrain, array, count, nsite, run, mols_todo)
            prediction = do_prediction(
                prediction,
                table,
                retrain,
                array,
                count,
                nsite,
                run,
                mols_todo)
            made_pred = True

        # 2. get best R and set run.best_pred
        best_pred = max(run.predictions, key=lambda x: x['R'])
        run.best_pred = best_pred
        print_title("best performing estimator: " +
                    best_pred['name'] +
                    " with R**2: " +
                    str(best_pred['R']), outline='l')

        # 3.1 log
        for molecule in mols_todo:
            print molecule.predictions
        # 3.2. json_log
        json_predictions(run.predictions)
        # 3.3. do plottings
        # plot_predictions(run.predictions)

        # 4. decide which molecules to calculate and which not
        if run.procedure == 'ga':
            GA = True
            # raise NotImplementedError('yet todo')
            # what it should do:
            # if pred>0.75. take all structures performing better than the current optimum.
            # if that number is smaller than the ntake. take the n other best structures
            #
            #
        if best_pred['R'] > 0.75 and run.ml:
            def get_ntake(R, n):
                fraction = 4. - 4. * R
                ntake = max(1, int(fraction * n))
                print "I will calculate only {:d} of the {:d} structures ;)".format(
                    ntake, n)
                return ntake
            ntake = get_ntake(best_pred['R'], len(mols_todo))
            print "prediction is good enough"
            mols_nocal, mols_tocal = (mols_nodo, [])
            mols_sorted = sorted(
                mols_todo, key=lambda x: x.predictions.get(
                    best_pred['name']), reverse=(
                    not run.optimum == 'minimum'))
            if GA:
                print 'there could be a problem in the GA when there are predicted values below the actual calculated ones!'
                # take here the ones better than optimum?
            for i, mol in enumerate(mols_sorted):
                print mol, mol.predictions.get(best_pred['name'], "empty")
                if i < ntake:
                    mols_tocal.append(mol)
                else:
                    mol.Pvalue = mol.predictions.get(best_pred['name'], None)
                    mol.predicted = True
                    mols_nocal.append(mol)
            print "mols_tocal:\n", mols_tocal, "\nmols_nocal:\n", mols_nocal

        else:
            mols_nocal = mols_nodo
            mols_tocal = mols_todo
    else:
        logging.info("    no predictions will be made    len(table)={}".format(
            len(table)))
        mols_nocal = mols_nodo
        mols_tocal = mols_todo

    if debug:
        print "data_nocal", mols_nocal
        try:
            #print "data_tocal",data_tocal
            print "indices_tocal", mols_tocal
        except NameError:
            print "NameError!"
    # return mols_nocal, mols_tocal, predict
    run.function = store_function
    return mols_nocal, mols_tocal, made_pred

