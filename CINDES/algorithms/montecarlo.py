""" MONTE CARLO PROCEDURE """
debug = 0

from CINDES.evaluation import construction as zcon
from CINDES.utils.writings import log_io, print_title, sprint, dump
from CINDES.molecule import Molecule
from CINDES.evaluation.predictor import learning_skl as learning
import random
random.seed(123)
from math import exp  # exp(x) returns e^x
#from operator import mul


def randomconf(subarray, maxconf=tuple(), nrandsites=0):
    '''this function makes a random configuration. choosing one sub for each site '''
    # version 20/01/2016

    arlen = len(subarray)  # is length of subarray
    if debug:
        sprint(10, subarray)
        print("number of changed sites:", nrandsites)
    while True:
        if nrandsites == 0:  # then choose a whole new configuration
            conf = []
            for i in range(arlen):
                conf.append(random.choice(subarray[i]))
        else:  # only change nrandsites
            conf = maxconf[:]  # start from same conf
            # this gives an error because range(arlen) seems to an integer.
            sitenumbers = list(range(arlen))
            rands = random.sample(sitenumbers, nrandsites)  # choose nrandsites
            for i in rands:
                newgroup = random.choice(subarray[i])
                # if True: # if we want to test if the group is really changed:
                #    pass
                conf[i] = newgroup

        if not conf == maxconf:
            break
        # note that here it is only tested that the configuration is not same as maxconf. not if really enough
        # sites were changed
    return conf


def generate(p):
    '''this generates zero or one on a probability of p'''
    return random.random() <= p


@log_io(print_time=True)
def ML_init(prediction, run, table, retrain=True, array=[], count=100, nsite=0, **kwargs):
    # prediction in: prediction, table, mols_todo, retrain, run, array, count, nsite

    # 0. log prediction:
    made_pred = True
    print_title(prediction['name'], outline='l')
    dump(prediction)

    # 1. initiate prediction experiment
    from predictions import get_experiment
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
    )

    return regressor


def ML_pred(regressor, conf, **kwargs):
    # 3. use model to predict
    mol = Molecule(conf)
    regressor.do_predict([mol], rstd=True, MC=True)
    #print mol.predictions
    return mol.predictions[regressor.name]


@log_io()
def montecarloprocedure(run, subarray, maxi, Dtable, **kwargs):
    # version 4/10/2015
    # MONTE CARLO PROCEDURE.
               #maxsite = montecarloprocedure(beta, array, maximum, table)
    # INPUT: beta - maximum - table - array
    # OUTPUT: maxsite
    fileparameters = run.__dict__

    print("Monte Carlo switched on!")

    # 1. set Metropolis criterium parameters
    T = fileparameters['montecarlo']
    kb = 8.6e-5  # boltzmann constant
    beta = 1.0 / (kb * T)
    print("intial temperature is: ", T)

    # 2. set additional initial parameters
    cmaximum = zcon.indtocon(maxi.index)  # maximum is Molecule instance. get conf without dihedral angles

    Tcount = 0  # temperature counter. to zero after increased.
    Rcount = 0  # number of random confs tested
    #Tcountmax = int(10 ** (float(1 + fileparameters['nrandsites'])))
    Tcount = 0
    Tcountmax = 1000
    print("Number of tested configurations per temperature:", Tcountmax)

    # 3. FOR ML
    if fileparameters['ml'] > 0 or fileparameters['predictions']:
        if not hasattr(run, 'best_pred'):
            run.best_pred = run.predictions[0]
        ml_instance = ML_init(table=Dtable, array=subarray, run=run, prediction=run.best_pred, **kwargs)

    # 4. select random configurations until one is accepted.
    while True:
        TAKEOTHER = False
        # 4.1 select a random conf
        rconf = randomconf(subarray, cmaximum, fileparameters['nrandsites'])  # make a total random configuration
        rind = zcon.contoind(rconf)

        # 4.2a predict property via difference algorithm
        if not (fileparameters['ml'] > 0 or fileparameters['predictions']):
            deltaetje = 0
            for i in range(
                    len(rconf)):  # now we want to have a value erandom for this configuration and test it with a certain probability
                if not rconf[i] == cmaximum[i]:
                    confje = cmaximum[0:i] + [rconf[i]] + cmaximum[i + 1:]
                    indje = zcon.contoind(confje)
                    #print "Dtable[indje]:", Dtable[indje]
                    try:
                        deltaetje += Dtable[indje] - maxi.Pvalue
                    except KeyError as e:
                        print("indje:", indje)
                        print("error:", e)
                        print("maybe an ignored structure is encountered")
                        TAKEOTHER = True

            if TAKEOTHER:
                print("check other random structure")
                continue
                #print "deltaetje:", deltaetje
            erandom = float(maxi.Pvalue + deltaetje)

            if debug:
                print("random conf:", rconf)
                print("random ind:", rind)
                print("indje    :", indje)
                print("Dtable[indje]", Dtable[indje])
                print("deltaetje:", deltaetje)
                print("erandom:", erandom)

        else:
            # 4.2b predict property via MACHINE LEARNING
            erandom = ML_pred(ml_instance, rconf, **kwargs)
            if debug:
                print("ml_instance:", ml_instance)
                print("random index:", rind)
                print("erandom_ML:", erandom)

        # 4.3 determine acceptance based on Delta-P
        # calculate the gradient energy. > resulttry
        p = exp(- beta * (abs(erandom - maxi.Pvalue)))
        acceptance = generate(p)

        # 4.4
        Rcount += 1
        if rind == maxi.index:
            print("rconf similar to maxconf. not accepted")
            acceptance = 0
        if acceptance == 1:
            print("configuration accepted")
            print(Rcount, " configurations tested")
            print("Random Conf:", rind)
            print("property_random:", erandom)
            print("chance of acceptance:", p)
            print("final temperature while acceptance:", T)
            break
        else:
            #print "configuration not accepted"
            Tcount += 1
            if Tcount >= Tcountmax:  # once in hundred configurations, the temperature is increased by 10%
                Tcount = 0
                T = T * 1.1
                beta = 1.0 / (kb * T)
    mol = Molecule(conf=rconf)
    mol.predicted = True
    mol.Pvalue = erandom
    return mol
