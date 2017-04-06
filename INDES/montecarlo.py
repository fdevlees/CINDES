""" MONTE CARLO PROCEDURE """
debug=0

import construction as zcon
from CINDES4.utils.writings import log_io, print_title, sprint
from CINDES4.predictor import learning_skl as learning
import random
from math import exp #exp(x) returns e^x
#from operator import mul

def randomconf(subarray, maxconf=[], nrandsites=0):
    '''this function makes a random configuration. choosing one sub for each site '''
    # version 20/01/2016

    arlen = len(subarray) #is length of subarray
    if debug:
        sprint(10, subarray)
        print "number of changed sites:", nrandsites
    while True:
        if nrandsites == 0: #then choose a whole new configuration
            conf = []
            for i in range(arlen):
                conf.append(random.choice(subarray[i]))
        else: #only change nrandsites
            conf = maxconf[:] #start from same conf
            #this gives an error because range(arlen) seems to an integer.
            sitenumbers = range(arlen)
            rands = random.sample(sitenumbers,nrandsites) #choose nrandsites
            for i in rands:
                newgroup=random.choice(subarray[i])
                #if True: # if we want to test if the group is really changed:
                #    pass
                conf[i] = newgroup

        if not conf==maxconf: break
        # note that here it is only tested that the configuration is not same as maxconf. not if really enough 
        # sites were changed
    return conf

def generate(p):
    '''this generates zero or one on a probability of p'''
    return random.random() <= p

@log_io()
def montecarloprocedure(fileparameters, subarray, maxi, table,**kwargs):
    #version 4/10/2015
    #MONTE CARLO PROCEDURE. 
               #maxsite = montecarloprocedure(beta, array, maximum, table)
    # INPUT: beta - maximum - table - array
    # OUTPUT: maxsite

    print "Monte Carlo switched on!"

    # 1. set Metropolis criterium parameters
    T = fileparameters['montecarlo']
    kb = 8.6e-5 #boltzmann constant
    beta = 1.0 / ( kb * T )
    print "intial temperature is: ",T

    # 2. set additional initial parameters
    cmaximum = zcon.indtocon(maxi[0]) # maximum is in index format. change to confformat
    if debug:
        print "maxi:", maxi
        print "len(table):", len(table)
        for item in table: print item

    Dtable = dict([ (item[0],item[1]) for item in table] )
    Tcount = 0 #temperature counter. to zero after increased.
    Rcount = 0 #number of random confs tested
    Tcountmax = int ( 10** ( float( 1 + fileparameters['nrandsites'] ) ) )
    print "Number of tested configurations per temperature:", Tcountmax

    # 3. FOR ML
    if fileparameters['ml']==1:
        ml_instance = learning.MC_init(table, **kwargs)

    # 4. select random configurations until one is accepted. 
    print "Temperatures:",
    while True:
        # 4.1 select a random conf
        rconf = randomconf(subarray,cmaximum,fileparameters['nrandsites']) # make a total random configuration
        rind = zcon.contoind(rconf)

        # 4.2a predict property via difference algorithm
        deltaetje = 0
        for i in range(len(rconf)): # now we want to have a value erandom for this configuration and test it with a certain probability
            if not rconf[i] == cmaximum[i]:
                confje = cmaximum[0:i] + [rconf[i]] + cmaximum[i+1:]
                indje= zcon.contoind(confje)
                #print "Dtable[indje]:", Dtable[indje]
                try:
                    deltaetje += Dtable[indje] - maxi[2]
                except KeyError as e:
                    print "indje:", indje
                    print "error:", e
                    raise
                #print "deltaetje:", deltaetje
        erandom = float (maxi[2] + deltaetje)

        if debug:
            print "random conf:", rconf
            print "random ind:",  rind
            print "indje    :",   indje
            print "Dtable[indje]",Dtable[indje]
            print "deltaetje:", deltaetje
            print "erandom:", erandom

        # 4.2b predict property via MACHINE LEARNING
        if fileparameters['ml']==1:
            print "ml_instance:", ml_instance
            print "indje:", indje
            erandom_ML = learning.MC_test_ind(ml=ml_instance, indices=[indje],**kwargs) 
            print "erandom_ML:", erandom_ML

        # 4.3 determine acceptance based on Delta-P
        # calculate the gradient energy. > resulttry
        p = exp(- beta * (abs( erandom - maxi[1] )))
        acceptance = generate(p)

        # 4.4 
        Rcount +=1
        if rind == maxi[0]:
            print "rconf similar to maxconf. not accepted"
            acceptance = 0
        if acceptance == 1:
            print maxi[0],'\n',rind
            print "configuration accepted"
            print Rcount, " configurations tested"
            print "Random Conf:" , rind
            print "property_random:", erandom
            print "chance of acceptance:",p
            print "final temperature while acceptance:", T
            break
        else:
            #print "configuration not accepted"
            Tcount +=1
            if Tcount >= Tcountmax: #once in hundred configurations, the temperature is increased by 10%
                Tcount = 0
                T = T * 1.1
                beta = 1.0 / ( kb * T )
    return [ rind,0.0, erandom ]



