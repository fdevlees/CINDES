from writings import print_title, log_io
import pickle
from CINDES.INDES import construction as zcon
from CINDES.INDES.procedures import Run, set_table
from CINDES.INDES import calculator
from CINDES.INDES.predictions import predictor
from CINDES.INDES.loggings import loggings

# procedure 7. Farthest Point Selection based on the diversity index. 
def database_construction(param,array):
    import numpy as np
    np.random.seed(41)
    debug=True
    myrun = Run(**param)
    myrun.divers_discardCH=False
    myrun.divers_minnch=True
    print(myrun) #this should print all the class elements via the __str__ function
    # the table with all the results of all calculated configs
    table = set_table(myrun)
    nmax = myrun.divers_nmax

    if not table:
        # run startind calculation
        table = None #make initial calculation. get table with length 1.

    count=1
    maxiter=100
    while True: #later while True
        print_title("COUNT: " + str(count),outline='l',signator="-")
        print "len table:", len(table)
        #1. get new structure(s) to calculate
        #if myrun.divindex==1:
            #1.1. run diversity on table and get occupancy per site. 
            #1.2. select for each site the least occuring group.
        mols = getdivers(array,table, myrun)
        with open('confsall','wb') as f:
            confsall = [ mol.index for mol in mols ]
            confsall.extend([item[0] for item in table])
            pickle.dump(confsall, f)
            print "new confsall written"
        if debug: print mols
        #1b check already in database
        mols_todo, mols_nodo = zcon.check_in_table(mols, table, myrun.props)
        #1c eventueel predictions
        mols_nocal, mols_tocal, made_pred = predictor(myrun, table, mols_todo,mols_nodo, count, array=array)
        #2. run new structure
        mols_all = calculator.procedure(myrun, mols_tocal, mols_nocal)
        #mols_all = submittingprocedure(mols_tocal,
        #                               mols_nocal,
        #                               myrun,
        #                             **myrun.TZmat     ) # here call submitting procedure
        #3. add structure to table
        if debug: print "after calculation:"
        if debug: print mols_all
        table = loggings(mols_all,table,count,1,1, made_pred=made_pred, tablename=myrun.tablename)
        if debug:
            for item in table: print item
        #4. stop if maxstructures is obtained. or other convergence criteria is met. 
        count += 1
        if len(table)>= nmax:
            print "desired number of samples reached!"
            break
        if count>maxiter:
            print "maxiterations reached"
            break
    return

@log_io()
def getdivers(array, table, myrun):
    '''give the next n most divers molecules

    This function should be replaced later to ... ?
    '''
    #0.
    index = myrun.divers_divindex
    batchsize = myrun.divers_batchsize
    discardCH = myrun.divers_discardCH
    minnch = myrun.divers_minnch

    #1.
    from CINDES.utils.molecule import Molecule
    from CINDES.utils.diversity1 import Diversifier
    from CINDES.utils.table import Tablebin
    diversifier = Diversifier(index=index)
    mols = []
    #confs = [ item[0].split('_') for item in table ]
    #confs = [ item.conf for item in table ]
    confs = table.keys()
    print "confs:", confs[:5]
    #seq = list(set([ group for conf in confs for group in conf ]))
    seq = [ 'CH','B','O','S','N','P',
                               'CNHH', 'CNOO','COH','CSH','CPh','CCHO','CSOOOH']
    print "seq:", seq
    once = True
    for _ in range(batchsize):
        if index==1 or index==2:
            occupancy_sum, occupancy = diversifier.get_occupancy12(seq, confs, discardCH=discardCH)
            print "in get divers: occupancy:", occupancy
            conf = make_molecule12(occupancy, seq, array, confs, index=myrun.divers_divindex)
        elif index==3:
            occupancy, occupancy_percentages = diversifier.get_occupancy3(seq, confs)
            if once:
                print "in get divers: occupancy:", occupancy
                once = False
            conf = make_molecule3(occupancy, seq, array, confs)
        elif index==30:
            occupancy, occupancy_percentages = diversifier.get_occupancy3(seq, confs)
            #print "in get divers: occupancy:", occupancy
            conf = make_molecule3_nsites(occupancy, seq, array, confs)
        if discardCH:
            from numpy.random import binomial, shuffle
            # adjust conf and place CH groups in it via a binomial distribution
            nsites=10
            while True:
                #nch = binomial(nsites+3,0.5) ###### here tuning factor. 
                nch = binomial(nsites,0.5) ###### here tuning factor. 
                if nch<nsites: break
            positions = range(nsites)
            shuffle(positions)
            conf = [ item if i<nch else 'CH' for item,i in zip(conf,positions)]
        elif minnch:
            # conf has to contain at least n ch groups.
            from numpy.random import shuffle
            nch = conf.count('CH')
            nmin= 4
            if nch<nmin:
                # place CH on certain groups
                #i_noch = [ i for i, group in enumerate(conf) if not group=='CH' ]
                i_noch = range(10)
                shuffle(i_noch)
                #print "i_noch", i_noch
                for i in range(5): conf[i_noch[i]]='CH'
                #conf = [ 'CH' if i in i_noch[:nmin-nch] else group for i, group in enumerate(conf) ]
        print conf
        fconf = zcon.indtocon( '_'.join(conf))
        mols.append(Molecule(conf=fconf))

        # now append to conf so div values can change.
        confs.append(conf)
    print "molecules:", mols
    return mols

def make_molecule12(occupancy, seq, array, confs, index=1):
    iarray = [ map("".join,item) for item in array ]
    if index==1:
        new_conf=[]
        for site_occ, site_array in zip(occupancy,iarray):
            #most_divers_group = min(zip(site_occ,seq))[1]
            # search for group that has the least occurance, but not zero because that indicates that is cannot occur on that site
            #most_divers_group = min(filter(lambda x:not x[0]==0.,zip(site_occ,seq)))[1]
            # see if in array for that site. 
            most_divers_group = min(filter(lambda x:x[1] in site_array,zip(site_occ,seq)))[1]

            new_conf.append(most_divers_group)
        if new_conf in confs:
            print "most divers mol already in table:"
            new_conf = take_nth(occupancy, seq, array, confs)
    return new_conf

def make_molecule3(occupancy, seq, array, confs):
    import numpy as np
    #print seq
    from collections import OrderedDict
    occD = dict()
    for group1, occ1 in zip(seq,occupancy):
        #print group1, occ1
        if group1 in ['S','O','CO']: continue #these are not possible so do not include
        for group2, occ12 in zip(seq,occ1):
            key = '{}_{}'.format(group1,group2)
            occD[key]=occ12
    occD = OrderedDict( sorted( occD.iteritems(), key=lambda x:x[1] ) )
    #print occD
    new_conf = ['']*len(confs[0])
    keys_visited=[]
    for ibond, gbond in zip(( (0,7),(1,5),(2,6),(3,9) ), occD.items()[:4]):
        keys_visited.append(gbond[0])
        i3,i2 = ibond
        g3,g2 = gbond[0].split('_')
        new_conf[i3]=g3
        new_conf[i2]=g2
    #now still two places have to be filled.
    # for bond4:
    # bond4 neighbors 1 and 3 so the tert group has to be g1 or g3
    for gbond in occD.items():
        if gbond[0] in keys_visited:continue
        g3,g2= gbond[0].split('_')
        if g3==new_conf[1] or g3==new_conf[3]:
            new_conf[4]=g2
            break
    else:
        raise StandardError
    #for bond8
    # bond 8 neighbors 2 and 3 so the tert group has to be g2 or g3
    for gbond in occD.items():
        if gbond[0] in keys_visited:continue
        g3,g2= gbond[0].split('_')
        if g3==new_conf[2] or g3==new_conf[3]:
            new_conf[8]=g2
            break
    else:
        raise StandardError
    #print new_conf
    assert not '' in new_conf, "one group not defined!"
    # now the change that a group appears on 4 8 is different from appearing on the others? so randomly symmetry permutation:
    Adasym = [
             [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10],
             [ 1, 3, 4, 2, 7, 8, 9,10, 5, 6],
             [ 1, 4, 2, 3, 9,10, 5, 6, 7, 8],
             [ 3, 2, 4, 1, 6, 7, 5, 9,10, 8],
             [ 4, 2, 1, 3, 7, 5, 6,10, 8, 9],
             [ 4, 1, 3, 2, 6,10, 8, 9, 7, 5],
             [ 2, 4, 3, 1,10, 5, 9, 7, 8, 6],
             [ 3, 1, 2, 4,10, 8, 6, 7, 5, 9],
             [ 2, 3, 1, 4, 9, 7, 8, 6,10, 5],
             [ 4, 3, 2, 1, 8, 9, 7, 5, 6,10],
             [ 3, 4, 1, 2, 5, 9,10, 8, 6, 7],
             [ 2, 1, 4, 3, 8, 6,10, 5, 9, 7] ]
    # randomly select one item from list
    new_order = Adasym[np.random.choice(range(10))]
    # permute new_conf according to new_order
    randomized = np.array(new_conf)[np.array(new_order)-1]
    return list(randomized)


def make_molecule3_nsites(occupancy, seq, array, confs, nsites=3):
    import numpy as np
    from collections import OrderedDict
    occD = dict()

    # make a dictionary of occupancies with siteT_siteS as keys
    for group1, occ1 in zip(seq,occupancy):
        if group1 in ['S','O','CO']: continue #these are not possible so do not include
        for group2, occ12 in zip(seq,occ1):
            key = '{}_{}'.format(group1,group2)
            occD[key]=occ12
    # order the dict so we have the least present combinations first
    occD = OrderedDict( sorted( occD.iteritems(), key=lambda x:x[1] ) )

    # make a new basic conf
    new_conf = ['CH']*len(confs[0])
    keys_visited=[]
    if nsites==3:
        # choose between 2-3-2 and 3-2-3
        # for adamantane ther are 4 2-3-2 possibilities and 6 3-2-3 possibilities so if we choose random:
        do_tert = np.random.rand()<0.4
        if do_tert:
            # choose two neighboring bonds 2-3-2
            bonds = ((0,5),(0,7))
        else:
            # choose two neighboring bonds 3-2-3
            bonds = ((0,5),(1,5))
        # take the first bond
        g3_1, g2_1 = occD.items()[0][0].split('_')
        i3_1, i2_1 = bonds[0]
        new_conf[i3_1]=g3_1
        new_conf[i2_1]=g2_1

        # now we take the second bond but or the second group or the first group has to be similar to the bond
        # that is already placed
        i=1 # the one but least occuring bond
        while True:
            g3_2, g2_2 = occD.items()[i][0].split('_')
            if do_tert:
                # the tertiary group has to match
                if g3_2==g3_1:
                    _, i2_2=bonds[1]
                    new_conf[i2_2]=g2_2
                    break
            else:
                # the secondary group has to match
                if g2_2==g2_1:
                    i3_2, _=bonds[1]
                    new_conf[i3_2]=g3_2
                    break
            i+=1
            print "i:", i, g3_2, g2_2
    else:
        raise SystemExit('stop no of sites not supported')

    print new_conf
    assert not '' in new_conf, "one group not defined!"
    # now we placed it on a specific 2-3-2 or 3-2-3 position after symmetry operations this is 
    # corrected
    Adasym = [
             [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10],
             [ 1, 3, 4, 2, 7, 8, 9,10, 5, 6],
             [ 1, 4, 2, 3, 9,10, 5, 6, 7, 8],
             [ 3, 2, 4, 1, 6, 7, 5, 9,10, 8],
             [ 4, 2, 1, 3, 7, 5, 6,10, 8, 9],
             [ 4, 1, 3, 2, 6,10, 8, 9, 7, 5],
             [ 2, 4, 3, 1,10, 5, 9, 7, 8, 6],
             [ 3, 1, 2, 4,10, 8, 6, 7, 5, 9],
             [ 2, 3, 1, 4, 9, 7, 8, 6,10, 5],
             [ 4, 3, 2, 1, 8, 9, 7, 5, 6,10],
             [ 3, 4, 1, 2, 5, 9,10, 8, 6, 7],
             [ 2, 1, 4, 3, 8, 6,10, 5, 9, 7] ]
    # randomly select one item from list
    new_order = Adasym[np.random.choice(range(10))]
    # permute new_conf according to new_order
    randomized = np.array(new_conf)[np.array(new_order)-1]
    return list(randomized)

def take_nth(occupancy, seq, array, confs):
    ''' it is possible that the previous function returns a molecule that already exists in the database
    so here i'll write a clever method to come up with a nth-but-most divers structure. '''
    iarray = [ map("".join,item) for item in array ]
    nsites = len(iarray)

    # the next line combines the seq and occupancy and filters only the ones that also are allowed on that site. i.e. they are in array
    sel_occseq = [ filter(lambda x:x[1] in iarray_site,zip(occ_site, seq)) for occ_site, iarray_site in zip(occupancy, iarray) ]
    # sort each element:
    sorted_occseq = [ sorted(item) for item in sel_occseq ]

    #make a list of possibly one-but-lasts. 
    nthbutbestconfs = []

    # loop over sites
    identity = [[0]*n + [1] + [0]*(nsites-n-1) for n in range(nsites) ]
    print "identity:", identity
    for row,sorted_occseq_site in zip(identity,sorted_occseq):
        divconf = [ sorted_occseq_site[i] for i in row ]
        nthbutbestconfs.append(divconf)
    print "nth but best confs:", nthbutbestconfs

    divvaluesnth = [ sum(zip(*item)[0]) for item in zip(*nthbutbestconfs) ]

    for item in sorted( zip(divvaluesnth, zip(*nthbutbestconfs))):
        conf = zip(*item[1])[1]
        if not conf in confs:
            return conf
    else:
        print("# every nth also in conf. not possible to find a good structure")
        SystemExit('stop')


