debug=False
#from writings import log_io, sprint
from CINDES4.utils.writings import log_io, print_title, sprint
from copy import deepcopy
import pickle
import pprint
import time

import numpy as np
from scipy import stats

def formatitem(item):
    index = '{:50s}'.format(item[0])
    abin  = ' {} '.format(str(item[1]))
    datas = ' '.join(( '{:15.8f}'.format(datatje) for datatje in item[2:] ) )
    itemstring = index + abin + datas
    return itemstring

def log_cyclesinfo(mols, count, k, l):
    with open('cyclesinfo','a') as cfid:
        #filedata=deepcopy(data[:])
        for molecule in mols:
            item = [ molecule.index ]
            item.append( int(not molecule.predicted) ) # 1 if really calculated 0 if only predicted 
            item.append( molecule.Pvalue)
            item.extend( molecule.boundaries )
            item.extend( molecule.infoline )
            item.extend([count,k,l])
            cfid.write(' '.join(pprint.pformat(i) for i in item)+'\n')
    #del filedata
    return

def log_table( mols, table):
    if debug:
        print "in log_table: mols:", mols
        print "table:", table
    # here move the new data to table except duplicates
    for molecule in mols:
        if molecule.predicted == False:
            if not molecule.index in [ item[0] for item in table ]:
                tableitem = [ molecule.index ]
                tableitem.append( molecule.Pvalue     )
                tableitem.extend( molecule.boundaries )
                tableitem.extend( molecule.infoline   )
                if debug: print "tableitem:", tableitem
                table.append(tableitem)
    # OLD:
    #for item in data:
    #    if item[1]==1:
    #        if not item[0] in [tja[0] for tja in table]:
    #            tableitem = [item[0]] + item[2:]
    #            table.append(tableitem)
    #            if debug: print "tableitem:", tableitem
    #    else:
    #        assert item[1]==0, "item[1] has to be 1 or 0 but is %s" % str(item[1])

    with open('tablebin','wb') as tfid: # write the table to a file
        pickle.dump(table,tfid)
        print "dumped tablebin"
    return table


def log_screen( mols ):
    for molecule in mols:
        item = [ molecule.index, molecule.Pvalue ]
        item.extend( molecule.boundaries )
        item.extend( molecule.infoline   )
        print formatitem(item)
    return

    #datadict = dict( ( (item[0], item[1:]) for item in data) )
    #HEADER
    #print "index, value", ' '.join( item['type'] for item in predict )
    #for key, values in datadict.iteritems():
    #    print key, values[1],
    #    for prediction in predict:
    #        try:
    #            print prediction['results'][key],
    #        except KeyError:
    #            pass
    #    print
    return

#def get_pred_info( data_dict, predict ):
#    # from the predict types that are present 
#    pred_info = []
#    pred_types = sorted( predict.keys() )
#    for index in sorted( predict[pred_types[0]].keys() ) : # getting the indices of the first prediction values in keys # for index in ['1D','2D','Dif','ML']
#        index_info = []
#        # add first the real data item. so it is the first item. 
#        index_info.extend( [ index, data_dict[index][1] ] )  #index 1 is normally the optimization property. 
#        for pred_type in pred_types:
#            index_info.append( predict[ pred_type ] [index])
#        pred_info.append(index_info)
#    return pred_info

def get_pred_info( mols ):
    pred_info = []
    for molecule in mols:
        mol_info=[]
        for prediction in molecule.predictions:
            mol_info.append( prediction )
        pred_info.append(mol_info)

    #for index in predictions[0]['results'].keys() :
    #    index_info = []
    #    index_info.extend( [ index, data_dict[index][1] ])
    #    for prediction in predictions:
    #        index_info.append( prediction['results'][index])
    #    pred_info.append(index_info)

    return pred_info

def log_pred_info(pred_info, count, k, l):
    with open('predinfo','a') as pfid:
        for item in pred_info:
            item.extend([count,k,l])
            pfid.write(' '.join(pprint.pformat(i) for i in item) + '\n')
    return

def pstats(pred_info):
    from CINDES4.utils import statistics
    #import statistics
    print "pred_info:"
    transp = zip(*pred_info)
    n = len(transp)-3
    print "n6?:", n
    ps = [ stats.pearsonr(transp[1],transp[i]) for i in range(2,n) ]
    # calculate sorting scores by sum(abs( x(i) - y(i) ) )
    order_scores = []
    real_order = np.argsort( np.asarray( transp[1] ) )
    for i in range(2,n):
        score = statistics.order_score2(transp[i], transp[1] )
        score2 = statistics.order_score3(transp[i], transp[1] )
        score3 = statistics.order_score4(transp[i], transp[1] )
        score4 = statistics.order_score5(transp[1], transp[i] ) # here the order is important!
        order_scores.append(score)
        order_scores.append(score2)
        order_scores.append(score3)
        order_scores.append(score4)
    ps.append( order_scores )

    with open('PRs','a') as p:
        p.write( ' '.join( [ ' '.join( [ str(i) for i in item ] ) for item in ps ] ) )
        p.write( '\n')
    return


@log_io()
def loggings(mols,table,count,k,l,predict=[]):

    #--- LOGGINGS: CYCLESINFO
    log_cyclesinfo(mols, count,k,l)

    #---- LOGGINGS: TABLEBIN
    table = log_table( mols, table )

    if debug:
        print "mols.Pvalue:", mols[0].Pvalue
        print "mols.index:", mols[0].index
        print "mols.infoline:", mols[0].infoline

    #---- LOGGINGS: to screen

    log_screen( mols )

    #---- NEW LOGGINGS: PREDICTIONS
    if not predict == []:
        pred_info = get_pred_info( mols )
        log_pred_info( pred_info, count, k, l )
        if True:
            pstats(pred_info)

    #-----
    print "TIME:", time.strftime("%d %B %Y %H:%M:%S")
    return table
