debug=False
#from writings import log_io, sprint
from CINDES4.utils.writings import log_io, print_title, sprint
from CINDES4.utils.utils import run_once
from copy import deepcopy
import pickle
import pprint
import time

import numpy as np
from scipy import stats
import pandas as pd

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

def log_table( mols, table, tablename='tablebin'):
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

    with open(tablename,'wb') as tfid: # write the table to a file
        pickle.dump(table,tfid)
        print "dumped table in {}".format(tablename)
    return table


def log_screen( mols ):
    for molecule in mols:
        item = [ molecule.index, molecule.Pvalue ]
        item.extend( molecule.boundaries )
        item.extend( molecule.infoline   )
        print formatitem(item)
    return

def log_screen_pred( mols ):
    preds = [ mol.predictions for mol in mols ]
    indices = [ mol.index for mol in mols ]
    pvalues = [ mol.Pvalue for mol in mols ]
    df = pd.DataFrame( preds, index = indices )
    df.insert(0,'pvalues', pvalues)
    print df
    return df
    #for molecule in mols:
    #    #pd.DataFrame( [ a.p, b.p, c.p ], index = [ a.name, b.name, c.name ] )

def log_pred_info(pred_info, count, k, l):
    @run_once
    def print_header(pfid, header):
        pfid.write( header )
        pfid.write( '\n' )
        return
    # add columns count, k, l 
    pred_info['count'], pred_info['site'], pred_info['nsite'] = ( count, k, l)

    # save dataframe
    with open('predinfo','a') as pfid:
            # do only once: print header
            print_header( pfid, ' '.join(pred_info.columns.values) )
            # print predictions
            pfid.write( pred_info.to_csv( sep=' ', header=None, mode='a'))
    return

def pstats(predinfo):
    from CINDES4.utils import statistics
    import pprint
    #import statistics
    #print "predinfo:\n", pprint.pformat(predinfo)

    #print predinfo['pvalues'].corr( predinfo['knn'])
    #print predinfo['pvalues'].corr( predinfo['knn'], method='spearman')

    # get all the pearson coefficients:
    pearsonr = [ predinfo['pvalues'].corr(predinfo[str(ml)]) for ml in map(str,predinfo.columns[1:-3]) ]
    print "pearsonr:", pearsonr

    with open('PRs','a') as p:
        p.write( ' '.join( map(str,pearsonr) + map(str, predinfo.iloc[0,-3:]) ) )
        p.write( '\n' )
    return


@log_io()
def loggings(mols,table,count,k,l, made_pred=False, tablename='tablebin'):

    #--- LOGGINGS: CYCLESINFO
    log_cyclesinfo(mols, count,k,l)

    #---- LOGGINGS: TABLEBIN
    table = log_table( mols, table, tablename=tablename )

    if debug:
        print "mols.Pvalue:", mols[0].Pvalue
        print "mols.index:", mols[0].index
        print "mols.infoline:", mols[0].infoline

    #---- LOGGINGS: to screen

    log_screen( mols )

    #---- NEW LOGGINGS: PREDICTIONS
    if made_pred:
        pred_frame = log_screen_pred( mols )
        log_pred_info( pred_frame, count, k, l )
        if True:
            pstats(pred_frame)

    #-----
    print "TIME:", time.strftime("%d %B %Y %H:%M:%S")
    return table
