debug = False
#from writings import log_io, sprint
from CINDES.utils.writings import log_io, print_title, sprint
from CINDES.utils.utils import run_once, jsonify
from copy import deepcopy
#import pickle
import json
import pprint
import time

import numpy as np
from scipy import stats
import pandas as pd


def formatitem(opt, item, maxlenconf=49):
    def formatter(item):
        if isinstance(item, bool):
            return str(item).rjust(15)
        try:  # float
            return '{:15.8f}'.format(item)
        except ValueError:  # str
            if item is None:
                return " "*15
            try:
                return item.rjust(15)
            except AttributeError:
                return "No single value"
            # return "      {}".format(item)

    index = '{opt} {conf:{width}s}'.format(opt=opt, conf=item[0], width=maxlenconf + 1)
    abin = ' {:5s} '.format(str(item[1]))
    try:
        #datas = ' '.join(( fmt.format('{:15.8f}', (datatje for datatje in item[2:] ))))
        datar = [formatter(datatje) for datatje in item[2:]]
        datas = ' '.join(datar) + "  |"
    except ValueError:
        print item[2:]
        raise
    itemstring = index + abin + datas
    return itemstring


def log_cyclesinfo(mols, count, k, l):
    with open('cyclesinfo', 'a') as cfid:
        # filedata=deepcopy(data[:])
        for molecule in mols:
            item = [molecule.index]
            item.append(int(not molecule.predicted))  # 1 if really calculated 0 if only predicted

            # new json style:
            #print "molecule.props:", molecule.props
            item.append(molecule.Pvalue)  # predictions have only this one?

            def isSingleValue(x): return isinstance(x, int) or isinstance(x, float) or isinstance(x, str)
            item.extend(filter(isSingleValue, molecule.props.values()))

            # old pickle style:
            #item.append( molecule.Pvalue)
            #item.extend( molecule.boundaries )
            #item.extend( molecule.infoline )

            item.extend([count, k, l])
            cfid.write(' '.join(pprint.pformat(i) for i in item) + '\n')
    #del filedata
    return


def log_table(mols, table, tablename='table'):
    ''' this function updates the table-dictionary used inside the program
    and updates the database.json file.
    '''
    def update_json(table, filename):
        # filename has to end in .json. else .json is added.
        # this prevents overwriting an pickle style formatted inputtable
        if not filename[-5:] == '.json':
            filename = '{}.json'.format(filename)

        # 1 read json object:
        try:
            with open(filename) as f:
                json_table = json.load(f)
        except ValueError:
            print "JSON file empty."
            json_table = dict()

        # 2 update (nested)
        for key, value in table.iteritems():
            if key in json_table:
                json_table[key].update(value)
            else:
                json_table[key] = value
        json_table = jsonify(json_table)

        # 3 write updated json object
        with open(filename, 'w') as f:
            json.dump(json_table, f, indent=None)
        print "dumped table in {} with {} of the {} molecules".format(filename, len(table), len(json_table))
        return
    # -------------

    # move the new data to table except duplicates
    for molecule in mols:
        if molecule.predicted == False:
            # if not molecule.index in [ item[0] for item in table ]:
            if not molecule.index in table:
                table[molecule.index] = molecule.props
                #tableitem = [ molecule.index ]
                #tableitem.append( molecule.Pvalue     )
                #tableitem.extend( molecule.boundaries )
                #tableitem.extend( molecule.infoline   )
                #if debug: print "tableitem:", tableitem
                # table.append(tableitem)

    update_json(table, tablename)
    return table


def log_screen(mols):
    def issinglevalued(x):
        return any([isinstance(x, t) for t in (str, int, float, bool)])

    # first print the properties of the first molecule
    try:
        print "first molecule:", mols[0].index
        for k, v in mols[0].props.iteritems():
            print "{:15s}:".format(k),
            if len(repr(v))>100:
                print
                pprint.pprint(v, indent=2, width=160)
            else:
                print v
    except IndexError:
        return

    # get property line.
    # get all the props that possibly have to be printed
    # but are not a list/dict
    # NB there are predicted confs that only have a Pvalue so they have no props attribute
    props = set()
    noprintprops = set()
    for mol in mols:
        try:
            for prop, value in mol.props.items():
                if prop in props:
                    continue
                if issinglevalued(value):
                    props.add(prop)
                else:
                    noprintprops.add(prop)
            #props.update(mol.props.keys())
        except AttributeError:
            pass
    if isinstance(props, set):
        props = list(props)
    print "Properties that cannot be represented as a single value are:", ", ".join(noprintprops)

    # try to find the target property by looking for the Pvalue in all props such that
    # the target property is not printed twice.
    # if the property is not found, the Pvalue is probably a unique value obtained by 
    # a function of other properties
    j = 0
    p = None
    while j < len(mols):
        try:
            keys, values = zip(*filter(
                lambda x:issinglevalued(x[1]), mols[j].props.items()
                ))
            # check if Pvalue is one of these singular props
            if mols[j].Pvalue in values:
                # so yes. Pvalue is one of the propvalues. but which one?
                # get index of prop
                #keys, values= zip(*mols[0].props.items())
                i = values.index(mols[j].Pvalue)
                p = keys[i]
                # remove that one from props
                props.remove(p)
                #print "property:", p
            else:
                #print "function"
                p = 'function'
            break
        except ValueError:
            j += 1
    print "props:", props
    # print everything:
    # .1 decide max conf lenght
    maxlenconf = max(map(lambda x: len(x.index), mols))

    # decide how many props per line and how many table need to be printed
    propsets=[]
    maxnpropsperline = 7
    nextralines = len(props) / maxnpropsperline
    if nextralines:
        npropsperline = len(props) / (nextralines+1)
        for i in range(nextralines+1):
            propsets.append(props[i*npropsperline:(i+1)*npropsperline])
        # correction
        if not props[-1] in propsets[-1]:
            propsets[-1].append(props[-1])
    else:
        propsets=[props]

    print "propsets:", propsets

    # then print a table for every propset
    for i, propset in enumerate(propsets):
        lenh = maxlenconf + 26 + len(propset) * 16
        if i==0:
            headers =  "| index" + (maxlenconf - 4) * " " + " pred?   {:15s} ".format(p) + \
                " ".join(('{:15s}'.format(prop) for prop in propset)) + "|"
        else:
            lenh -= 16
            headers =  "| index" + (maxlenconf - 4) * " " + " pred?   ".format(p) + \
                " ".join(('{:15s}'.format(prop) for prop in propset)) + "|"
            print "    &"
        print "+{}+".format(lenh * "-")
        print headers
        print "}}{}{{".format(lenh * "-")
        # print data
        for molecule in mols:
            opt = "|"
            if molecule.opt:
                opt = "+"

            if i and hasattr(molecule, 'smiles'):
                item = [molecule.smiles, molecule.predicted]
            else:
                item = [molecule.index, molecule.predicted]
            if not i:
                item += [ molecule.Pvalue ]

            propvals = [molecule.props.get(prop, 'unknown') for prop in propset]
            item.extend(propvals)
            print formatitem(opt, item, maxlenconf)
        print "+{}+".format(lenh * "-")
    return p


def log_screen_pred(mols):
    print
    preds = [mol.predictions for mol in mols]
    indices = [mol.index for mol in mols]
    pvalues = [mol.Pvalue for mol in mols]
    df = pd.DataFrame(preds, index=indices)
    df.insert(0, 'pvalues', pvalues)

    # ----- pretty print df -----
    s = df.to_string().split('\n')
    ls = len(s[0])
    print "+{}+".format((ls + 2) * "-")
    print "| {} |".format(s[0])
    print "+{}+".format((ls + 2) * "-")
    for item in s[1:]:
        print "| {} |".format(item)
    print "+{}+".format((ls + 2) * "-")
    # ----- end pretty print-----

    return df
    # for molecule in mols:
    #    #pd.DataFrame( [ a.p, b.p, c.p ], index = [ a.name, b.name, c.name ] )


def log_pred_info(pred_info, count, k, l):
    @run_once
    def print_header(pfid, header):
        pfid.write(header)
        pfid.write('\n')
        return
    # add columns count, k, l
    pred_info['count'], pred_info['site'], pred_info['nsite'] = (count, k, l)

    # save dataframe
    with open('predinfo', 'a') as pfid:
            # do only once: print header
        print_header(pfid, ' '.join(pred_info.columns.values))
        # print predictions
        pfid.write(pred_info.to_csv(sep=' ', header=None, mode='a'))
    return


def pstats(predinfo):
    from CINDES.utils import statistics
    import pprint
    #import statistics
    #print "predinfo:\n", pprint.pformat(predinfo)

    #print predinfo['pvalues'].corr( predinfo['knn'])
    #print predinfo['pvalues'].corr( predinfo['knn'], method='spearman')

    # get all the pearson coefficients:
    pearsonr = [predinfo['pvalues'].corr(predinfo[str(ml)]) for ml in map(str, predinfo.columns[1:-3])]
    print "pearsonr:", pearsonr

    with open('PRs', 'a') as p:
        p.write(' '.join(map(str, pearsonr) + map(str, predinfo.iloc[0, -3:])))
        p.write('\n')
    return


@log_io()
def loggings(mols, table, count, k, l, made_pred=False, tablename='table'):

    # --- LOGGINGS: CYCLESINFO
    log_cyclesinfo(mols, count, k, l)

    # ---- LOGGINGS: TABLE.JSON
    table = log_table(mols, table, tablename=tablename)

    # ---- LOGGINGS: to screen
    log_screen(mols)

    # ---- NEW LOGGINGS: PREDICTIONS
    if made_pred:
        pred_frame = log_screen_pred(mols)
        log_pred_info(pred_frame, count, k, l)
        if True:
            pstats(pred_frame)

    # -----
    print "TIME:", time.strftime("%d %B %Y %H:%M:%S")
    return table
