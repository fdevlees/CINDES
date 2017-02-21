''' module for getting data.xyz and descriptors '''

from CINDES4.utils.converter import Converter
from CINDES4.INDES import construction as zcon

import numpy as np
from copy import deepcopy
from descriptor import get_X_1D

def get_XY(table, core,active,passive, tableindex=2, descriptor='BoB', identify=''):
    ''' calculte X and y '''
    def get_y(table, tableindex):
        y =  np.fromiter((item[tableindex] for item in table ),np.float)
        return y

    #### MAKE Y
    # get input from inputfile table
    y = get_y(table, tableindex)

    #### MAKE X
    ## X.1: get indices from table
    indices = (item[0] for item in table)
    if '1D' in descriptor:
        X = get_X_1D(indices, descriptor, identify)
    else:
        X = get_X( indices, core=core, active=active, passive=passive)
    return X,y

def get_X(indices, descriptor='BoB', **TZmat):

    # 1. convert new indices to confs to ZMAT
    converter = Converter()
    mats = tuple( contozma(zcon.indtocon(item),**TZmat) for item in indices)
 
    # 2. convert new ZMATs to XYZs
    converter = Converter()
    try:
        xyzs = [ zmatoxyz(converter,item) for item in mats ]
    except KeyError:
        print "Error: with:", item
        i = mats.index(item)
        print "index:", table[i]
        raise
 
    ## X.3: convert cartesian coordinates to descriptor
    if descriptor=='BoB':
        X = np.asarray( tuple( BoB(item) for item in xyzs) )
    else:
        X = np.asarray( tuple( coulomb(item) for item in xyzs) )
    return X



def contozma(conf,core,active,passive):
    c = deepcopy(core)
    a = deepcopy(active)
    p = deepcopy(passive)
    mat = zcon.constructor2(conf,c,a,p)
    return mat

def zmatoxyz(converter,mat):
    ''' convert a zmat to xyz coordinates via the Converter instance '''
    zmat = converter.read_zmalist(mat)
    return converter.zmatrix_to_cartesian()

def BoB(xyz):
    d = False
    from collections import OrderedDict
    l = len(xyz)
    Bag = []
    # max no of atoms per type
    typef = OrderedDict((( 'H' , 36 ),
                         ( 'C' , 20 ),
                         ( 'O' , 20 ),
                         ( 'N' , 10 ),
                         ( 'F' , 30 ),
                         ( 'S' , 10 ),
                         ( 'Cl', 10 ) ) )
    types = typef.keys()
    # make a list of typef with max no of combination of atom1 with atom2 
    def trianglen(typef,key1,key2):
        if key1==key2: return int( .5 * typef[key1] * ( typef[key1] - 1 ) )
        else:          return typef[key1]*typef[key2]
    ncombis = OrderedDict( ( (''.join(sorted((key1,key2),key = lambda x: typef.keys().index(x)) ), trianglen(typef,key1, key2) ) for key1 in
        types for key2 in types ) )
    ntypes= len(types)
    monos = OrderedDict(( (key,[]) for key in types ))
    # the next line makes dicts of every possible atom combination with combined keys. combined in order as in typef!
    duos  = OrderedDict( ( (''.join(sorted((key1,key2),key = lambda x: types.index(x)) ), [] ) for key1 in types for key2 in types ) )
    if d:
        print "monos, duos:", monos, duos
        print "typef, ncombis", typef, ncombis
    for i in xrange(l): # for every atom
        for j in range(i,l): #so for every combination with that atom not yet visited
            if i==j:
                nuclear = 0.5*xyz[i][2]**(2.4)
                monos[xyz[i][0]].append(nuclear)
            else:
                t = xyz[i][2] * xyz[j][2]
                n = np.sqrt(
                        np.sum(
                            np.square(
                                xyz[i][1] - xyz[j][1] ) ) )
                force = t/n
                duo_indices = ( xyz[i][0], xyz[j][0] )
                duo_key = ''.join(sorted(duo_indices ,key=lambda x: types.index(x) ))
                duos[duo_key].append(force)
    # now sorted every item in the dictionaries and pad with zeros
    if d:
        print "monos, duos:", monos, duos
    for dictio, ntypes in ( (monos, typef  ),
                            (duos , ncombis)):
        for key,value in dictio.iteritems():
            N = ntypes[key]
            dictio[key] = sorted(dictio[key])[::-1] + [0.0] * ( N - len(value) )
    if d:
        print "monos, duos:", monos, duos
    # merge everything together orderly
    monos_flat = [ x for v in monos.itervalues() for x in v ]
    duos_flat  = [ x for v in duos.itervalues()  for x in v ]
    Bag = monos_flat + duos_flat
    print "B",
    if d:
        print "Bag:", Bag
        print "len(Bag):", len(Bag)
        raise SystemExit('stop')
    return Bag


def coulomb(xyz, ctype='norm4'):
    ''' make a coulomb matrix '''
    l = len(xyz)
    C = np.zeros([l,l])
    for i in xrange(l):
        for j in range(l): #changed this from i+1 to j
            if i==j:
                C[i][i]= 0.5*xyz[i][2]**(2.4)
            else:
                t = xyz[i][2] * xyz[j][2]
                n = np.sqrt(
                        np.sum(
                            np.square(
                                xyz[i][1] - xyz[j][1] ) ) )
                C[i][j] = t/n
                # let not make it symmetric. because we don't use these elements                  
    print "&",
    if ctype=='norm1': #return a sorted Coulomb matrix based on norm
        return symsort(C)
    elif ctype=='norm2':
        return np.sort(np.linalg.norm(C,axis=0))[::-1]
    elif ctype=='norm3':
        #make a sorted Coulomb matrix
        C = symsort(C)
        if True:
            maxn = 56   # max no of atoms for the adamantane and diamantane derivatives. 
            C = padzeros(C,maxn=maxn)
            return C[ np.tril_indices(maxn) ]
        else:
            # return a lower triangular matrix of the coulomb matrix. 
            return C[ np.tril_indices(l) ]
    elif ctype=='norm4':
        # here no triangularization
        C = symsort(C)
        C = padzeros(C, maxn=56)
        #plotmat(C,log=False)
        return C.flatten()
    else:
        assert ctype=='normal'
        return C

def symsort(mat):
    indexlist = np.argsort(np.linalg.norm(mat,axis=1))[::-1]
    return mat[indexlist][:,indexlist]

def padzeros(M, maxn=76 ):
    ''' pads the matrix M until size is size*size '''
    npad = maxn - M.shape[0]
    padM = np.pad(M, (0, npad), 'constant', constant_values=0.0)
    return padM
