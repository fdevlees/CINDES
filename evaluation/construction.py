""" subfunctions for inputfile construction """
import pprint
from pprint import pformat
from itertools import izip, islice
from re import findall
import re
import logging
from CINDES.utils.writings import log_io
#logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

#
debug = False


pp = pprint.PrettyPrinter(indent=4, width=100)


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))


def indtocon(index):
    conf = []
    for item in index.split('_'):
        splitted = re.findall(r"[a-zA-Z]+|\d+", item)
        site = findall('[A-Z0-9][^A-Z1-9]*', splitted[0])
        if len(splitted) == 2:
            dihedral = splitted[1]
            site.append(dihedral)
        conf.append(site)
    return conf

def contoind(conf):
    # return '_'.join([''.join(item) for item in conf])
    # without dihedrals:
    return '_'.join([''.join(filter(lambda x:str(x).isalpha(), item)) for item in conf])


def contoint(conf, array):
    ''' makes an integer list representation of conf '''
    inconf = [site.index(group) for site, group in zip(array, conf)]
    return inconf


def intocon(inconf, array):
    ''' transforms integer list back to normal conf representation '''
    conf = [site[index] for site, index in zip(array, inconf)]
    return conf


def demethyl(passive, defaultgroups):
    """here is now a quite simple operations but i here have
    an open option to fix some other groups later on
    now for each site the methyl group is changed for an H.
    """
    newpassive = []

    # default groups for sites that will not be -H
    if defaultgroups is None:
        defaultgroups = {}

    for i, passivesite in enumerate(passive):
        #coreindex = int(passivesite[0][1])
        try:
            siteindex = int(passivesite[1][1])
        except IndexError:
            print passivesite

        if siteindex in defaultgroups:
            defaultconf = indtocon(defaultgroups[siteindex])[0]
            ldefconf = len(defaultconf)
            if ldefconf==1:
                # no substitution at all even no H
                continue
            else:
                passivesite = passivesite[:ldefconf-1]
                for j, atom in enumerate(defaultconf[1:]):
                    passivesite[j][0]=atom
                newpassive.append(passivesite)
            print "final passivesite:", passivesite
        else:
            # remove the hydrogens from the methyl groups
            passivesite = [passivesite[0]]
            # change carbons to hydrogens
            passivesite[0][0] = 'H'
            newpassive.append(passivesite)
    return newpassive


def hydrogenizer(totalmat):
    ''' This function is meant to be a kind of geom optimizer.
    All standard bond lengths belong to a C-C bond therefore:
    Each -H bond is moved to an average C-H bond length
    Each -N bond except for the core, is decreased
    Each -O bond except for the core, is decreased
    '''
    for i in range(len(totalmat)):
        if totalmat[i][0] == 'H' and 1.6 > float(totalmat[i][2]) > 1.2:
            totalmat[i][2] = '1.1'
        elif totalmat[i][0] == 'O' and i > 9 and (float(totalmat[i][2]) > 1.5 or float(totalmat[i][2]) < 1.0):
            totalmat[i][2] = '1.35'
        elif totalmat[i][0] == 'N' and i > 9 and float(totalmat[i][2]) < 2.0 and (float(totalmat[i][2]) > 1.5 or float(totalmat[i][2]) < 1.2):
            totalmat[i][2] = '1.4'
    return totalmat

# @profile


def matrixmerger2(core, active, passive):
    ''' this has to be translated to numpy '''
    import numpy as np
    V = np.vstack  # is a tweaked form of concatenate
    C = np.concatenate
    try:
        actpas = C((V(active), V(passive)))
    except (UnboundLocalError, ValueError):
        if active==[] and passive==[]:
            actpas = []
        elif passive==[]:
            actpas = V(active)
        elif active == []:
            actpas = V(passive)
        else:
            print "passive:", passive
            print "active:", active
            raise

    core.extend(actpas)

    #totalmat = []
    # for i in range(len(core)):
    #    totalmat.append(core[i])
    # for j in range(len(active)):
    #    # THIS PART IS TIME CONSUMING (all 3 next lines)
    #    for l in range(len(active[j])):
    #        totalmat.append(active[j][l])
    # for k in range(len(passive)):
    #    totalmat.append(passive[k][0])
    # with open('TOTALMAT','w') as tmfid:
    #    tmfid.write(pprint.pformat(totalmat))
    # return totalmat
    return core


def get_configurations(startconf, array, k, run=None):
    'select on site k all the configurations with the different functionalizations for that site present in array'
    configurations = [startconf[0:k] + [array[k][i]] + startconf[k + 1:] for i in range(len(array[k]))]
    return configurations

def get_molecules(startconf, array, k, run):
    from CINDES.molecule import Molecule, Population
    confs = get_configurations(startconf, array, k, run)
    individuals = [Molecule(conf=conf) for conf in confs]  # list of molecules
    return individuals

def get_molecules_SD(startconf, array, run):
    ''' same as get_molecules but now every single mutation '''
    confs = []
    for i in run.restingsites:
        confs.extend(get_configurations(startconf, array, i, run=run))

    # remove duplicates by sorting and subsequently only adding when the previous one is not similar
    sortedconfs = sorted(confs)
    confs = [sortedconfs[i] for i in xrange(len(sortedconfs)) if i == 0 or sortedconfs[i] != sortedconfs[i - 1]]

    individuals = [Molecule(conf=conf) for conf in confs]  # list of molecules
    mols_todo, mols_nodo = check_in_table(individuals, table, run.props)
    return mols_todo, mols_nodo


def check_in_table(individuals, table, props=set(), check_ignored=False):
    ''' this function checks if molecules are already present in the table
    if not it is optionally checked if molecules are already flagged as discarded or
    if the molecule is already ignored earlier, i.e., the molecule name is present in the IGNORED file
    '''
    mols_todo = individuals[:]
    mols_nodo = []
    neglect_dihedrals = True
    i = 0
    digits = '0123456789'
    if table:
        for individual in individuals:
            if individual.index in table:
                # check if the right properties are given for this property
                if all(prop in table[individual.index] for prop in props):
                    individual.props = table[individual.index]
                    individual.predicted = False
                    # individual.Pvalue=0.0 #function!
                    #raise SystemExit('implement furter')
                    mols_todo.remove(individual)
                    mols_nodo.append(individual)
                    i += 1
                    #print "already calculated:", individual.index
        if i:
            logging.info("{} molecules are already in database".format(i))
    if debug:
        print "mols_todo:", mols_todo, "mols_nodo:", mols_nodo

    if check_ignored:
        try:
            # open the file 
            with open('IGNORED', 'r') as f:
                ignored = set(f.read().split('\n'))

            # check for every mol if:
                # - already in IGNORED file
                # - already flagged as discarded
            for mol in mols_todo[:]:
                if mol.index in ignored:
                    print mol, 'in IGNORE'
                    mol.ignore = True
                    mol.IsDiscarded = True
                    mols_nodo.append(mol)
                    mols_todo.remove(mol)
                elif mol.IsDiscarded is True:
                    mols_nodo.append(mol)
                    mols_todo.remove(mol)
        except IOError:
            pass

    return mols_todo, mols_nodo  # indicesfull are all the indices.


def constructor2(conf, core, active, passive, links=None, defaultgroups=None):
    ''' This is the main constructor of the zmatrix for a given configuration using
    the core, active and passive zmatrices. Also symmetry links can be given

        >demethyl
        >doper2
        >substituter2
        >matrixmerger2
        >hydrogenizer
    '''
    def extend_conf(conf, links):
        #print "Symmetry applied!",
        # allocate room for new sites:
        conf.extend([[] for _ in range(len(links))])
        # for all links: fill site with same group as the linked site.
        for i, j in links:
            conf[j - 1] = conf[i - 1]
        return conf

    # 1. if links change conf to extended conf
    if links:
        conf = extend_conf(conf, links)

    #print "CONFIGURATION:",conf
    #if core is nitrogen passive = total remove
    passive = demethyl(passive, defaultgroups)

    # set counter for nth atom in new zmat
    count = 0
    # first the core part
    count += len(core)

    # start with the first site in active
    indices_empty_active = []
    for i in range(len(conf)):
        # read the first element of conf
        # if (not conf[i][0]=='C' or conf[i] == ['C','O']):
        if True:
            if debug:
                print "in constructor 2. before doper:", conf[i], "len passive:", len(passive)
            core, passive = doper2(conf[i], active[i], core, passive)
            if debug:
                print "in constructor 2. AFTER  doper:", conf[i], "len passive:", len(passive)
        # now for EACH! one goes to the substituter
        active[i], count = substituter2(conf[i], active[i], count)
        if not active[i]:
            indices_empty_active.append(i)
    # remove empty elemements from active:
    for i in indices_empty_active[::-1]:
        del active[i]
    mat = matrixmerger2(core, active, passive)
    mat = hydrogenizer(mat)
    return mat


def doper2(group, geom, core, passive):
    debug = False
    if debug:
        print "group:", group
        print "geom:", geom
        print "core"
        for item in core:
            print item
        print "passive:", passive
    # the actual doping command. taking care of index difference Gaussian/Python
    coreindex = int(geom[0][1])
    if not group[0] == 'C':
        # this assumes standard a C is present.
        # this also prevents other groups from overwriting dopants on this site
        # only one site is allowed to have the dopants in that case btw!
        core[coreindex - 1][0] = group[0]
    if debug:
        print "coreindex:", coreindex, "group:", group
        print "core after doping:",
        for item in core:
            print item
    #print "coreindex is: ", coreindex
    if group in [['O'], ['S'], ['C', 'O'], 'O', 'S', 'CO']:
        # we remove the hydrogen at the passive site on that location
        for i in range(len(passive)):
            # number [0][1] is the former C index. it bond length index is the index of the number
            # in the core where it is attached to.
            # het is om het even of: 1. int(str)==int or 2. str(int)==int
            if passive[i][0][1] == str(coreindex):
                del passive[i]
                break
    return core, passive


def geomfiller(zma, geom, count):
    ''' this function fills the zma of a functionalisation into the -methyl geometry of that site '''
    #print "geom1:", geom
    nagroup = len(zma)
    nageom = len(geom)
    delta = nagroup - nageom
    if delta > 0:  # make geom of same length as zma
        for _ in range(delta):
            geom.insert(-1, geom[-1][:])  # now len(geom) is 5
    # the first entry of geom is always the carbon atom.
    # its bond length and type can change. angle and dihedral stay fixed
    # type = geom[0][0] is zma[0][0]
    geom[0][0] = zma[0][0]
    if len(zma[0]) >= 2:  # test if uberhaupt a bond length is given. if so
        geom[0][2] = zma[0][2]  # replace bond length
    # from the first atom in geom i take its bond index and angle index. dihedral index not needed
    bi = geom[0][1]  # this are strings
    ai = geom[0][3]  # string
    # afspraak1: if in zma an index is 0 it will become the angle index
    # afspraak2: if in zma an index is 1 it will become the bond index
    # afspraak3: if in zma an index is empty?
    # for all the other entries

    # THIS PART IS TIME CONSUMING:
    #print "zma:", zma
    #print "geom:", geom
    #print "count:", count, "bi", bi, "ai", ai

    for i in range(1, len(zma)):  # loop over zma except first entry
        # test entries of i and if they exist, fill geom with the right thing
        for j in [0, 2, 4, 6]:  # just replacements of strings
            geom[i][j] = zma[i][j]
        for j in [1, 3, 5]:  # the indexjes
            if zma[i][j] == 0:
                geom[i][j] = ai
            elif zma[i][j] == 1:
                geom[i][j] = bi
            else:
                geom[i][j] = str(zma[i][j] + count - 1)
    #print "final geom:", geom
    #raise SystemExit('stop')
    return geom


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#----- SUBSTITUTER2 IS WORKING NORMALLY ---#


def substituter2(group, geom0, count):
    # aantal te verwijderen H is gerelateerd aan de lengte
    # 5 - 0 / 4 - 1 etc
    # ndelh = 5 - length(group)
    from copy import copy
    geom = copy(geom0)
    dihedral = None
    if is_float(group[-1]):
        dihedral = float(group[-1])
        group = group[:-1]
    if len(group) == 7:
        # C N H C H H H
        # C S O C H H H
        if group == ['C', 'N', 'H', 'C', 'H', 'H', 'H']:
            zma = [['N', 1, '1.457'],
                   ['C', 2, '1.457', 1, '140.5', 0, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '60.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '300.0'],
                   ['H', 2, '1.018', 3, '109.5', 4, '-60.1']]
        elif group == ['C', 'S', 'O', 'C', 'H', 'H', 'H']:
            zma = [['S', 1, '1.457'],
                   ['C', 2, '1.457', 1, '140.5', 0, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '60.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '300.0'],
                   ['O', 2, '1.500', 3, '109.5', 4, '-60.1']]
        else:
            print "group:", group
            raise SystemExit('group not recognized')
        if dihedral:
            zma[1][6] = dihedral
        geom = geomfiller(zma, geom, count)
        count += len(zma)
    elif len(group) == 6:  # this is for the -OCH3 or -SCH3 group for example
        # i have to introduce a new line in geom.

        if group in [['C', 'O', 'C', 'H', 'H', 'H'], ['C', 'S', 'C', 'H', 'H', 'H']]:
            geom.insert(-1, geom[-1][:])  # now len(geom) is 5
            zma = [['O', 1, '1.457'],
                   ['C', 2, '1.457', 1, '140.5', 0, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '60.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '180.0'],
                   ['H', 3, '1.095', 2, '109.5', 1, '300.0']]
            if group == ['C', 'S', 'C', 'H', 'H', 'H']:
                zma[0][0] = 'S'
            # change C to O or S
            if group[1] == 'S':
                geom[0][0] = 'S'
            elif group[1] == 'O':
                geom[0][0] = zma[0][0]
            # change one H to C
            if True:
                if dihedral:
                    zma[1][6] = dihedral
                geom = geomfiller(zma, geom, count)
            else:
                geom[1][0] = zma[1][0]
                geom[1][1] = str(zma[1][1] + count)  # count + 1 index of the N
                geom[1][2] = zma[1][2]  # bond length
                geom[1][4] = zma[1][4]  # angle. (geom[2][3] stays the same) as do 5 and 6
                # now the 3 H atoms
                for i in [2, 3, 4]:
                    geom[i][0] = zma[i][0]
                    geom[i][1] = str(zma[i][1] + count)
                    geom[i][2] = zma[i][2]  # bond length
                    geom[i][3] = str(zma[i][3] + count)
                    geom[i][4] = zma[i][4]
                    geom[i][6] = zma[i][6]
        elif group == ['C', 'S', 'O', 'O', 'O', 'H']:
            #print "sulfonyl oid group"
            zma = [['S', 1, '1.79'],
                   ['O', 2, '1.46', 1, '110.0', 0, '51.0'],
                   ['O', 2, '1.46', 1, '110.0', 3, '134.7'],
                   ['O', 2, '1.65', 3, '108.8', 4, '-126.25'],
                   ['H', 5, '0.97', 2, '107.4', 3, '66.3']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
        count += 5

    elif len(group) == 5:
        # it is CCHHH or CCFFF or CCOOH
        geom[0][0] = group[1]  # until now this stays just C
        geom[1][0] = group[2]
        geom[2][0] = group[3]
        geom[3][0] = group[4]
        if group in [['C', 'C', 'F', 'F', 'F'], ['C', 'C', 'H', 'H', 'H'],
                     ['N', 'C', 'H', 'H', 'H'], ['C', 'N', 'H', 'H', 'H']]:
            geom[1][1] = str(count + 1)
            geom[2][1] = str(count + 1)  # if also attached to that one
            geom[3][1] = str(count + 1)
            if group in [['C', 'N', 'H', 'H', 'H']]:
                print "WARNING: charged group: +1"
        elif group in [['C', 'C', 'O', 'O', 'H']]:
            #print "carbonic acid"
            geom[1][1] = str(count + 1)  # this is the =O
            geom[2][1] = str(count + 1)  # this is the -O-

            geom[1][2] = 1.3  # bond length C=O
            geom[1][4] = 120.1  # angle coreC-C=O
            if dihedral:
                geom[1][6] = dihedral
            else:
                geom[1][6] = 0.1  # dihedral with one of the core
            geom[2][2] = 1.3  # bond length C-O
            geom[2][4] = 120.1  # angle coreC-C=O
            if dihedral:
                geom[2][6] = dihedral - 180.0
            else:
                geom[2][6] = 180.1  # dihedral with one of the core

            geom[3][1] = str(count + 3)  # attached to -O-
            geom[3][3] = str(count + 1)  # angled with -C
            geom[3][5] = str(count + 2)  # dihedraled with =O

            geom[3][2] = 0.9
            geom[3][4] = 109.5
            geom[3][6] = 2.1  # just to make it changebla i don't make it zero
        else:
            print group
            raise Exception('Group does not exist')

        count += 4
    # when length == 4 then it is an amine or nitro group
    elif len(group) == 4:
        # need to remove one of the hydrogens
        del geom[3]
        geom[0][0] = group[1]
        geom[1][0] = group[2]
        geom[2][0] = group[3]

        geom[1][1] = str(count + 1)
        geom[2][1] = str(count + 1)  # if also attached to that one

        if group == ['C', 'C', 'C', 'H']:
            # zma = [['C',1, '1.4554490', 0, '111.0852500', -1, '121.9277283'],
            #       ['C',1, '2.6636450', 0, '111.0852500', -1, '121.9277283'],
            #       ['H',1, '3.7296390', 0, '111.0852500', -1, '121.9277283']]
            # geom of second ethyn carbon is same as first:
            geom[1] = geom[0][:]
            # also the H atom is same as first
            geom[2] = geom[0][:]
            # r of C
            geom[1][2] = '2.6636'
            # r of H
            geom[2][2] = '3.7296'

            # geom[0][0] = group[1] # change C to O or also C
            # geom[1][0] = group[2] # change H to N or also H
            geom[2][0] = group[3]  # change H to N or also H
        elif group == ['C', 'N', 'O', 'O']:
            zma = [['N', 1, '1.51'],
                   ['O', 2, '1.227', 1, '117.01', 0, '30.1'],
                   ['O', 2, '1.227', 1, '117.01', 3, '178.5']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
        elif group == ['C', 'C', 'O', 'O']:
            print "WARNING: charged group: -1"
            zma = [['C', 1, '1.51'],
                   ['O', 2, '1.227', 1, '117.01', 0, '30.1'],
                   ['O', 2, '1.227', 1, '117.01', 3, '178.5']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
        elif group == ['C', 'C', 'H', 'H']:
            zma = [['C', 1, '1.51'],
                   ['H', 2, '0.9', 1, '117.01', 0, '30.1'],
                   ['H', 2, '0.9', 1, '117.01', 3, '178.5']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
        elif group in [['C', 'O', 'O', 'H'], ['C', 'S', 'O', 'H']]:
            zma = [['O', 1, '1.4'],
                   ['O', 2, '1.45', 1, '104.0', 0, '-61.0'],
                   ['H', 3, '0.96', 2, '98.5', 1, '170.0']]
            if group == ['C', 'S', 'O', 'H']:
                zma[0][0] = 'S'
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
        elif group == ['C', 'N', 'H', 'H']:
            geom[0][2] = '1.465'
            geom[1][2] = '1.02'
            geom[2][2] = '1.02'
            if dihedral:
                geom[1][6] = dihedral
                geom[2][6] = dihedral - 120.0
        elif group in [['C', 'C', 'H', 'O'], ['C', 'C', 'O', 'H']]:
            geom[1][1] = str(count + 1)  # this is the =O
            geom[2][1] = str(count + 1)  # this is the -H

            geom[1][2] = 1.18  # bond length C=O
            geom[1][4] = 125.1  # angle coreC-C=O
            geom[1][6] = 0.1  # dihedral with one of the core

            geom[2][2] = 1.11  # bond length C-O
            geom[2][4] = 115.1  # angle coreC-C=O
            geom[2][6] = 180.1  # dihedral with one of core

        count += 3  # three atoms added to activemat
    elif len(group) == 3:
        # when length is 3 it is a COH, NOH? or CCN group
        # we need to remove two hydrogens
        del geom[3]
        del geom[2]
        if group in [['C', 'N', 'C'], ['C', 'C', 'N']]:
            #print "cyano"
            geom[1] = geom[0][:]
            geom[1][2] = '2.62'
        else:
            # because the N now directly binds to the core
            # geom[1][1] verwijst naar the position of geom[0]
            geom[1][1] = str(count + 1)  # the bondlength index has to count+1 denk ik
        geom[0][0] = group[1]  # change C to O or also C
        geom[1][0] = group[2]  # change H to N or also H
        count += 2  # because two atoms are added to the activemat
        # shorten the C-O bond a bit
        if group == ['C', 'O', 'H']:
            geom[0][2] = '1.42'
            geom[1][2] = '0.97'
        if dihedral:
            geom[1][6] = dihedral
    elif len(group) == 2:
        # it is an C-H,C-F,Si-H or C-Cl group. # or it is C Ph
        if group == ['C', 'Ph']:
            #print "benzene group"
            zma = [['C', 1, '1.51'],
                   ['C', 2, '1.401', 1, '120.1', 0, '30.1'],
                   ['C', 3, '1.395', 2, '120.9', 1, '178.5'],
                   ['C', 4, '1.395', 3, '120.1', 2, '0.1'],
                   ['C', 5, '1.395', 4, '119.5', 3, '0.1'],
                   ['C', 6, '1.395', 5, '120.1', 4, '0.1'],
                   ['H', 3, '1.087', 2, '119.3', 7, '180.1'],
                   ['H', 4, '1.086', 3, '119.8', 2, '180.1'],
                   ['H', 5, '1.086', 4, '119.9', 3, '180.1'],
                   ['H', 6, '1.086', 5, '120.1', 4, '180.1'],
                   ['H', 7, '1.087', 6, '119.7', 5, '180.1']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
            count += len(zma)

        elif group == ['C', 'Th']:
            #print "thiophene group"
            zma = [['C', 1, '1.4952220'],
                   ['C', 2, '1.3695029', 1, '128.4633279', 0, '0.0023405'],
                   ['S', 2, '1.7399783', 3, '110.1483490', 1, '179.9964292'],
                   ['C', 3, '1.4258839', 2, '113.6486669', 4, '-0.0147807'],
                   ['H', 3, '1.0851597', 2, '122.4471034', 4, '179.9989798'],
                   ['C', 5, '1.3649242', 3, '112.6274587', 2, '0.0275216'],
                   ['H', 5, '1.0843806', 3, '124.0375540', 2, '-179.9905779'],
                   ['H', 7, '1.0814341', 5, '128.6304156', 3, '179.9872962']]
            if dihedral:
                zma[1][6] = dihedral
            geom = geomfiller(zma, geom, count)
            count += len(zma)

        else:
            # all the three hydrogens need to be removed

            #print "len 2!"
            del geom[3]
            del geom[2]
            del geom[1]
            # and thange the C to H or Cl or F
            geom[0][0] = group[1]
            if group == ['C', 'O']:  # changes angle. analyses the different cases
                #print "keton"
                geom[0][2] = '1.215'
                geom[0][4] = '120.5'
                if -64.0 <= float(geom[0][6]) <= -60.0 or -170.0 <= float(geom[0][6]) <= -164.0:
                    geom[0][6] = '-117.1'
                if 62.0 <= float(geom[0][6]) <= 65.0 or 170.0 <= float(geom[0][6]) <= 178.0:
                    geom[0][6] = '122.0'
            count += 1
            # the first atom is directly attached to the core so the zmat indices don't change
    elif len(group) == 1:
        # when the length of the group is one, that means it is doped with N,O,S,B,P
        del geom[3]
        del geom[2]
        del geom[1]
        del geom[0]
        return None, count
    # so whole geom is deleted actualy
    # there will not appear any of this ones in the activemat so count is not changed
    return geom, count


if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    multcharge = re.compile('^\-?[0-9]\s[0-9]')
    with open(filename) as fid:
        for line in fid:
            if multcharge.match(line):
                print "match!"
                print multcharge.match(line).group()
                zmat = []
                line = next(fid)
                while not line == '\n':
                    zmat.append(line.split())
                    line = next(fid)
    ncore = 12
    for pos in positions:
        hornot = maker(zmat, pos, ncore)
        if not hornot == 1:
            maker2(zmat, pos, ncore)
        del zmatnew
        del ind
        del hline
        del item
    print "DONE"
