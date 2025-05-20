""" reader of zmat """
debug = False

import pprint
import logging
#logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
# compu chem. library


# NEW IDEA:
class ZMatrix(object):
    def __init__(self, zmatfile, **kwargs):
        self.zmatfile = zmatfile
        self.read()
        return

    def read(self):
        self.zmat, fileid = self.zmatread()
        self.zmatdic = self.zmatvalues(fileid)
        fileid.close()
        return

    def split(self):
        # FORMATTING AND SPLITTING OF ZMATRIX
        self.zmat = zmatprinter(self.zmat, self.zmatdic)
        logging.debug("zmat:\n" + pprint.pformat(zmat))
        (coremat, activemat, passivemat) = sitesplitter(zmat, param['ncore'], param['sites'])
        return

def get_cartesian(xyzfile='XYZ', **param):
    with open(xyzfile, 'r') as f:
        xyz_raw = [ line.split() for line in f ]
    xyz = [ [line[0], list(map(float, line[1:]))] for line in xyz_raw ]
    print('xyz:', xyz)
    return xyz


def geometry(zmatfile, **param):
    '''reads the zmat from a file and splits it'''
    # note that zmatfile is now in **param
    zmat, fileid = zmatread(zmatfile)
    zmatdic = zmatvalues(fileid)
    logging.debug(pprint.pformat(zmatdic))
    fileid.close()
    # FORMATTING AND SPLITTING OF ZMATRIX
    zmat = fill_zmat(zmat, zmatdic)
    logging.debug("zmat:\n" + pprint.pformat(zmat))
    coremat, activemat, passivemat, corresp = sitesplitter(zmat, param['ncore'], param['sites'])
    # now i save here the matrices for later use, and then the others are allowed to change for each molecule
    logging.info("----- END FORMATTING & SPLITTING -----")
    Total_Zmat = {'core': coremat, 'active': activemat, 'passive': passivemat}
    return Total_Zmat, corresp


def zmatread(filename):
    """
    This function reads the Z matrix
    """
    fid = open(filename, 'r')
    # create a variable list where to put the zmatrix in
    czmat = []
    while True:
        line = fid.readline().split()
        # print line,
        czmat.append(line)
        # stop when an empty line is encountered
        if line == '\n':
            break
        if not line:
            break
    return (czmat, fid)


def zmatvalues(fileid):
    """
    This function will make a dictionary of all the zmat values
    NOTE: this function read until end of file.
          A leading empty line is not allowed
    """
    dictio = {}
    for line in fileid:
        splitted = line.split()
        dictio[splitted[0]] = splitted[1]
    return dictio


def fill_zmat(czmat, dictio):
    """
    filling in the values by replacing the variable names B#, A# and D#
    """
    # print "czmat[1][2]:", czmat[1][2]
    # search the dictionary
    for i in range(1, len(czmat) - 1):
        for j in [2, 4, 6]:  # 2 bondlengt #4 angle #6 dihedral
            # now the first two don't have 4 and 6 giving IndexError so
            try:
                # look if key matches with object
                czmat[i][j] = dictio[czmat[i][j]]
            except IndexError:
                pass
    return czmat


def sitesplitter(zmatrix, natomscore, index):
    """this module splits zmatrix in core part and part for sites
       next the sites are split up in an active and passive part
    """
    i=natomscore
    coremat = zmatrix[0:i]
    # pprint.pprint(coremat)
    # here i split the rest, that are all ch3 groups

    # OLD IMPLEMENTATION
    # sitemat = []
    #for i in range(nsites):
    #    sitemat.append(zmatrix[natomscore + 4 * i:natomscore + 4 * (i + 1)])
    #logging.debug('sitemat' + pprint.pformat(sitemat))
    # now split in passive and active part
    #activeindex = [((x - natomscore + 3) / 4) - 1 for x in index]
    #logging.debug('len(sitemat)' + str(len(sitemat)))
    #active = []
    #passive = []
    # here mistake!!! and inefficient better to loop over activesites
    #for i in activeindex:
    #    active.append(sitemat[i])
    #for i in range(len(sitemat)):
    #    if not i in activeindex:
    #        passive.append(sitemat[i])
    #return (coremat, active, passive)

    # NEW IMPLEMENTATION
    corresp= dict()
    active = dict()
    passive = []
    line = zmatrix[i]
    while True:
        if line==[]:
            break
        group = [line]
        corresp[i+1]=int(line[1])
        if i+1 in index:
            active[i+1]=group
            #active.append(group)
        else:
            passive.append(group)
        while True:
            i+=1
            line = zmatrix[i]
            if line==[]:
                break
            if line[0]=='H':
                group.append(line)
            else:
                break
    active = [ active[i] for i in index ]

    print("corresp:", corresp)
    #pprint.pprint(active)
    #pprint.pprint(passive)
    return coremat, active, passive, corresp


if __name__ == "__main__":
    import sys
    file = sys.argv[1]
    # from here, this still has to be included in mail file
    (zmat, fileid) = zmatread(file)
    pprint.pprint(zmat)
    zmatdic = zmatvalues(fileid)
    fileid.close()
    pprint.pprint(zmatdic)
    zmatprinter(zmat, zmatdic)
    # zmat now contains the values for bonds/angles/dihedrals
#pprint.pprint( zmat)
    # now i want to split them in:
    #  1 core part, first 10 atoms
    # 16 sites of 4 lines
    # remove the last line
    # hard coding that the core are the first 10 atoms
    activesites = [39, 23, 31, 55, 47]
    (coremat, activematrix, passivematrix) = sitesplitter(zmat, 10, activesites)
    # the index of the carbons of the methyl groups for the active sites
    # later on we can try to derive this matrix from the inputfile
    # now it is hard coded
    # the corresponding index in the sitemat
    # now we make two matrices, one with the active sites
    # and one with the sites that stay fixed
    print("activematrix: ")
    pprint.pprint(activematrix)
    print("passivematrix: ")
    pprint.pprint(passivematrix)
    # now we will change the passivematrix such that it are only single
    # hydrogen atoms and no methyl groups
