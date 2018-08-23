""" reader of zmat """
debug = False

import pprint
import logging
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
# compu chem. library


# NEW IDEA:
class ZMatrix(object):
    def __init__(self, zmatfile='ZMAT', **kwargs):
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
        (coremat, activemat, passivemat) = sitesplitter(zmat, param['ncore'], param['sites'], param['nch3'])
        return


def geometry(zmatfile='ZMAT', **param):
    '''reads the zmat from a file and splits it'''
    # note that zmatfile is now in **param
    zmat, fileid = zmatread(zmatfile)
    zmatdic = zmatvalues(fileid)
    logging.debug(pprint.pformat(zmatdic))
    fileid.close()
    # FORMATTING AND SPLITTING OF ZMATRIX
    zmat = zmatprinter(zmat, zmatdic)
    logging.debug("zmat:\n" + pprint.pformat(zmat))
    (coremat, activemat, passivemat) = sitesplitter(zmat, param['ncore'], param['sites'], param['nch3'])
    # now i save here the matrices for later use, and then the others are allowed to change for each molecule
    logging.info('activemat:' + pprint.pformat(activemat))
    logging.info('passivemat:' + pprint.pformat(passivemat))
    logging.info("----- END FORMATTING & SPLITTING -----")
    Total_Zmat = {'core': coremat, 'active': activemat, 'passive': passivemat}
    return Total_Zmat


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
        lijntje = line.split()
        dictio[lijntje[0]] = lijntje[1]
    return dictio


def zmatprinter(czmat, dictio):
    """
    The aim of this function is to make a more compact zmatrix
    printing the values directly by replacing the variable names

    I think i just have to print zmat, but checking for each element
    if it is a value from the dictionary. If so print value instead of key

    maybe, i need to check only the values
    [1][2]
    [2][2] & [2][4]
    [3][2] & [3][4] & [3][6]
    same for the others
    """
    # w = write (will erase existing file) a = append (to the end of file)
    # you can only write strings to a file
    # use join to mute the objects in a list of strings together. separated
    # from each other by a whitespace.
    # print "czmat[1][2]:", czmat[1][2]
    # search the dictionary
    for i in range(1, len(czmat) - 1):  # first row already printed, contains no value
        for key in dictio.keys():
            for j in [2, 4, 6]:  # 2 bondlengt #4 angle #6 dihedral
                # now the first two don't have 4 and 6 giving IndexError so
                try:
                    # look if key matches with object
                    if key == czmat[i][j]:
                        # actually i prefer not to change it, only to print it
                        czmat[i][j] = dictio[key]
                except IndexError:  # that is the name of the Error
                    pass  # go on to the next one
    return czmat


def sitesplitter(zmatrix, natomscore, index, nsites):  # here nsites is number of possible sites so nch3
    """this module splits zmatrix in core part and part for sites
    next the sites are split up in an active and passive part
    """
    coremat = zmatrix[0:natomscore]
    # pprint.pprint(coremat)
    sitemat = []
    # here i split the rest, that are al ch3 groups
    for i in range(nsites):
        sitemat.append(zmatrix[natomscore + 4 * i:natomscore + 4 * (i + 1)])
    logging.debug('sitemat' + pprint.pformat(sitemat))
    # now split in passive and active part
    activeindex = [((x - natomscore + 3) / 4) - 1 for x in index]
    logging.info('activeindex: ' + str(activeindex))
    logging.debug('len(sitemat)' + str(len(sitemat)))
    active = []
    passive = []
    # here mistake!!! and inefficient better to loop over activesites
    for i in activeindex:
        active.append(sitemat[i])
    for i in range(len(sitemat)):
        if not i in activeindex:
            passive.append(sitemat[i])
    return (coremat, active, passive)


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
    print "activematrix: "
    pprint.pprint(activematrix)
    print "passivematrix: "
    pprint.pprint(passivematrix)
    # now we will change the passivematrix such that it are only single
    # hydrogen atoms and no methyl groups
