""" subfunctions for inputfile construction """
import pprint
from pprint import pformat
from itertools import izip, islice
from re import findall
import re
import logging
from CINDES4.utils.writings import log_io
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)

#
simple = 1
semiempirical = 1



pp = pprint.PrettyPrinter(indent=4, width=100)


def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

def indtocon(index):
    conf = []
    for item in index.split('_'):
        splitted = re.findall(r"[a-zA-Z]+|\d+", item)
        site = findall('[A-Z0-9][^A-Z1-9]*', splitted[0] )
        if len(splitted)==2:
            dihedral=splitted[1]
            site.append(dihedral)
        conf.append(site)
    return conf

    #return  [ findall('[A-Z0-9][^A-Z1-9]*',item) for item in index.split('_') ]
    #return  [ findall('[A-Z0-9][^A-Z1-9]*',item) for item in index.split('_') ]
    #return [list(item) for item in index.split('_')]
    #return  [ findall('[A-Z][^A-Z]*',item) for item in index.split('_') ]

def contoind(conf):
    return '_'.join([''.join(item) for item in conf])

def contoint(conf,array):
    ''' makes an integer list representation of conf '''
    inconf = [ site.index(group) for site,group in zip(array,conf) ]
    return inconf

def intocon(inconf,array):
    ''' transforms integer list back to normal conf representation '''
    conf = [ site[index] for site,index in zip(array,inconf) ]
    return conf

def demethyl(passive):
    """here is now a quite simple operations but i here have
    an open option to fix some other groups later on 
    now for each site the methyl group is changed for an H. 
    """
    for i in range(len(passive)):
        # remove the hydrogens from the methyl groups
        del passive[i][1:4]
        # change carbons to hydrogens
        passive[i][0][0]='H'
    return passive

def hydrogenizer(totalmat):
    ''' This function is meant to be a kind of geom optimizer. 
    All standard bond lengths belong to a C-C bond therefore: 
    Each -H bond is moved to an average C-H bond length
    Each -N bond except for the core, is decreased
    Each -O bond except for the core, is decreased
    '''
    for i in range(len(totalmat)): 
        if totalmat[i][0] == 'H' and float(totalmat[i][2]) > 1.2:
            totalmat[i][2] = '1.1'
        elif totalmat[i][0] == 'O' and i>9 and ( float(totalmat[i][2]) > 1.5 or float(totalmat[i][2]) < 1.0): 
            totalmat[i][2] = '1.35'
        elif totalmat[i][0] == 'N' and i>9 and float(totalmat[i][2]) < 2.0 and (float(totalmat[i][2]) > 1.5 or float(totalmat[i][2]) < 1.2):
            totalmat[i][2] = '1.4'
    return totalmat

def matrixmerger2(core,active,passive):
    totalmat = []
    nact = 0 
    for i in range(len(core)):
        totalmat.append(core[i])
    for j in range(len(active)):
        for l in range(len(active[j])):
            totalmat.append(active[j][l])
            nact += 1
    for k in range(len(passive)):
        totalmat.append(passive[k][0])
    #print "************************************************************************************"
    with open('TOTALMAT','w') as tmfid:
        tmfid.write(pprint.pformat(totalmat))
    return totalmat       

def get_configurations(startconf,array,k, run=[]):
    'select on site k all the configurations with the different functionalizations for that site present in array'
    configurations =  [ startconf[0:k] + [array[k][i]] + startconf[k+1:] for i in range(len(array[k]))]
    if hasattr(run,'nlinks'):
        print "type(run)", type(run)
        for link in run.symlinks:
            (i,j) = (link[0]-1,link[1]-1)
            print "link is:", i, " ",j
            for conf in configurations:
                if not conf[i]==conf[j]:
                    conf[j]=conf[i]
    logging.debug(pprint.pformat(configurations))
    return configurations

def indexmaker3(startconf,array,k,table,run=[]): #for CINDES2.3.py for the new symmetry feature
    '''checks for confs already calculated'''
    print "IN INDEXMAKER3", type(run)
    confs = get_configurations(startconf,array,k,run=run)
    data=[]
    indices = []
    for i in range(len(confs)):
        index = contoind(confs[i])
        indices.append(index)
        #pp.pprint(confs[i])
    indicesfull = indices[:]
    if not table == []:
        for item in table:
            for index,confje in izip(indices[:],confs[:]):
                if item[0] == index:
                    # remove that from the configurations
                    indices.remove(index)
                    confs.remove(confje)
                    # add that item from table to data
                    if item[1]==1:
                        #raise SystemExit('elements in tablebin shouldnt be one')
                        data.append(item)
                    else:
                        new_item = item[:]
                        new_item.insert(1,1)
                        data.append(new_item)
        if not data == []:
            logging.info('filled data with ones already calced:' + pprint.pformat(data))
    return indices,data,confs,indicesfull #indicesfull are all the indices. 

def indexmaker4(table,indices, confs):
    '''checks for confs already calculated this one is used in GA.py'''
    indicesfull = indices[:]
    data = []
    if not table == []:
        for item in table:
            for index,confje in izip(indices[:],confs[:]):
                if item[0] == index:
                    # remove that from the configurations
                    indices.remove(index)
                    confs.remove(confje)
                    # add that item from table to data
                    if item[1]==1:
                        raise SystemExit('elements in tablebin shouldnt be one')
                        data.append(item)
                    else:
                        new_item = item[:]
                        new_item.insert(1,1)
                        data.append(new_item)
        if not data == []:
            logging.info('filled data with ones already calced:' + pprint.pformat(data))
    return indices,data,confs #indicesfull are all the indices. 

def constructor2(conf,core,active,passive):
    '''another constructor now with a counter
    I hope it needs less functions but. yes
        >demethyl
    '''
    #print "CONFIGURATION:",conf
    passive = demethyl(passive)
    count = 0
    #first the core part
    count += len(core)
    #then the active part
    # len(conf)==len(active)
    # start with the first site in active
    for i in range(len(conf)):
        #read the first element of conf
        if (not conf[i][0]=='C' or conf[i] == ['C','O']):
            core,passive=doper2(conf[i],active[i],core,passive)
        # now for EACH! one goes to the substituter
        active[i],count = substituter2(conf[i],active[i],count)
        logging.debug('active' + str(i))
        logging.debug(pprint.pformat(active[i]))
        #here we have to dope and we need to know the correct position in the core
    mat = matrixmerger2(core,active,passive)      
    mat = hydrogenizer(mat) 
    return mat

def doper2(group, geom, core, passive):
    # the actual doping command. taking care of index difference Gaussian/Python
    coreindex = int(geom[0][1])
    core[coreindex-1][0]=group[0]
    #print "coreindex is: ", coreindex
    if group in [['O'],['S'],['C','O']]:
        # we remove the hydrogen at the passive site on that location
        for i in range(len(passive)):
            # number [0][1] is the former C index. it bond length index is the index of the number
            # in the core where it is attached to. 
            # het is om het even of: 1. int(str)==int or 2. str(int)==int
            if passive[i][0][1] == str(coreindex):
                del passive[i]
                break
    return core,passive


def geomfiller(zma,geom,count):
    logging.debug( "geomfiller")
    nagroup = len(zma)
    nageom = len(geom)
    delta = nagroup - nageom
    if delta > 0: # make geom of same length as zma
        for _ in range(delta):
            geom.insert(-1,geom[-1][:]) #now len(geom) is 5
    # the first entry of geom is always the carbon atom.
    # its bond length and type can change. angle and dihedral stay fixed
    # type = geom[0][0] is zma[0][0]
    geom[0][0] = zma[0][0]
    if len(zma[0]) >= 2: # test if uberhaupt a bond length is given. if so
        geom[0][2]= zma[0][2] # replace bond length
    # from the first atom in geom i take its bond index and angle index. dihedral index not needed
    bi = geom[0][1] #this are strings
    ai = geom[0][3] # string
    # afspraak1: if in zma an index is 0 it will become the angle index
    # afspraak2: if in zma an index is 1 it will become the bond index
    # afspraak3: if in zma an index is empty? 
    # for all the other entries
    for i in range(1,len(zma)): # loop over zma except first entry
        # test entries of i and if they exist, fill geom with the right thing
        logging.debug("geom item")
        logging.debug(pprint.pformat(geom[i]))
        logging.debug("zma item")
        logging.debug(pprint.pformat(zma[i]))
        for j in [0,2,4,6]: # just replacements of strings
            geom[i][j] = zma[i][j]
        for j in [1,3,5]: # the indexjes
            if zma[i][j] == 0:
                geom[i][j] = ai
            elif zma[i][j] == 1:
                geom[i][j] = bi
            else:
                geom[i][j] = str(zma[i][j] + count-1)
    return geom

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

#---- NEW MAKERS ARE INTRODUCED HERE! -------#
#---- END MAKERS ----------#
#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter2(zmat,index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat 
    The name of the file contains the index in the name
    Still a lot of hardcoded information:
        - charge
        - multiplicity
        - memory
        - number of processors
        - calculation procedure:
            > geometry optimization
            > unrestricted DFT - B3LYP functional
            > basisset: 6-31G(d)    '''
    #------------
    # this function uses globals: identify, path, gaussianline, gaussianline2, twojob
    # maybe something like:
    # filewriter(zmat,*args,**kwargs):
    # and then calling it with
    # filewriter(zmat, gaussianline, gaussianline2, path=path, identify=identify, charge=0, mult=1)
    # or with
    # filewriter(zmat, **filedic)
    # with filedic is:
    # filedic = {"charge":0,"mult":1,"identify":identify}
    #------------
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    #fid.write("%chk=" + identify + str(index) + ".chk\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    if not paras['nprocs']==1:
        fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
    fid.write(paras['gaussianline'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    # here the zmat
    # core
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    if paras['polar'] == 1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianlinepolar'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
    if paras['ip'] == 1 or paras['aip']==1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
        fid.write(str(paras['charge']+1) + " " + str(paras['mult']+1) + "\n")
        fid.write("\n")
    if paras['ea'] == 1 or paras['aea']==1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
        fid.write(str(paras['charge']-1) + " " + str(paras['mult']+1) + "\n")
        fid.write("\n")
    print "---- FILE PRINTED SUCCESFULLY -----"
    return
#----- END FILEWRITER ----#
#----- BEGIN FILEWRITER JOB TYPE A BDE-MODEL ----#
def filewriterA(zmat,index,**paras): #paras is short for fileparameters
    paras['charge']=0
    paras['mult']=1
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    # gaussianline now is:
    # "# opt=calcfc ub3lyp/6-31g(d) pop=npa freq'
    fid.write(paras['gaussianline1'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    # here the zmat
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    # JOB2 larger basis for E(A) to calculate energy for E(A)BDE
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nproc=2\n")
    #gausline3 is:
    # '# geom=check guess=read b3p86/6-311+G(d,p)'
    fid.write(paras['gaussianline3'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    fid.write("\n")
    # JOB3 E0 E(A)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    #gausline2 is:
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    fid.write("\n")
    # JOB4 IP E(A+)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    #gausline4 is same as gausline2:
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']+1) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    # JOB5 EA E(A-)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    #gausline5 is same as gausline2
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']-1) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
    fid.close()
    return

def filewriterAH(zmat,pos,index,**paras): #paras is short for fileparameters
    paras['charge']=0
    paras['mult']=1
    filename = paras['identify'] + str(index) + "_" + pos + ".com"
    print "in filewriter AH; filename:",filename
    fid=open(paras['path'] + '/' + index + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + "_" + pos + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    fid.write(paras['gaussianline1'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    # zmat:
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + "_" + pos + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    fid.write(paras['gaussianline3'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
    fid.close()
    return

#----- SUBSTITUTER2 IS WORKING NORMALLY ---#
def substituter2(group,geom0,count):
    # aantal te verwijderen H is gerelateerd aan de lengte
    # 5 - 0 / 4 - 1 etc 
    # ndelh = 5 - length(group)
    geom = geom0
    dihedral = None
    if is_float(group[-1]):
        dihedral = float(group[-1])
        group = group[:-1]
    if len(group) == 7:
        # C N H C H H H
        # C S O C H H H
        zma = [ ['N',1,'1.457'],
                ['C',2,'1.457',1,'140.5',0,'180.0'],
                ['H',3,'1.095',2,'109.5',1, '60.0'],
                ['H',3,'1.095',2,'109.5',1,'180.0'],
                ['H',3,'1.095',2,'109.5',1,'300.0'],
                ['H',2,'1.018',3,'109.5',4,'-60.1'] ]
        if dihedral: zma[1][6]=dihedral
        geom=geomfiller(zma,geom,count)
        count += len(zma)       
    elif len(group) == 6: # this is for the -OCH3 or -SCH3 group for example
        # i have to introduce a new line in geom. 

        if group in [['C','O','C','H','H','H'],['C','S','C','H','H','H']]:
            geom.insert(-1,geom[-1][:]) #now len(geom) is 5
            zma = [         ['O',1,'1.457'],
                            ['C',2,'1.457',1,'140.5',0,'180.0'],
                            ['H',3,'1.095',2,'109.5',1,'60.0'],
                            ['H',3,'1.095',2,'109.5',1,'180.0'],
                            ['H',3,'1.095',2,'109.5',1,'300.0']  ]
            # change C to O or S
            if group[1] == 'S':
                geom[0][0] = 'S'
            elif group[1] == 'O':
                geom[0][0] = zma[0][0]
            # change one H to C
            if True:
                if dihedral: zma[1][6]=dihedral
                geom=geomfiller(zma,geom,count)
            else:
                geom[1][0] = zma[1][0]
                geom[1][1] = str(zma[1][1] + count) # count + 1 index of the N
                geom[1][2] = zma[1][2] # bond length
                geom[1][4] = zma[1][4] # angle. (geom[2][3] stays the same) as do 5 and 6
                # now the 3 H atoms
                for i in [2,3,4]:
                    geom[i][0] = zma[i][0]
                    geom[i][1] = str(zma[i][1] + count)
                    geom[i][2] = zma[i][2] # bond length
                    geom[i][3] = str(zma[i][3] + count)
                    geom[i][4] = zma[i][4]
                    geom[i][6] = zma[i][6]
        elif group == ['C','S','O','O','O','H']:
            #print "sulfonyl oid group"
            zma = [ ['S', 1, '1.79'],
                    ['O', 2, '1.46', 1, '110.0', 0, '51.0'   ],
                    ['O', 2, '1.46', 1, '110.0', 3, '134.7'  ],
                    ['O', 2, '1.65', 3, '108.8', 4, '-126.25'],
                    ['H', 5, '0.97', 2, '107.4', 3, '66.3'   ] ]
            if dihedral: zma[1][6]=dihedral
            geom = geomfiller(zma,geom,count)
        count += 5

    elif len(group) == 5:
        # it is CCHHH or CCFFF or CCOOH
        #geom[0][0] = group[1] #until now this stays just C
        geom[1][0] = group[2]
        geom[2][0] = group[3]
        geom[3][0] = group[4]
        
        if group in [['C','C','F','F','F'],['C','C','H','H','H'],['N','C','H','H','H']]:
            geom[1][1] = str(count+1)
            geom[2][1] = str(count+1) #if also attached to that one
            geom[3][1] = str(count+1)
        elif group in [['C','C','O','O','H']]:
            #print "carbonic acid"
            geom[1][1] = str(count+1) # this is the =O
            geom[2][1] = str(count+1) # this is the -O-

            geom[1][2] = 1.3 # bond length C=O
            geom[1][4] = 120.1 # angle coreC-C=O
            if dihedral:
                geom[1][6]=dihedral
            else:
                geom[1][6] = 0.1 # dihedral with one of the core
            
            geom[2][2] = 1.3 # bond length C-O
            geom[2][4] = 120.1 # angle coreC-C=O
            if dihedral:
                geom[2][6]=dihedral-180.0
            else:
                geom[2][6] = 180.1 # dihedral with one of the core

            geom[3][1] = str(count+3) # attached to -O-
            geom[3][3] = str(count+1) # angled with -C
            geom[3][5] = str(count+2) # dihedraled with =O

            geom[3][2] = 0.9
            geom[3][4] = 109.5
            geom[3][6] = 2.1 # just to make it changebla i don't make it zero

        count += 4
    # when length == 4 then it is an amine or nitro group
    elif len(group) == 4:
        # need to remove one of the hydrogens
        del geom[3]
        # CNHH or CNOO
        geom[0][0] = group[1]
        geom[1][0] = group[2]
        geom[2][0] = group[3]
        
        geom[1][1] = str(count+1)
        geom[2][1] = str(count+1) #if also attached to that one
        if group == ['C','N','O','O']:
            #geom[0][2] = '1.5'
            #geom[1][2] = '1.227'
            #geom[2][2] = '1.227'
            #geom[1][4] = '117.02'
            #geom[2][4] = '117.02'
            #geom[2][6] = '179.5'
            zma = [ ['N',1,'1.51'],
                    ['O',2,'1.227',1,'117.01',0, '30.1'],
                    ['O',2,'1.227',1,'117.01',3,'178.5']]
            if dihedral: zma[1][6]=dihedral
            geom=geomfiller(zma,geom,count)
            #if dihedral: 
            #    geom[1][6]=dihedral
            #    geom[2][6]=dihedral-180.0
        if group == ['C','N','H','H']:
            geom[0][2] = '1.465'
            geom[1][2] = '1.02'
            geom[2][2] = '1.02'
            if dihedral: 
                geom[1][6]=dihedral
                geom[2][6]=dihedral-120.0
        if group in [['C','C','H','O'],['C','C','O','H']]:
            geom[1][1] = str(count+1) # this is the =O
            geom[2][1] = str(count+1) # this is the -H

            geom[1][2] = 1.18 # bond length C=O
            geom[1][4] = 125.1 # angle coreC-C=O
            geom[1][6] = 0.1 # dihedral with one of the core
            
            geom[2][2] = 1.11 # bond length C-O
            geom[2][4] = 115.1 # angle coreC-C=O
            geom[2][6] = 180.1 # dihedral with one of core

        count += 3 # three atoms added to activemat
    elif len(group) == 3:
        # when length is 3 it is a COH, NOH? or CCN group
        # we need to remove two hydrogens
        del geom[3]
        del geom[2]
        if group == ['C','C','N']:
            #print "cyano"
            geom[1] = geom[0][:]
            geom[1][2] = '2.62'
        else:
            # because the N now directly binds to the core
            # geom[1][1] verwijst naar the position of geom[0]
            geom[1][1] = str(count+1) #the bondlength index has to count+1 denk ik
        geom[0][0] = group[1] # change C to O or also C
        geom[1][0] = group[2] # change H to N or also H
        count += 2 # because two atoms are added to the activemat
        # shorten the C-O bond a bit
        if group == ['C','O','H']:
            geom[0][2] = '1.42'
            geom[1][2] = '0.97'
            if dihedral: geom[1][6]=dihedral
    elif len(group) == 2:
        # it is an C-H,C-F,Si-H or C-Cl group. # or it is C Ph
        if group == ['C','Ph']:
            #print "benzene group"
            zma = [ ['C',1,'1.51'],
                    ['C',2,'1.401',1,'120.1',0, '30.1'],
                    ['C',3,'1.395',2,'120.9',1,'178.5'],
                    ['C',4,'1.395',3,'120.1',2,  '0.1'],
                    ['C',5,'1.395',4,'119.5',3,  '0.1'],
                    ['C',6,'1.395',5,'120.1',4,  '0.1'],
                    ['H',3,'1.087',2,'119.3',7,'180.1'],
                    ['H',4,'1.086',3,'119.8',2,'180.1'],
                    ['H',5,'1.086',4,'119.9',3,'180.1'],
                    ['H',6,'1.086',5,'120.1',4,'180.1'],
                    ['H',7,'1.087',6,'119.7',5,'180.1'] ]
            logging.debug('zma benzene:')
            logging.debug(pprint.pformat(zma))
            if dihedral: zma[1][6]=dihedral
            geom=geomfiller(zma,geom,count)
            count += len(zma)
        else:
            # all the three hydrogens need to be removed
            del geom[3]
            del geom[2]
            del geom[1]
            # and thange the C to H or Cl or F 
            geom[0][0] = group[1]
            if group == ['C','O']: # changes angle. analyses the different cases
                #print "keton"
                geom[0][2] = '1.215'
                geom[0][4] = '120.5'
                if -64.0 <= float(geom[0][6]) <= -60.0  or -170.0 <= float(geom[0][6]) <= -164.0:
                    geom[0][6] = '-117.1'
                if 62.0 <= float(geom[0][6]) <= 65.0  or 170.0 <= float(geom[0][6]) <= 178.0:
                    geom[0][6] = '122.0' 
            count += 1
            # the first atom is directly attached to the core so the zmat indices don't change
    elif len(group) == 1:
    # when the length of the group is one, that means it is doped with N,O,S,B,P
        del geom[3]
        del geom[2]
        del geom[1]
        del geom[0]
    # so whole geom is deleted actualy
    # there will not appear any of this ones in the activemat so count is not changed
    return geom,count


if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    multcharge = re.compile('^\-?[0-9]\s[0-9]')
    with open(filename) as fid:
        for line in fid:
            if multcharge.match(line):
                print "match!"
                print multcharge.match(line).group()
                zmat=[]
                line=next(fid)
                while not line == '\n':
                    zmat.append(line.split())
                    line=next(fid)
    ncore = 12
    for pos in positions:
        hornot = maker(zmat,pos,ncore)
        if not hornot == 1:
            maker2(zmat,pos,ncore)
        del zmatnew
        del ind
        del hline
        del item   
    print "DONE"

