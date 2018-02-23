""" a conveniece Molecule Class """
import numpy as np
from converter import Converter
from collections import MutableSequence

debug=0

try:
    pass
    #import pybel
except:
    print "pybel not installed"

def contoind(conf):
    ''' to convert a configuration to an index-string '''
    #return '_'.join([''.join(str(item)) for item in conf])
    return '_'.join([ ''.join(filter(lambda x:str(x).isalpha(), item)) for item in conf ])

def contoindD(conf):
    ''' to convert a configuration to an index-string '''
    #return '_'.join([''.join(str(item)) for item in conf])
    return '_'.join([ ''.join(item) for item in conf ])

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

class Population(MutableSequence):
    def __init__(self, population):
        self.population = population
        self.indices = [ individual.index for individual in self.population ]
        return

    def index(self, key):
        return self.indices.index(key)

    def __iter__(self):
        return iter(self.population)

    def __len__(self):
        return len(self.population)

    def __repr__(self):
        ret = ""
        for individual in self.population:
            ret += '    '
            ret += individual.__repr__()
            ret += '\n'
        return "Population class with:\n" + ret +  "####### END POPULATION #######"

    def __getitem__(self,value):
        '''this gives the class dict and list behavior calling it like:
            Population['index'] or Population[i] both will work '''
        if type(value)==str:
            ii = self.get_indices().index(value)
            return self.population[ii]
        else:
            return self.population[value]

    def get_indices(self):
        return [ individual.index for individual in self.population ]

    def log(self):
        return dict( [ [ I.index, I.Pvalue ] for I in self.population ] )

    def __delitem__(self,i):
        del self.population[i]

    def insert(self,ii, val):
        self.population.insert(ii, val)

    def __add__(self, Pop):
        for ind in Pop:
            self.population.append(ind)
        return self

class BaseMolecule(object):
    def __init__(self):
        self.converter = Converter()
        self.Pvalue = None # for storing the principal properties
        self.boundaries = [] # for storing the boundary condition properties
        self.infoline = [] # for storing additional properties
        self.predictions = {}
        self.predicted = None
        self.opt = False
        # for jsonification:
        self.props = {}
        return

    def __str__(self):
        return "BaseMolecule"

    def __repr__(self):
        return self.__str__()

    def addjob(self, jobname):
        if hasattr(self, 'jobs'): self.jobs.append(jobname)
        else: self.jobs=[jobname]

    def deletejobs(self):
        try: delattr(self, 'jobs')
        except AttributeError: pass

    def submit_jobs(self):
        pass

class SmiMolecule(BaseMolecule):
    def __init__(self, smiles):
        super(SmiMolecule, self).__init__()
        self.smiles=smiles
        self.set_index()
        return

    def __str__(self):
        return "SmiMolecule: " + self.smiles

    def set_index(self):
        replacements = {
                '=':'a',
                '(':'d', ')':'e',
                '[':'g', ']':'i',
                '\\':'j','/':'k',
                '@':'m','-':'q', '+':'r',
                '.':'t'
                }
        #ireplacements = {v: k for k, v in replacements.iteritems()}
        index = "".join([replacements.get(c, c) for c in self.smiles])
        self.index = index
        return

class Molecule(BaseMolecule):
    def __init__(self, conf, dihedral=False):
        super(Molecule, self).__init__()

        # Molecule instances can be initiated with both their index or conf but conf is preferred
        if type(conf)==list:
            self.conf = conf
            if dihedral:
                self.index=contoindD(self.conf)
            else:
                self.index= contoind(self.conf)
        elif type(conf)==str:
            self.conf = indtocon(conf)
            self.index= conf
        self.mat = None
        return

    def __getitem__(self,key):
        return self.conf[key]

    def __contains__(self,value):
        return value in self.conf

    def __len__(self):
        return len(self.conf)

    def __str__(self):
        return "Molecule: " + self.index

    def __eq__(self, other):
        return self.index==other.index

    def copy(self):
        new_mol = Molecule(self.conf)
        new_mol.Pvalue = self.Pvalue
        new_mol.props = self.props
        new_mol.predicted = self.predicted
        return new_mol

    def log(self):
        ret = [ self.index ]
        ret.append( int( not self.predicted ) )
        ret.append( self.Pvalue     )
        ret.extend( self.boundaries )
        #ret.extend( self.infoline   )
        return ret

    def set_path(self, path, extension='.com'):
        pass

    def make_zmat(self, core, active, passive):
        pass

    def set_framework(self,core,active,passive):
        self.core   =core
        self.active =active
        self.passive=passive
        zmat= []
        zmat.extend( core )
        rest = sorted( active+passive , key=lambda x:x[1][1] )
        zmat.extend( [ item for sublist in rest for item in sublist ] )
        self.zmat = zmat
        print "SMILES of framework:", self.get_format()
        return

    def set_xyz(self,xyz):
        # format example: with numpy arrays
        # [['C', array([ 0.729536,  1.425402, -0.080846]), 12.011], 
        # ['H', array([ 0.723016,  1.759466, -1.118944]), 1.00794], ['H', array([ 1.721619, 1.577878,  0.343836]), 1.00794], 
        # ['H', array([ 0.009375,  2.007087,  0.494851]), 1.00794] ] 
        self.xyz=xyz
        return

    def set_zmat(self,zmat):
        self.zmat=zmat
        return

    def xyztozmat(self):
        xyz = self.converter.read_xyzlist(self.xyz)
        self.zmat = self.converter.cartesian_to_zmatrix()
        return

    def zmatoxyz(self):
        zmat = self.converter.read_zmalist(self.zmat)
        self.xyz =  self.converter.zmatrix_to_cartesian()
        return

    def set_OBMol(self):
        self.natoms= len(self.xyz)
        try:
            new_format_xyz = [ [item[0]]+ list(map(str,item[1]))  for item in self.xyz ]
        except IndexError:
            new_format_xyz = self.xyz
        OBxyz = []
        OBxyz.append([str(self.natoms)])
        OBxyz.append([self.index])
        for item in new_format_xyz:
            OBxyz.append(item)
        if debug:
            print "OBxyz:",
            for item in OBxyz: print item
        OBxyz_str = '\n'.join([' '.join(item) for item in OBxyz])
        self.OBMol = pybel.readstring('xyz',OBxyz_str)
        return

    def get_format(self,form='smiles'):
        if not hasattr(self,'OBMol'):
            if not hasattr(self,'xyz'):
                if not hasattr(self,'zmat'):
                    raise SystemExit('molecule object has no data')
                self.zmatoxyz()
            if debug: print "self.xyz",self.xyz
            self.set_OBMol()
        if debug: print "self.OBxyz:", self.OBMol
        return self.OBMol.write(form)

    def optimize(self,set=False):
        self.OBMol.localopt(forcefield='mmff94',steps=100)
        new = self.OBMol.write(format='xyz',filename=None, overwrite=False)
        #print "new:", new
        #print "type new string?:", type(new)
        xyz = self.changeformat(new)
        #print xyz
        if set: self.xyz = xyz
        return xyz

    def changeformat(self,new):
        '''change from pure xyz to converter format with numpy arrays'''
        #print repr(new)
        newl = new.split('\n')[2:-1]
        #print "newl", newl
        xyzs = []
        for line in newl: #for each atom
            assert not line=='\n'
            splitted = line.split()
            xyztje = np.zeros([3])
            #print "xyztje", xyztje
            for j in xrange(3): #for x,y,z
                xyztje[j] = splitted[j+1]
            xyzs.append([splitted[0],xyztje,self.converter.masses[splitted[0]]])
        return xyzs

try:
    from qml import compound
except ImportError:
    print "QML not imported!"
else:
  class my_Compound(compound.Compound):
    '''an inherited class of Compound which is different only in the fact that it reads from list input instead of filename input'''

    def read_xyz(self, lines):
        '''this function is overwritten from the inherited class and sees instead of a filename xyzcoordinates'''

        # - f = open(filename, "r")
        # - lines = f.readlines()
        # - f.close()

        # - self.natoms = int(lines[0])
        # + :
        self.natoms = len(lines)
        self.atomtypes = []
        self.nuclear_charges = []
        self.coordinates = []

        self.name = filename

        for line in lines[2:]:
            tokens = line.split()

            if len(tokens) != 4:
                break

            self.atomtypes.append(tokens[0])
            # - self.nuclear_charges.append(NUCLEAR_CHARGE[tokens[0]])
            self.nuclear_charges.append(compound.NUCLEAR_CHARGE[tokens[0]])
            x = float(tokens[1])
            y = float(tokens[2])
            z = float(tokens[3])

            self.coordinates.append(np.array([x, y, z]))

        self.coordinates = np.array(self.coordinates)

'''
Format of converter cartesian
['C', array([ -1.10135287e+00,   1.52693429e+00,   3.31106160e-04]), 12.011]
['C', array([  5.24031346e-01,   1.52693429e+00,   3.31106160e-04]), 12.011]
['C', array([ -1.61480506e+00,  -1.50096181e-02,   3.31106160e-04]), 12.011]
['C', array([-1.08043513, -0.77828059, -1.33143038]), 12.011]
['C', array([ 1.07969834,  0.77922403, -1.33119709]), 12.011]
['C', array([ 0.54487177, -0.75565322, -1.31836097]), 12.011]
'''



