""" a conveniece Molecule Class """
import numpy as np
from converter import Converter
debug=0

if True:
    import pybel

class Molecule(object):
    def __init__(self,name='name'):
        self.converter = Converter()
        self.name = name
        return

    def __str__(self):
        return "I am a molecule"

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
        OBxyz.append([self.name])
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
                if debug: 
                    print "self.zmat",
                    for item in self.zmat: print item
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



'''
Format of converter cartesian
['C', array([ -1.10135287e+00,   1.52693429e+00,   3.31106160e-04]), 12.011]
['C', array([  5.24031346e-01,   1.52693429e+00,   3.31106160e-04]), 12.011]
['C', array([ -1.61480506e+00,  -1.50096181e-02,   3.31106160e-04]), 12.011]
['C', array([-1.08043513, -0.77828059, -1.33143038]), 12.011]
['C', array([ 1.07969834,  0.77922403, -1.33119709]), 12.011]
['C', array([ 0.54487177, -0.75565322, -1.31836097]), 12.011]
'''



