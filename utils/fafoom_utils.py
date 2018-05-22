
import fafoomga
from rdkit import Chem
import os, sys

def GetLowestXYZ(mol):
    print "in GetLowestXYZ"
    for filename in ['kill.dat', 'optimized_structures.sdf']:
        if os.path.exists(filename): os.remove(filename)
    try:
        fafoomga.main(mol.smiles)
    except SystemExit as e:
        print e
    suppl = Chem.SDMolSupplier('optimized_structures.sdf', removeHs=False)
    mols = [ mol for mol in suppl ]
    lowmol = sorted(mols, key=lambda mol:mol.GetDoubleProp('Energy'))[0]
    xyz = XYZFromMol(lowmol)
    print "xyz:", xyz
    return xyz

def GetConformers(mol):
    print "in GetConformers"
    for filename in ['kill.dat', 'optimized_structures.sdf']:
        if os.path.exists(filename): os.remove(filename)
    try:
        fafoomga.main(mol.smiles)
    except SystemExit as e:
        print e
    suppl = Chem.SDMolSupplier('optimized_structures.sdf', removeHs=False)
    mols = [ mol for mol in suppl ]
    mols = sorted(mols, key=lambda mol:mol.GetDoubleProp('Energy'))
    for mol in mols:
        mol.SetProp('xyz', XYZFromMol(mol))
    return mols

def XYZFromMol(mol):
    molAddH = Chem.AddHs(mol)
    molStr = Chem.MolToMolBlock(molAddH)

    #parse mol string file into list
    molCoords = []

    molStr2List = molStr.split('\n')

    #get number of atoms
    num = int(molStr2List[3].split()[0])

    for i in xrange(4,4+num):
        coords = molStr2List[i].split()[0:4]
        molCoords.append(coords)

    reorder   = lambda xyz:"{3} {0} {1} {2}".format(*xyz)
    cartesian = map(reorder, molCoords)
    xyz       = "\n".join(cartesian)
    xyz      += '\n\n'

    return xyz

