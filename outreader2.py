from pprint import pprint
import numpy
class prettyfloat(float):
    def __repr__(self):
        return "%-0.4f" % self

class Molecule():

    def __init__(self,name):
        self.data = []
        self.fid = open(name) # + '.log')
        self.name = name
    def extract(self):
      for line in self.fid:
        if line[1:23] == "Optimization completed":

            if not hasattr(self, 'optdone'):
                self.optdone = True

        # Catch message about stopped optimization (not converged).
        if line[1:21] == "Optimization stopped":
            if not hasattr(self, "optdone"):
                self.optdone = False

        # THERMAL CORRECTION TO ENTHALPY
        if line[1:32] == "Thermal correction to Enthalpy=":
            #if not hasattr(self, "Hcorr"):
                #self.Hcorr = []
            #self.Hcorr.append(float(line.split()[4]))
            self.Hcorr = float(line.split()[4])
   
        # SPIN DENSITIES
        if line[1:32] ==  "Mulliken atomic spin densities:":
            print "Remco"
            if not hasattr(self, "spindensities"):
                self.spindensities = []
            line = next(self.fid)
            line = next(self.fid)
            while line[1:4] != "Sum":
                broken = line.split()
                self.spindensities.append((int(broken[0]),broken[1],float(broken[2])))
                line = next(self.fid)

        # Note: this needs to follow the section where 'SCF Done' is used
        #   to terminate a loop when extracting SCF convergence information.
        if line[1:9] == 'SCF Done':

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []

            self.scfenergies.append(line.split()[4])
        # gmagoon 5/27/09: added scfenergies reading for PM3 case
        # Example line: " Energy=   -0.077520562724 NIter=  14."
        # See regression Gaussian03/QVGXLLKOCUKJST-UHFFFAOYAJmult3Fixed.out
        if line[1:8] == 'Energy=':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(float(line.split()[1]))
        if line.strip('* \n') == 'Alpha spin orbitals':
            line = next(self.fid)
            while not line.strip() == 'Summary of Natural Population Analysis:':
                line = next(self.fid)
            for _ in xrange(6):
                line= next(self.fid)
            if not hasattr(self, "npa"):
                self.npa = []
            while not '=' in line:
                self.npa.append(float(line.split()[2]))
                line = next(self.fid)
        if line.strip('* \n') == 'Beta  spin orbitals': #note 2 spaces!
            line = next(self.fid)
            while not line.strip() == 'Summary of Natural Population Analysis:':
                line = next(self.fid)
            for _ in xrange(6):
                line= next(self.fid)
            if not hasattr(self, "npab"):
                self.npab = []
            while not '=' in line:
                self.npab.append(float(line.split()[2]))
                line = next(self.fid)
             
# Summary of Natural Population Analysis:                  
#                                                          
#                                       Natural Population 
#                Natural  -----------------------------------------------
#    Atom  No    Charge         Core      Valence    Rydberg      Total
# -----------------------------------------------------------------------
#      C    1    0.18075      0.99939     1.80808    0.01177     2.81925
#
        if line[2:22] == 'Exact polarizability':
            if not hasattr(self,'polex'):
                self.polex = [] 
            regel = line.split()
            for item in line.split()[2:8]:
                self.polex.append(float(item))
        if line[1:22] == 'Approx polarizability':
            if not hasattr(self,'polprox'):
                self.polprox = [] 
            regel = line.split()
            for item in line.split()[2:8]:
                self.polprox.append(float(item))
#  Exact polarizability:  11.170   0.000  11.170   0.000   0.000  11.170
# Approx polarizability:   8.585   0.000   8.585   0.000   0.000   8.585
        if line.strip('* :\n') == 'Diagonal vibrational polarizability':
            line = next(self.fid)
            self.polvibr = [ float(item) for item in line.split() ]
        if line.strip('* :\n') == 'Diagonal vibrational hyperpolarizability':
            line = next(self.fid)
            self.hypolvibr = [ float(item) for item in line.split() if not float(item)==0.0 ]
# Diagonal vibrational polarizability:
#        0.3451187       0.3451187       0.3451187
# Diagonal vibrational hyperpolarizability:
#        0.0000000       0.0000000       0.0000000

def gausread(filename):
    mymol = Molecule(filename)
    mymol.extract()
    if not hasattr(mymol,'optdone'):
        print "program did not do optimization or crashed in earlier state"
    if hasattr(mymol,'polex'):
        print "polar exact densities:"
        pprint(mymol.polex)
        pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
        # next line calculates in one line the Radical delocalisation value
        # RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
        print "polarisability value is: ", pola
    return pola

if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    mymol = Molecule(filename)
    print mymol , "mymol"
    print mymol.name , "mymol.name"
    print mymol.fid ,'mymol.fid'
    mymol.extract()
    print "extract done" 
    if hasattr(mymol,'optdone'):
        print "opt done?:" , mymol.optdone
    if hasattr(mymol,'scfenergies'):
        print "last two of scfenergies"
        pprint(mymol.scfenergies[-3:-1])
    if hasattr(mymol,"Hcorr"): print "Thermalcorrectionenthalpy:", mymol.Hcorr
    if hasattr(mymol,'spindensities'):
        print "Mulliken spin densities:"
        pprint(mymol.spindensities)
        spiden= mymol.spindensities
        # next line calculates in one line the Radical delocalisation value
        RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
        print "RDV value is: ", RDV
    if hasattr(mymol,'npa'):
        print "npa:"
        pprint(mymol.npa)
    if hasattr(mymol,'npab'):
        print "npab:"
        pprint(mymol.npa)
        snpa = [ a-b for a,b in zip(mymol.npa,mymol.npab) ]
        print "----- npa spin densities -----"
        for item in snpa: print '{:>8.5f}'.format(float(item))
    if hasattr(mymol,'polex'):
        print " Exact polarisability:", mymol.polex
    if hasattr(mymol,'polprox'):
        print "Approx polarisability:", mymol.polprox
    if hasattr(mymol,'polvibr'):
        print "Diagonal vibrational polarisability:", mymol.polvibr
    if hasattr(mymol,'hypolvibr'):
        if not mymol.hypolvibr==[]: print "Diagonal vibrational hyperpolarisability:", mymol.hypolvibr
     
  
