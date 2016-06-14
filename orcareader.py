# -*- coding: utf-8 -*-
from pprint import pprint
import re
import numpy
import utils
import math
import time

import construction 

corresp = {2:46,6:18,7:42,9:34,11:22,12:30}

class prettyfloat(float):
    def __repr__(self):
        return "%-0.4f" % self

class OrcaMolecule():

    def __init__(self,name):
        self.data = []
        self.optdone = []
        self.fid = open(name) # + '.log')
        self.name = name

    def skip_lines(self, sequence):
        """Read trivial line types and check they are what they are supposed to be.
        This function will read len(sequence) lines and do certain checks on them,
        when the elements of sequence have the appropriate values. Currently the
        following elements trigger checks:
            'blank' or 'b'      - the line should be blank
            'dashes' or 'd'     - the line should contain only dashes (or spaces)
            'equals' or 'e'     - the line should contain only equal signs (or spaces)
            'stars' or 's'      - the line should contain only stars (or spaces)
        """
        expected_characters = {
            '-' : ['dashes', 'd'],
            '=' : ['equals', 'e'],
            '*' : ['stars', 's'],
        }
        lines = []
        for expected in sequence:

            # Read the line we want to skip.
            line = next(self.fid)

            # Blank lines are perhaps the most common thing we want to check for.
            if expected in ["blank", "b"]:
                try:
                    assert line.strip() == ""
                except AssertionError:
                    frame, fname, lno, funcname, funcline, index = inspect.getouterframes(inspect.currentframe())[1]
                    parser = fname.split('/')[-1]
                    msg = "In %s, line %i, line not blank as expected: %s" % (parser, lno, line.strip())
                    self.logger.warning(msg)

            # All cases of heterogeneous lines can be dealt with by the same code.
            for character, keys in expected_characters.items():
                if expected in keys:
                    try:
                        assert all([c == character for c in line.strip() if c != ' '])
                    except AssertionError:
                        frame, fname, lno, funcname, funcline, index = inspect.getouterframes(inspect.currentframe())[1]
                        parser = fname.split('/')[-1]
                        msg = "In %s, line %i, line not all %s as expected: %s" % (parser, lno, keys[0], line.strip())
                        self.logger.warning(msg)
                        continue

            # Save the skipped line, and we will return the whole list.
            lines.append(line)

        return lines

    skip_line = lambda self, expected: self.skip_lines( [expected])

    def set_attribute(self, name, value, check=True):
        """Set an attribute and perform a check when it already exists.

        Note that this can be used for scalars and lists alike, whenever we want
        to set a value for an attribute. By default we want to check that
        the value does not change if the attribute already exists, and this function
        is a good place to add more tests in the future.        
        """
        if check and hasattr(self, name):
            try:
                assert getattr(self, name) == value
            except AssertionError:
                self.logger.warning("Attribute %s changed value (%s -> %s)" % (name, getattr(self, name), value))
        setattr(self, name, value)

    def extract(self,coords=0,orbens=0):
      for line in self.fid:

#    def extract(self, self.fid, line):
        """Extract information from the file object self.fid."""

        if line[0:15] == "Number of atoms":
            if not hasattr(self, "natoms"):
                self.natoms = []
            natom = int(line.split()[-1])
            self.natoms.append(natom) 
            #self.set_attribute('natom', natom)

        if line[1:13] == "Total Charge":
            if not hasattr(self, "charges"):
                self.charges = []
            if not hasattr(self, "mults"):
                self.mults = []
            charge = int(line.split()[-1])
            #self.set_attribute('charge', charge)
            self.charges.append(charge)

            line = next(self.fid)

            mult = int(line.split()[-1])
            self.mults.append(mult)
            #self.set_attribute('mult', mult)

        if 'JOBS TO BE PROCESSED THIS RUN' in line:
            njobs = int(line.split()[3])
            self.set_attribute('njobs', njobs)
        #         $ THERE ARE  2 JOBS TO BE PROCESSED THIS RUN $

        # SCF convergence output begins with:
        #
        # --------------
        # SCF ITERATIONS
        # --------------
        # 
        # However, there are two common formats which need to be handled, implemented as separate functions.
        if "SCF ITERATIONS" in line:

            self.skip_line('dashes')

            line = next(self.fid)
            colums = line.split()
            if colums[1] == "Energy":
                self.parse_scf_condensed_format(colums)
            elif colums[1] == "Starting":
                self.parse_scf_expanded_format(colums)

        # Information about the final iteration, which also includes the convergence
        # targets and the convergence values, is printed separately, in a section like this:
        #
        #       *****************************************************
        #       *                     SUCCESS                       *
        #       *           SCF CONVERGED AFTER   9 CYCLES          *
        #       *****************************************************
        #
        # ...
        #
        # Total Energy       :         -382.04963064 Eh          -10396.09898 eV
        #
        # ...
        #
        # -------------------------   ----------------
        # FINAL SINGLE POINT ENERGY     -382.049630637
        # -------------------------   ----------------
        #
        # We cannot use this last message as a stop condition in general, because
        # often there is vibrational output before it. So we use the 'Total Energy'
        # line. However, what comes after that is different for single point calculations
        # and in the inner steps of geometry optimizations.
        if "SCF CONVERGED AFTER" in line:

            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            #if not hasattr(self, "scfvalues"):
            #    self.scfvalues = []
            #if not hasattr(self, "scftargets"):
            #    self.scftargets = []

            while not "Total Energy       :" in line:
                line = next(self.fid)
            energy = float(line.split()[5])
            self.scfenergies.append(energy)

            #self._append_scfvalues_scftargets( line)

        # Sometimes the SCF does not converge, but does not halt the
        # the run (like in bug 3184890). In this this case, we should
        # remain consistent and use the energy from the last reported
        # SCF cycle. In this case, ORCA print a banner like this:
        #
        #       *****************************************************
        #       *                     ERROR                         *
        #       *           SCF NOT CONVERGED AFTER   8 CYCLES      *
        #       *****************************************************
        if "SCF NOT CONVERGED AFTER" in line:
            print "SCF NOT CONVERGED!"
            #if not hasattr(self, "scfenergies"):
            #    self.scfenergies = []
            #if not hasattr(self, "scfvalues"):
            #    self.scfvalues = []
            #if not hasattr(self, "scftargets"):
            #    self.scftargets = []
            #energy = self.scfvalues[-1][-1][0]
            #self.scfenergies.append(energy)
            #self._append_scfvalues_scftargets(line)  

#         *************************************************************
#         *                GEOMETRY OPTIMIZATION CYCLE   3            *
#         *************************************************************
#---------------------------------
#CARTESIAN COORDINATES (ANGSTROEM)
#---------------------------------
#  C      0.000000    0.000000   -0.003840
#  O      0.000000    0.000000    1.133840
#
#----------------------------
#CARTESIAN COORDINATES (A.U.)
#----------------------------
#  NO LB      ZA    FRAG    MASS        X           Y           Z
#   0 C     6.0000    0    12.011          0.000000000000000          0.000000000000000         -0.007256002283802
#   1 O     8.0000    0    15.999          0.000000000000000          0.000000000000000          2.142646533614816
#
#--------------------------------
#INTERNAL COORDINATES (ANGSTROEM)
#--------------------------------
# C      0   0   0   0.000000     0.000     0.000
# O      1   0   0   1.137679     0.000     0.000
#
#---------------------------
#INTERNAL COORDINATES (A.U.)
#---------------------------
# C      0   0   0   0.000000     0.000     0.000
# O      1   0   0   2.149903     0.000     0.000

        if line[0:32] == "CARTESIAN COORDINATES (ANGSTROEM)":
            if not hasattr(self, "atomcoords"):
                self.atomcoords = []
            self.skip_line('dashes')
            atomnos = []
            atomcoord = []
            line = next(self.fid)
            line = next(self.fid)
            while len(line) > 1:
                broken = line.split()
                atomnos.append(broken[0])
                atomcoord.append(list(map(float, broken[1:4])))
                line = next(self.fid)
            #self.set_attribute('natom', len(atomnos))
            self.set_attribute('atomnos', atomnos)
            self.atomcoords.append(atomcoord)

        if "GEOMETRY OPTIMIZATION CYCLE" in line:
            # Keep track of the current cycle jsut in case, because some things
            # are printed differently inside the first/last and other cycles.
            self.gopt_cycle = int(line.split()[4])
            self.skip_lines(['s', 'd', 'text', 'd'])
            if not hasattr(self,"atomcoords"):
                self.atomcoords = []
            atomnos = []
            atomcoords = []
            for i in range(self.natoms[-1]):
                line = next(self.fid)
                broken = line.split()
                atomnos.append(broken[0])
                atomcoords.append(list(map(float, broken[1:4])))
            self.atomcoords.append(atomcoords)
            self.set_attribute('atomnos', atomnos)

#                 *******************************************************
#                 *** FINAL ENERGY EVALUATION AT THE STATIONARY POINT ***
#                 ***               (AFTER    3 CYCLES)               ***
#                 *******************************************************
#---------------------------------
#CARTESIAN COORDINATES (ANGSTROEM)
#---------------------------------
#  C      0.000000    0.000000   -0.003851
#  O      0.000000    0.000000    1.133851
#
#----------------------------
#CARTESIAN COORDINATES (A.U.)
#----------------------------
#  NO LB      ZA    FRAG    MASS        X           Y           Z
#   0 C     6.0000    0    12.011          0.000000000000000          0.000000000000000         -0.007276828886230
#   1 O     8.0000    0    15.999          0.000000000000000          0.000000000000000          2.142667360217245
#
#--------------------------------
#INTERNAL COORDINATES (ANGSTROEM)
#--------------------------------
# C      0   0   0   0.000000     0.000     0.000
# O      1   0   0   1.137701     0.000     0.000
#
#---------------------------
#INTERNAL COORDINATES (A.U.)
#---------------------------
# C      0   0   0   0.000000     0.000     0.000
# O      1   0   0   2.149944     0.000     0.000
        if line[21:68] == "FINAL ENERGY EVALUATION AT THE STATIONARY POINT":

            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))
            self.skip_lines(['text', 's', 'd', 'text', 'd']) #CARTESIAN ANGSTROEM
            atomcoords = []
            for i in range(self.natoms[-1]):
                line = next(self.fid)
                broken = line.split()
                atomcoords.append(list(map(float, broken[1:4])))

            self.atomcoords.append(atomcoords)

#FINAL SINGLE POINT ENERGY      -703.230614805462
        if 'FINAL SINGLE POINT ENERGY' in line:
            if not hasattr(self,'scfenergies'):
                self.scfenergies = []
            self.scfenergies.append(float(line.split()[4]))

        if "The optimization did not converge" in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []


        atomcharge = {}
        atomspin = {}
        mul=0
        low=0
        has_spins=0
        if line[:23] == "MULLIKEN ATOMIC CHARGES":
            mul = 1
            has_spins = "AND SPIN POPULATIONS" in line
            if not hasattr(self, "atomcharges"):
                self.atomcharges = []
            if has_spins and not hasattr(self, "atomspins"):
                self.atomspins = []
            self.skip_line( 'dashes')
            charges = []
            if has_spins:
                spins = []
            line = next(self.fid)
            while line[:21] != "Sum of atomic charges":
                charges.append(float(line[8:20]))
                if has_spins:
                    spins.append(float(line[20:]))
                line = next(self.fid)
            atomcharge["mulliken"] = charges
            if has_spins:
                atomspin["mulliken"] = spins    
        # Things are the same for Lowdin populations, except that the sums
        #   are not printed (there is a blank line at the end).
        if line[:22] == "LOEWDIN ATOMIC CHARGES":
            low = 1
            has_spins = "AND SPIN POPULATIONS" in line
            if not hasattr(self, "atomcharges"):
                self.atomcharges = []
            if has_spins and not hasattr(self, "atomspins"):
                self.atomspins = []
            self.skip_line('dashes')
            charges = []
            if has_spins:
                spins = []
            line = next(self.fid)
            while line.strip():
                charges.append(float(line[8:20]))
                if has_spins:
                    spins.append(float(line[20:]))
                line = next(self.fid)
            atomcharge["lowdin"] = charges
            if has_spins:
                atomspin["lowdin"] = spins
        if mul==1 or low==1:
            self.atomcharges.append(atomcharge)
            if has_spins:
                self.atomspins.append(atomspin)

        if line[0:16] == "ORBITAL ENERGIES":
            self.skip_lines(['d', 'text', 'text'])
            self.moenergies = [[]]
            self.homos = [[0]]
            line = next(self.fid)
            while len(line) > 20: #restricted calcs are terminated by ------
                info = line.split()
                self.moenergies[0].append(float(info[3]))
                if float(info[1]) > 0.00: #might be 1 or 2, depending on restricted-ness
                    self.homos[0] = int(info[0])
                line = next(self.fid)
            line = next(self.fid)
            #handle beta orbitals
            if line[17:35] == "SPIN DOWN ORBITALS":
                text = next(self.fid)
                self.moenergies.append([])
                self.homos.append(0)
                line = next(self.fid)
                while len(line) > 20: #actually terminated by ------
                    info = line.split()
                    self.moenergies[1].append(float(info[3]))
                    if float(info[1]) == 1.00:
                        self.homos[1] = int(info[0])
                    line = next(self.fid)

        if line.strip() == "DIPOLE MOMENT":
            self.skip_lines(['d', 'XYZ', 'electronic', 'nuclear', 'd'])
            total = next(self.fid)
            assert "Total Dipole Moment" in total
            reference = [0.0, 0.0, 0.0]
            dipole = numpy.array([float(d) for d in total.split()[-3:]])
            dipole = utils.convertor(dipole, "ebohr", "Debye")
            if not hasattr(self, 'moments'):
                self.moments = [reference, dipole]
            else:
                try:
                    assert numpy.all(self.moments[1] == dipole)
                except AssertionError:
                    self.moments = [reference, dipole]


    def parse_scf_condensed_format(self, line):
        """ Parse the SCF convergence information in condensed format """

        # This is what it looks like
        # ITER       Energy         Delta-E        Max-DP      RMS-DP      [F,P]     Damp
        #                ***  Starting incremental Fock matrix formation  ***
        #   0   -384.5203638934   0.000000000000 0.03375012  0.00223249  0.1351565 0.7000
        #   1   -384.5792776162  -0.058913722842 0.02841696  0.00175952  0.0734529 0.7000
        #                                ***Turning on DIIS***
        #   2   -384.6074211837  -0.028143567475 0.04968025  0.00326114  0.0310435 0.0000
        #   3   -384.6479682063  -0.040547022616 0.02097477  0.00121132  0.0361982 0.0000
        #   4   -384.6571124353  -0.009144228947 0.00576471  0.00035160  0.0061205 0.0000
        #   5   -384.6574659959  -0.000353560584 0.00191156  0.00010160  0.0025838 0.0000
        #   6   -384.6574990782  -0.000033082375 0.00052492  0.00003800  0.0002061 0.0000
        #   7   -384.6575005762  -0.000001497987 0.00020257  0.00001146  0.0001652 0.0000
        #   8   -384.6575007321  -0.000000155848 0.00008572  0.00000435  0.0000745 0.0000
        #          **** Energy Check signals convergence ****
        assert line[2] == "Delta-E"
        assert line[3] == "Max-DP"

        if not hasattr(self, "scfvalues"):
            self.scfvalues = []

        #self.scfvalues.append([])

        # Try to keep track of the converger (NR, DIIS, SOSCF, etc.).
        diis_active = True
        while not line == []:

            if 'Newton-Raphson' in line:
                diis_active = False
            elif 'SOSCF' in line:
                diis_active = False
            elif line[0].isdigit() and diis_active:
                energy = float(line[1])
                #deltaE = float(line[2])
                #maxDP = float(line[3])
                #rmsDP = float(line[4])
                #self.scfvalues[-1].append([deltaE, maxDP, rmsDP])
                #self.scfvalues[-1].append([deltaE])
            elif line[0].isdigit() and not diis_active:
                energy = float(line[1])
                #deltaE = float(line[2])
                #maxDP = float(line[5])
                #rmsDP = float(line[6])
                #self.scfvalues[-1].append([deltaE, maxDP, rmsDP])
                #self.scfvalues[-1].append([deltaE])
            line = next(self.fid).split()

    def parse_scf_expanded_format(self, line):
        """ Parse SCF convergence when in expanded format. """


# The following is an example of the format
# -----------------------------------------
#
#               ***  Starting incremental Fock matrix formation  ***
#
#                         ----------------------------
#                         !        ITERATION     0   !
#                         ----------------------------
#   Total Energy        :    -377.960836651297 Eh
#   Energy Change       :    -377.960836651297 Eh
#   MAX-DP              :       0.100175793695
#   RMS-DP              :       0.004437973661
#   Actual Damping      :       0.7000
#   Actual Level Shift  :       0.2500 Eh
#   Int. Num. El.       :    43.99982197 (UP=   21.99991099 DN=   21.99991099)
#   Exchange            :   -34.27550826
#   Correlation         :    -2.02540957
#
#
#                         ----------------------------
#                         !        ITERATION     1   !
#                         ----------------------------
#   Total Energy        :    -378.118458080109 Eh
#   Energy Change       :      -0.157621428812 Eh
#   MAX-DP              :       0.053240648588
#   RMS-DP              :       0.002375092508
#   Actual Damping      :       0.7000
#   Actual Level Shift  :       0.2500 Eh
#   Int. Num. El.       :    43.99994143 (UP=   21.99997071 DN=   21.99997071)
#   Exchange            :   -34.00291075
#   Correlation         :    -2.01607243
#
#                               ***Turning on DIIS***
#
#                         ----------------------------
#                         !        ITERATION     2   !
#                         ----------------------------
# ....
#
        if not hasattr(self, "scfvalues"):
            self.scfvalues = []

        #self.scfvalues.append([])
        line = "Foo" # dummy argument to enter loop
        while line.find("******") < 0:
            line = next(self.fid)
            info = line.split()
            if len(info) > 1 and info[1] == "ITERATION":
                dashes = next(self.fid)
                energy_line = next(self.fid).split()
                energy = float(energy_line[3])
                #deltaE_line = next(self.fid).split()
                #deltaE = float(deltaE_line[3])
                #if energy == deltaE:
                #    deltaE = 0
                #maxDP_line = next(self.fid).split()
                #maxDP = float(maxDP_line[2])
                #rmsDP_line = next(self.fid).split()
                #rmsDP = float(rmsDP_line[2])
                #self.scfvalues[-1].append([deltaE, maxDP, rmsDP])

        return

    # end of parse_scf_expanded_format

    def _append_scfvalues_scftargets(self, line):
        # The SCF convergence targets are always printed after this, but apparently
        # not all of them always -- for example the RMS Density is missing for geometry
        # optimization steps. So, assume the previous value is still valid if it is
        # not found. For additional certainty, assert that the other targets are unchanged.
        while not "Last Energy change" in line:
            line = next(self.fid)
        deltaE_value = float(line.split()[4])
        deltaE_target = float(line.split()[7])
        line = next(self.fid)
        if "Last MAX-Density change" in line:
            maxDP_value = float(line.split()[4])
            maxDP_target = float(line.split()[7])
            line = next(self.fid)
            if "Last RMS-Density change" in line:
                rmsDP_value = float(line.split()[4])
                rmsDP_target = float(line.split()[7])
            else:
                rmsDP_value = self.scfvalues[-1][-1][2]
                rmsDP_target = self.scftargets[-1][2]
                assert deltaE_target == self.scftargets[-1][0]
                assert maxDP_target == self.scftargets[-1][1]
            self.scfvalues[-1].append([deltaE_value, maxDP_value, rmsDP_value])
            self.scftargets.append([deltaE_target, maxDP_target, rmsDP_target])

def orcaread(filename,prop,multiplejobs=2,rdvindex=1):
    mymol = OrcaMolecule(filename)
    mymol.extract()
    if mymol.optdone == []:
        print "program did not do optimization or crashed in earlier state"
    #if prop=='polar':
    #    if hasattr(mymol,'polex'):
    #        print "polar exact densities:"
    #        pprint(mymol.polex)
    #        pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
    #        # next line calculates in one line the Radical delocalisation value
    #        # RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
    #        print "polarisability value is: ", pola
    #        return pola
    #    else:
    #        print "NO POLAR DATA FOUND IN FILE!"
    #        return
    if prop=='stabA':
        pprint(mymol.scfenergies[-5:])
        Eopt= mymol.scfenergies[-4]
        E0= mymol.scfenergies[-3]
        IP= mymol.scfenergies[-2] - E0
        EA= E0 - mymol.scfenergies[-1]
        return (Eopt,E0,IP,EA) #returns in EV?
    elif prop=='energy':
        ESCFs = mymol.scfenergies
        return ESCFs[-1]
    elif prop in ['dipole','mu']:
        dipole = mymol.moments[1]
        total = math.sqrt( sum( [ item**2 for item in dipole ] ) )
        return total
    elif prop=='gap':
        from cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        Ehomo= myfile.moenergies[0][myfile.homos[0]] #NOT IMPLEMENTED NEED CCLIB
        Elumo= myfile.moenergies[0][myfile.homos[0]+1] #NOT IMPLEMENTED NEED CCLIB
        Egap = Elumo - Ehomo
        return Egap
    elif prop in ['IP','ip']:
        ESCFs = mymol.scfenergies
        E0 = ESCFs[-multiplejobs]
        EIP = ESCFs[-multiplejobs+1]
        I = EIP - E0
        print "I:",I
        return I
    elif prop in ['EA','ea']:
        ESCFs = mymol.scfenergies
        E0 = ESCFs[-multiplejobs]
        EA = ESCFs[-1]
        A = E0-EA
        print A, "= A"
        return A
    elif prop in ['natoms','natom']:
        return len(mymol.atomnos)
    elif prop in ['rdv','RDV']:
        print "Mulliken spin densities:"
        type1 = 'mulliken'
        type2 = 'lowdin'
        if type1 in mymol.atomspins[-1]:
            spiden = mymol.atomspins[-1][type1]
        elif type1 in mymol.atomspins[-2]:
            spiden = mymol.atomspins[-2][type1]
        else:
            print "no" + type1 + "found in output switch to" + type2
            if type2 in mymol.atomspins[-1]:
                spiden = mymol.atomspins[-1][type2]
            elif type2 in mymol.atomspins[-2]:
                spiden = mymol.atomspins[-2][type2]
            else:
                print "no RDV found in output!"
                return 0
        pprint(spiden)
        spiden= mymol.spindensities[rdvindex]
        # next line calculates in one line the Radical delocalisation value
        RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
        print "RDV value is: ", RDV
        return RDV
    else:
        print "property not recognised!"
        return

def datareader(indices,jobids,path,data,fileparameters):
    files=[]
    end = '.out'
    for i in range(len(indices)): #make a list of paths from which the data has to be extracted
        files.append(path + '/' + fileparameters['identify'] + indices[i] + end)
        if fileparameters['stab']==1:
            for pos in fileparameters['positions']: #extract al AH energies and take the lowest
                files.append(path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + end)
    normaltermination(files) #test normal termination of all the files
    for i in range(len(indices)):
        file1 = path + '/' + fileparameters['identify'] + indices[i] + end
        #------
        if fileparameters['stab']:#if so we are sure we have to do the following
            #---- some parameters needed
            bde_a = -12.68 #kJ/mol/eV^2
            bde_b = -218.1 #kJ/mol
            stab_h = 235.8 #kJ/mol
            Dw_h = 0.063 #eV
            chi_h = 2.20 
            chi_c = 2.60
            chi_n = 3.05
            H_h = -0.516817233 #a.u.
            kJmol = 2625.5
            eV = 27.2113838
            avtc = -28.1290706 #kJ/mol #average thermal correction for 5 random structures kJ/mol
            chi_term = bde_b*(chi_h-3)*(chi_n-3) #term is independent of the molecule itself. ongeveer 8.4 kJ/mol?
            #---- 
            (Eopt,E0,I,A)=orcaread(file1,'stabA')
            EAHs=[] #all EAHs from all different positions in here
            Npos=[] #positions with a nitrogen in here
            for pos in fileparameters['positions']: #extract al AH energies and take the lowest
                file2= path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + end
                EAH= orcaread(file2,'energy')
                #---- HERE THE electronegativity part of the stab a bit tricky
                confje = construction.indtocon(indices[i]) # change index to conf list using the construction module
                corresp = {2:46,6:18,7:42,9:34,11:22,12:30} # map the alpha positions to methyl indices 
                siteindex = corresp[pos] #find for each position the methyl index
                if siteindex in fileparameters['line1']:# look if that index is used as a site
                    # fileparameters['line1'].index(siteindex) is the place. this is the same as in confje
                    if confje[fileparameters['line1'].index(siteindex)]==['N']:
                        print "electronegativity correction for nitrogen: ", chi_term
                        Npos.append(pos)
                        #EAH -= chi_term
                #-----
                EAHs.append([EAH,pos])
            print "EAHs:"
            pprint(EAHs)
            print "min EAHs:", min(EAHs)
            print "Npos:",Npos
            #gasconstant = 8.3144621
            E_ah=min(EAHs)
            omega= ( ((I+A)**2 )/(8*(I-A)) )*eV #in eV
            Domega = omega - 2
            BDE_ah = (Eopt + H_h - E_ah[0])*kJmol + avtc #avtc is AVerage Thermal Correction. 
            #----
            RDV = orcaread(file1,'RDV',2)
            #----
            if E_ah[1] in Npos:
                stabA= BDE_ah - stab_h - bde_a * Domega * Dw_h - chi_term 
            else:
                stabA= BDE_ah - stab_h - bde_a * Domega * Dw_h 
            #----------
            if 'bcprop' in fileparameters:#decide how to put the data in the datalist
                if fileparameters['property']=='stab': #optimize stab and use another prop as bc
                     if fileparameters['bcprop'] in ['ip','IP','I']: 
                         propy = I
                     elif fileparameters['bcprop'] in ['ea','EA','A']: 
                         propy = A
                     else:
                         propy = orcaread(file1,fileparameters['bcprop'])
                     data.append([indices[i],stabA,propy,omega,RDV])
                else: #so bcprop is stab so propx is the other property to optimize
                     if fileparameters['property'] in ['ip','IP','I']: 
                         propx = I
                     elif fileparameters['property'] in ['ea','EA','A']: 
                         propx = A
                     else:
                         propx = orcaread(file1,fileparameters['property'])
                     data.append([indices[i],propx,stabA,omega,RDV])
            else: #just simple single stab property optimization
                data.append([indices[i],stabA,BDE_ah,I,A,RDV,E_ah[1]])
        #------
        else:
            propx = orcaread(file1,fileparameters['property']) # later this has to change to EHOMO and ELUMO etc
            if 'bcprop' in fileparameters:
                propy = orcaread(file1,fileparameters['bcprop']) # later this has to change to EHOMO and ELUMO etc
                data.append([indices[i],propx,propy])
            else:
                data.append([indices[i],propx])
    return data

def normaltermination(filepaths):
    #-----
    def termination(filepath):
       with open(filepath,'r') as fid:
           text = fid.readlines()[-3:]
#                             ****ORCA TERMINATED NORMALLY****
           if re.search('TERMINATED NORMALLY',''.join(text)):
               fid.close()
               return 1
    #-----
    for path in filepaths: #test all for information which jobs crashed
        if not termination(path)==1:
            print "Error termination:",path
            errortermination(path)
    extratime = 0
    once = 0
    for path in filepaths: #test one by one waiting for normal termination
        while True:
           if termination(path)==1:
               break
           else:
               print "no normal termination for: ",path
           time.sleep(300) # wait 5 minudtes
           extratime += 300
           print "extra waittime/h:", extratime/3600, "||",
    return 

def errortermination(path):
    mymol=OrcaMolecule(path)
    mymol.extract(coords=1)
    import utils
    t=utils.PeriodicTable()
    if hasattr(mymol,'atomcoords'): 
        for sym,xyz in zip(mymol.atomnos,mymol.atomcoords[-1]):
            xyz.insert(0,sym) 
        print "atomcoords and added elements:"
        for item in mymol.atomcoords[-1]:
            print ' '.join(map(str,item)) 
        if False:
            #test if file contains message about not convergence
            fid = open(path[:-3]+'com','r')
            multcharge = re.compile('^\-?[01]\s[12]')
            newfile=[]
            once=0
            for line in fid:
                if multcharge.match(line) and once==0:
                    once+=1
                    while True:
                        line= next(fid)
                        if line=='\n':
                            newfile.extend([' '.join(map(str,item))+'\n' for item in mymol.atomcoords[-1]])
                            newfile.extend(['\n'])
                            break
                else:
                    newfile.append(line)
            print "="*20
  
            #for line in newfile: print line,
            print newfile
            open('zzz','w').writelines(newfile)
    #import construction see top of file
    #construction.filewriter(mat,index,paras)
    return

if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    mymol = OrcaMolecule(filename)
    print mymol , "mymol"
    print mymol.name , "mymol.name"
    print mymol.fid ,'mymol.fid'
    mymol.extract()
    print "extract done" 
    print "attributes:"
    pprint(vars(mymol))

    if hasattr(mymol,'optdone'):
        print "opt done?:" , mymol.optdone
    if hasattr(mymol,'scfenergies'):
        print "last two of scfenergies"
        pprint(mymol.scfenergies[-3:])
    
    if hasattr(mymol,'atomspins'):
        print "atomspins:"
        pprint(mymol.atomspins)
 
    if hasattr(mymol,"Hcorr"): print "Thermalcorrectionenthalpy:", mymol.Hcorr
    if hasattr(mymol,'spindensities'):
        print "Mulliken spin densities:"
        pprint(mymol.spindensities)
        spiden= mymol.spindensities[0]
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
     
  
