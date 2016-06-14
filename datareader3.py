from pprint import pprint
import re
import numpy
import time

import construction 

corresp = {2:46,6:18,7:42,9:34,11:22,12:30}

class prettyfloat(float):
    def __repr__(self):
        return "%-0.4f" % self

class Molecule():

    def __init__(self,name):
        self.data = []
        self.fid = open(name) # + '.log')
        self.name = name
    def extract(self,coords=0):
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
            self.Hcorr = float(line.split()[4])
   
        # SPIN DENSITIES
        if line[1:32] ==  "Mulliken atomic spin densities:":
            if not hasattr(self, "spindensities"):
                self.spindensities = []
            spinstance = []
            line = next(self.fid)
            line = next(self.fid)
            while line[1:4] != "Sum":
                broken = line.split()
                spinstance.append((int(broken[0]),broken[1],float(broken[2])))
                line = next(self.fid)
            self.spindensities.append(spinstance)
        if line[1:37] == "Mulliken charges and spin densities:":
            if not hasattr(self, "spindensities"):
                self.spindensities = []
            spinstance = []
            line = next(self.fid)
            line = next(self.fid)
            while line[1:4] != "Sum":
                broken = line.split()
                spinstance.append((int(broken[0]),broken[1],float(broken[3])))
                line = next(self.fid)
            self.spindensities.append(spinstance)

        # Note: this needs to follow the section where 'SCF Done' is used
        #   to terminate a loop when extracting SCF convergence information.
        if line[1:9] == 'SCF Done':
            if not hasattr(self, "scfenergies"):
                self.scfenergies = []
            self.scfenergies.append(float(line.split()[4]))
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

        if coords==1:#only extract the molecular coordinates when asked for
            if line.strip() == "Standard orientation:":
                if not hasattr(self, "atomcoords"):
                    self.atomcoords = []
                for _ in range(4): line=next(self.fid)
                atomnos = []
                atomcoords = []
                line = next(self.fid)
                while list(set(line.strip())) != ["-"]:
                    broken = line.split()
                    atomnos.append(int(broken[1]))
                    atomcoords.append(list(map(float, broken[-3:])))
                    line = next(self.fid)
                self.atomcoords.append(atomcoords)
                self.natom=len(atomnos)
                self.atomnos=atomnos
             
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
        if line[1:13] == 'Molar volume':
            self.volume = float(line.split()[3])
# Molar volume =  230.840 bohr**3/mol ( 20.600 cm**3/mol)
# Recommended a0 for SCRF calculation =  2.72 angstrom (  5.13 bohr)
# Dipole moment
#    X=             -1.4604    Y=             -0.1878    Z=              1.4210  Tot=              2.0462
        if line[1:14] == 'Dipole moment':
            line = next(self.fid)
            self.dipole = float(line.split()[7])
      self.fid.close()

def datareader(indices,jobids,path,data,fileparameters):
    files=[]
    for i in range(len(indices)): #make a list of paths from which the data has to be extracted
        files.append(path + '/' + fileparameters['identify'] + indices[i] + '.log')
        if fileparameters['stab']==1:
            for pos in fileparameters['positions']: #extract al AH energies and take the lowest
                files.append(path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.log')
    normaltermination(files,fileparameters['debug']) #test normal termination of all the files
    for i in range(len(indices)):
        file1 = path + '/' + fileparameters['identify'] + indices[i] + '.log'
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
            (Eopt,E0,I,A)=gausread(file1,'stabA')[0]
            EAHs=[] #all EAHs from all different positions in here
            Npos=[] #positions with a nitrogen in here
            for pos in fileparameters['positions']: #extract al AH energies and take the lowest
                file2= path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.log'
                EAH= gausread(file2,'energy')[0]
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
            RDV = gausread(file1,'RDV',2)[0]
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
                         propy = gausread(file1,fileparameters['bcprop'])[0]
                     data.append([indices[i],stabA,propy,omega,RDV])
                else: #so bcprop is stab so propx is the other property to optimize
                     if fileparameters['property'] in ['ip','IP','I']: 
                         propx = I
                     elif fileparameters['property'] in ['ea','EA','A']: 
                         propx = A
                     else:
                         propx = gausread(file1,fileparameters['property'])[0]
                     data.append([indices[i],propx,stabA,omega,RDV])
            else: #just simple single stab property optimization
                data.append([indices[i],stabA,omega,RDV])
        #------
        else:
            #propx,extra = gausread(file1,fileparameters['property']) # later this has to change to EHOMO and ELUMO etc
            datax = gausread(file1,fileparameters['property']) # later this has to change to EHOMO and ELUMO etc
            (propx,extradata) = (datax[0],datax[1:]) #if no extradata = []
            if 'bcprop' in fileparameters:
                datay = gausread(file1,fileparameters['bcprop']) # later this has to change to EHOMO and ELUMO etc
                (propy,extradata)=(datay[0],datay[1:])
                data.append([indices[i],propx,propy]+extradata) #extradata may be an empty list
            else:
                data.append([indices[i],propx]+extradata)
    return data

def gausread(filename,prop,multiplejobs=0,rdvindex=1):
    mymol = Molecule(filename)
    extra = []
    if prop in ['natom','natoms']:
        mymol.extract(coords=1)
    else:
        mymol.extract()
    if not hasattr(mymol,'optdone'):
        print "program did not do optimization or crashed"
    if prop=='polar':
        if hasattr(mymol,'polex'):
            print "polar exact densities:"
            pprint(mymol.polex)
            pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
            # next line calculates in one line the Radical delocalisation value
            # RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
            print "polarisability value is: ", pola
            data = pola,[]
        else:
            print "NO POLAR DATA FOUND IN FILE!"
    elif prop=='stabA':
        pprint(mymol.scfenergies[-5:])
        Eopt= mymol.scfenergies[-4]
        E0= mymol.scfenergies[-3]
        IP= mymol.scfenergies[-2] - E0
        EA= E0 - mymol.scfenergies[-1]
        data = (Eopt,E0,IP,EA) #returns in EV?
    elif prop=='energy':
        ESCFs = mymol.scfenergies
        data = ESCFs[-1]
    elif prop in ['homo','HOMO']:
        from cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        Ehomo= myfile.moenergies[0][myfile.homos[0]] #NOT IMPLEMENTED NEED CCLIB
        Elumo= myfile.moenergies[0][myfile.homos[0]+1] #NOT IMPLEMENTED NEED CCLIB
        data = Ehomo
        extra.extend([Elumo])
    elif prop in ['lumo','LUMO']:
        from cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        Ehomo= myfile.moenergies[0][myfile.homos[0]] #NOT IMPLEMENTED NEED CCLIB
        Elumo= myfile.moenergies[0][myfile.homos[0]+1] #NOT IMPLEMENTED NEED CCLIB
        data = Elumo
        extra.extend([Ehomo])
    elif prop=='gap':
        from cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        Ehomo= myfile.moenergies[0][myfile.homos[0]] #NOT IMPLEMENTED NEED CCLIB
        Elumo= myfile.moenergies[0][myfile.homos[0]+1] #NOT IMPLEMENTED NEED CCLIB
        Egap = Elumo - Ehomo
        #return Egap, (Ehomo,Elumo)
        data = Egap
        extra.extend([Ehomo,Elumo])
    elif prop in ['IP','ip']:
        ESCFs = mymol.scfenergies
        E0 = ESCFs[-(1+multiplejobs)]
        EIP = ESCFs[-multiplejobs]
        I = EIP - E0
        print "I:",I
        data = I
        extra.extend([E0,EIP])
    elif prop in ['EA','ea']:
        ESCFs = mymol.scfenergies
        E0 = ESCFs[-(1+multiplejobs)]
        EA = ESCFs[-1]
        A = E0-EA
        print A, "= A"
        data = A
        extra.extend([E0,EA])
    elif prop in ['natoms','natom']:
        data = len(mymol.atomnos)
    elif prop=='volume':
        data = mymol.volume
    elif prop=='poldens':
        pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
        print "polarisability value is: ", pola
        data = pola/mymol.volume
    elif prop in ['rdv','RDV']:
        print "Mulliken spin densities:"
        pprint(mymol.spindensities)
        spiden= mymol.spindensities[rdvindex]
        # next line calculates in one line the Radical delocalisation value
        RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
        print "RDV value is: ", RDV
        data =  RDV
    elif prop in ['dipole','Dipole']:
        data= mymol.dipole
    else:
        print "property not recognised!"
    return [data]+extra

def normaltermination(filepaths,debug=False):
    #-----
    def termination(filepath):
       with open(filepath,'r') as fid:
           text = fid.readlines()[-3:]
           if re.search('Normal termination',''.join(text)):
               fid.close()
               return 1
    #-----
    copyfilepaths = filepaths[:] #copy to be able to append to it while looping over it
    for path in copyfilepaths: #test all for information which jobs crashed
        if not termination(path)==1:
            print "Error termination:",path
            bnewfile = errortermination(path,debug)
            #if bnewfile and debug:
            #    filepaths.append(path[:-4]+'zzz.com')
    extratime = 0
    once = 0
    for path in filepaths: #test one by one waiting for normal termination
        while True:
           if termination(path)==1:
               break
           else:
               print "no normal termination for: ",path
           #time.sleep(300) # wait 5 minudtes
           time.sleep(300) # wait 5 minudtes
           extratime += 300
           print "extra waittime/h:", extratime/3600, "||",
    return 

def errortermination(path,debug=False):
    mymol=Molecule(path) #read outputfile
    mymol.extract(coords=1) #extract file with also the coordinates
    import utils
    t=utils.PeriodicTable()
    if hasattr(mymol,'atomcoords'): 
        for sym,xyz in zip(mymol.atomnos,mymol.atomcoords[-1]):
            xyz.insert(0,t.element[sym]) 
        print "atomcoords and added elements:"
        for item in mymol.atomcoords[-1]:
            print ' '.join(map(str,item)) 
        if debug==True:
            import submitter
            if hasattr(mymol,'optdone'):
                if mymol.optdone==False:
                    print "Optimizations not converged!"
                elif mymol.optdone==True:
                    print "Optimization is converged!"
            fid = open(path[:-3]+'com','r') #change .log in .com extension and read input file
            multcharge = re.compile('^\-?[01]\s[12]') #a regex for the mult charge line
            newfile=[]
            once=0 #only find that line once
            for line in fid: #copy file exept for the zmat found in the inputfile
                if multcharge.match(line) and once==0: #when found 
                    once+=1
                    newfile.append(line) #the line with the match itself has to be included in the newfile
                    while True: 
                        line= next(fid) #take al new lines
                        if line=='\n': #end of zmat
                            #now instead of this zmat that is now completely skipped place in newfile
                            #the last coordinates of the crashed run
                            newfile.extend([' '.join(map(str,item))+'\n' for item in mymol.atomcoords[-1]])
                            newfile.extend(['\n'])
                            break
                else:
                    newfile.append(line) #copy that line because it is not the zmat found in the inputfile
            print "="*20
            #for line in newfile: print line,
            #print newfile
            open(path[:-4]+'zzz.com','w').writelines(newfile)
            print "newfile written in: ", path[:-4] + 'zzz.com'

            #---- preparation for submit command ---
            splitpath = path.split('/')
            filename = splitpath[-1]
            folder = '/'.join(splitpath[:-1])
            filenamesplit = filename[:-4].split('_')
            identify= filenamesplit[0]+'_'
            index = '_'.join(filenamesplit[1:])+'zzz.com'
            print "folder", folder
            print "index:", index
            print "identi", identify
            #----- keywords constructed so:
            submitter.submit(folder,index,identify)
            return True
    return False

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
        pprint(mymol.scfenergies[-3:])
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
     
  
