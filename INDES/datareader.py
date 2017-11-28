from pprint import pprint
#from writings import log_io
from CINDES4.utils.writings import log_io, print_title, sprint
import re
import numpy
import time
import logging
logging.basicConfig(level=logging.DEBUG)

import construction

corresp = {2:46,6:18,7:42,9:34,11:22,12:30}

def round_sig( x, sig=8):
    from math import log10, floor
    try:
        return round(x, sig-int(floor(log10(abs(x))))-1)
    except ValueError:
        if not x==0.0: print "ValueError:", x
        return x

def rm_duplicates(seq, nsig=8):
    seen = set()
    seen_add = seen.add
    new_seq = []
    if True: # try to round to numerical precision errors:
        print "    the sequence(scfenergies?) is rounded to max 10 significant digits"
        seq = [ round_sig( item, sig=10 ) for item in seq ]
    for x in seq:
        if x in seen:
            print "    multiples found in sequence. name probably scfenergies! | value: ", x
        else:
            new_seq.append(x)
            seen_add(x)
    return new_seq

class prettyfloat(float):
    def __repr__(self):
        return "%-0.4f" % self

class Logfile():

    def __init__(self,name, afile=True):
        self.data = []

        # here a trick. later i will remove the afile part but for now i want to keep 
        # both functionalities. so in the new version. afile=False and self.name is just a longstring.
        if afile:
            self.fid = open(name) # + '.log')
        else:
            self.fid = self.name.split('\n')
        self.name = name
        return

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


def extract_stab(file1 , molecule, fileparameters):
    index = molecule.index
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

    propsA =gausread(file1,'stabA')
    # propsA = { Eopt:..., E0:..., IP:..., EA:... }

    EAHs=[] #all EAHs from all different positions in here
    Npos=[] #positions with a nitrogen in here
    for pos in fileparameters['positions']: #extract al AH energies and take the lowest
        file2= fileparameters['path'] + '/' + index + '/' + fileparameters['identify'] + index + '_' + str(pos) + '.log'

        EAH= gausread(file2,'energy',multiplejobs=0)['energy']

        #---- HERE THE electronegativity part of the stab a bit tricky
        #confje = construction.indtocon(index) # change index to conf list using the construction module
        confje = index.split('_')
        #corresp = {2:46,6:18,7:42,9:34,11:22,12:30} # map the alpha positions to methyl indices 
        try:
            siteindex = fileparameters['corresp'][pos] #find for each position the methyl index
        except KeyError:
            print "position of H atom is not a possible site. Therefore the program assumes H is attached to a nitrogen atom!"
            print "electronegativity correction for nitrogen: ", chi_term
            Npos.append(pos)
        else:
            if siteindex in fileparameters['sites']:# look if that index is used as a site
                if confje[ fileparameters['sites'].index(siteindex) ]==['N']:
                    print "electronegativity correction for nitrogen: ", chi_term
                    Npos.append(pos)
        #-----
        EAHs.append([EAH,pos])
    print "EAHs:"
    pprint(EAHs)
    print "min EAHs:", min(EAHs)
    print "Npos:",Npos
    #gasconstant = 8.3144621
    E_ah =min(EAHs)

    I = propsA['IP']
    A = propsA['EA']
    propsA['omega'] = ( ((I+A)**2 )/(8*(I-A)) )*eV #in eV
    Domega = propsA['omega'] - 2
    propsA['BDE_ah']  = ( propsA['Eopt'] + H_h - E_ah[0])*kJmol + avtc #avtc is AVerage Thermal Correction. 
    #----
    propsA.update(gausread(file1,'rdv',2) )
    #----
    if E_ah[1] in Npos:
        propsA['stab'] = propsA['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    else:
        propsA['stab'] = propsA['BDE_ah'] - stab_h - bde_a * Domega * Dw_h
    propsA['H_pos'] = E_ah[1]
    return propsA

def get_paths( mols, fileparameters):
    ''' get all paths that need to be examined later '''
    files=[]
    path = fileparameters['path']
    for molecule in mols:
        files.append(path + '/' + fileparameters['identify'] + molecule.index + '.log')

        if fileparameters['stab']==1:
            for pos in fileparameters['positions']: #extract al AH energies and take the lowest
                files.append(path + '/' + molecule.index + '/' + fileparameters['identify'] + molecule.index + '_' + str(pos) + '.log')

    return files


@log_io()
def datareader( mols_tocal, fileparameters):
    # 1. get all paths
    paths = get_paths( mols_tocal, fileparameters)

    # 2. test normal termination
    normaltermination( paths, fileparameters['debug'] )

    # 3. get a list of properties that need to be extracted for each molecule
    # THIS IS ALREADY DONE AT INPUTREADER > run.props
    uni_props_set = fileparameters['props']

    # 4. obtain data for each molecule
    for molecule in mols_tocal:
        print "><"*10, molecule, "><"*10
        # make a copy of props_dict
        props_set = uni_props_set.copy()

        file1 = fileparameters['path'] + '/' + fileparameters['identify'] + molecule.index + '.log'

        # start by looking if stab is one of the crucial properties because it contains many others
        if 'stab' in props_set:
            X_stab_props = extract_stab( file1, molecule, fileparameters )
            # expect to get something like: { 'stab': value, 'I':..., 'A':...,'omega':...,'RDV'....}

            # fill props_dict
            #props_dict.update(X_stab_props)
            molecule.props.update(X_stab_props)

        # check which properties are still necessary to obtain:
        #to_read_props = [ key for key, value in props_dict.iteritems() if value==None ]
        # new: test which in fileparameters['props'] but not in molecule.props.viewkeys()
        # NB: - is here a set operator! returns a set!
        to_read_props = props_set - molecule.props.viewkeys()
        print "to read props:", to_read_props

        # extract them
        # assume all properties can be easily obtained by gausread
        if to_read_props:
            if fileparameters['jobs']: # new JSON / cclib style
                readings = new_style_reader( file1, to_read_props, fileparameters )
            else: # old pickle. Logfile-gausread style
                readings = gausread( file1, to_read_props, multiplejobs=fileparameters['multiplejobs'])
            molecule.props.update( readings )

        molecule.predicted = False

    return mols_tocal

def new_style_reader( file1, to_read_props, fileparameters ):
    # for every jobfile do a cclib extraction. faking the separate jobs as if it were single files
    jobslines = open(file1).read().split('termination')[:-1]
    from cStringIO import StringIO
    jobfiles = map(StringIO, jobslines)
    datadict = dict()
    for jobfile, job in zip(jobfiles, fileparameters['jobs']):
        from CINDES4.cclib.parser.gaussianparser import Gaussian
        job_data = Gaussian(jobfile).parse()
        for inf in job['info']:
            if inf[0]=='e': # so it concerns an energy!:
                print inf, "scfenergies:", job_data.scfenergies
                datadict[inf]=job_data.scfenergies[-1]
            elif inf in ['homo','lumo']:
                datadict['homo']=job_data.moenergies[-1][job_data.homos[0]]
                datadict['lumo']=job_data.moenergies[-1][job_data.homos[0]+1]
            elif inf=='dipole': datadict[inf]=job_data.dipole
            elif inf=='polar':
                datadict['polar'] = sum([job_data.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
            elif inf=='mw':
                datadict['mw'] = float( sum( job_data.atomnos) )
            else:
                print "value not recognized:", inf

    print "datadict:", datadict

    # so now datadict should have all energy keys + homo/lumo + dipole
    # but not yet omega/solv/gap so:
    results=datadict.copy()
    if 'gap' in to_read_props:
        results['gap']=datadict['lumo']-datadict['homo']
    if 'solv' in to_read_props:
        results['solv']= (datadict['e1_solv']-datadict['e0_solv'])*627.5
    if any(i in to_read_props for i in ['IP','omega']):
        results['IP']= (datadict['eIP']-datadict['e0'])
        print "results:IP", results['IP']
    if any(i in to_read_props for i in ['EA','omega']):
        results['EA']= (datadict['e0']-datadict['eEA'])
    if 'omega' in to_read_props:
        results['omega'] = ( ( results['IP'] + results['EA'] )**2 ) / ( 8 * ( results['IP'] - results['EA'] ))

    print "results:", results
    return results

def gausread(filename,props,multiplejobs=1,rdvindex=1, afile=True):
    ''' props is a set of props to extract '''
    if afile:
        mymol = Logfile(filename)
    else:
        mymol = Logfile(filename, afile=False)
    results = {}

    # extract
    if any( prop in ['natom','natoms'] for prop in props ):
        mymol.extract(coords=1)
    else:
        mymol.extract()

    # check opt 
    if not hasattr(mymol,'optdone'):
        print "program did not do optimization or crashed"

    # set props
    if 'polar' in props:
        if hasattr(mymol,'polex'):
            print "polar exact densities:"
            pprint(mymol.polex)
            pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
            # next line calculates in one line the Radical delocalisation value
            # RDV = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
            print "polarisability value is: ", pola
            results['polar']=pola
        else:
            print "NO POLAR DATA FOUND IN FILE!"
    if 'stabA' in props:
        pprint(mymol.scfenergies[-5:])
        #try to remove duplicates from list
        if True:
            mymol.scfenergies = rm_duplicates(mymol.scfenergies)
        results['Eopt'] = mymol.scfenergies[-4]
        results['E0']   = mymol.scfenergies[-3]
        results['IP']   = mymol.scfenergies[-2] - results['E0']
        results['EA']   = results['E0'] - mymol.scfenergies[-1]
    if 'omega' in props:
        pprint(mymol.scfenergies[-5:])

        #try to remove duplicates from list
        if True:
            mymol.scfenergies = rm_duplicates(mymol.scfenergies)

        if multiplejobs > 1:
            print 'i am here'
            scfenergies = mymol.scfenergies[:-(multiplejobs-1)]
        else:
            scfenergies = mymol.scfenergies[:]

        results['Eopt'] = scfenergies[-4]
        results['E0']   = scfenergies[-3]
        results['IP']   = scfenergies[-2] - results['E0']
        results['EA']   = results['E0'] - scfenergies[-1]
        #omega = lambda I,A: ( (I+A)**2 ) / ( 8 * (I-A) )
        #results['omega']= omega(results['IP'], results['EA'])
        results['omega'] = ( ( results['IP'] + results['EA'] )**2 ) / ( 8 * ( results['IP'] - results['EA'] ) ) * 27.2113838
        if results['omega']<0.0:
            print "    FAULTY VALUE FOR ELECTROPHILICITY: cannot be a negative value:"
            print "    VALUE is set to None"
            results['omega']=None
    if 'energy' in props:
        ESCFs = mymol.scfenergies
        results['energy'] = ESCFs[-(multiplejobs+1)]
    if any( prop in ['homo','lumo'] for prop in props):
        from CINDES4.cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        results['homo'] = myfile.moenergies[0][myfile.homos[0]] #NOT IMPLEMENTED NEED CCLIB
        results['lumo'] = myfile.moenergies[0][myfile.homos[0]+1] #NOT IMPLEMENTED NEED CCLIB
    if 'gap' in props:
        from CINDES4.cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        results['homo'] = myfile.moenergies[0][myfile.homos[0]]
        results['lumo'] = myfile.moenergies[0][myfile.homos[0]+1]
        results['gap']  = results['lumo'] - results['homo']
    if 'ip' in props:
        ESCFs = mymol.scfenergies

        #try to remove duplicates from list
        if True:
            mymol.scfenergies = rm_duplicates(mymol.scfenergies)

        results['energy']  = ESCFs[-(multiplejobs+1)]
        results['ecation'] = ESCFs[-multiplejobs]
        ip                 = results['ecation'] - results['energy']

        # now i want to test if it is not too small.
        if True:
            if ip < 0.001:
                print "Ionization Potentential smaller than expected range!."
                print "Program will take one scfenergy earlier!"
                results['energy']  = ESCFs[-(multiplejobs+2)]
                ip                 = results['ecation'] - results['energy']

        results['ip']      = ip
    if 'ea' in props:
        ESCFs = mymol.scfenergies

        #try to remove duplicates from list
        if True:
            mymol.scfenergies = rm_duplicates(mymol.scfenergies)

        results['energy'] = ESCFs[-(1+multiplejobs)]
        results['eanion'] = ESCFs[-1]
        results['ea']     = results['energy'] - results['eanion']
    if 'natoms' in props:
        results['natoms'] = len(mymol.atomnos)
    if 'volume' in props:
        results['volume'] = mymol.volume
    if 'poldens' in props:
        pola = sum([ mymol.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
        results['polar']   = pola
        results['poldens'] = pola/mymol.volume
    if 'rdv' in props:
        print "Mulliken spin densities:"
        pprint(mymol.spindensities[rdvindex])
        spiden= mymol.spindensities[rdvindex]
        # next line calculates in one line the Radical delocalisation value
        results['rdv'] = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
    if 'dipole' in props:
        results['dipole'] = mymol.dipole
    if 'solv' in props:
        ESCFs = mymol.scfenergies
        print "ESCFs:", ESCFs
        results['e0_solv'] = ESCFs[-3]
        results['e1_solv'] = ESCFs[-2]
        results['solv']    = - ( results['e0_solv'] - results['e1_solv'] ) * 627.5
    if 'mw' in props:
        from CINDES4.cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        results['mw'] = float( sum( myfile.atomnos) )

    # assert that all props are filled
    #print 'results:', results
    #print "props:", props
    #if not 'stabA' in props:
    #    assert all( prop in results for prop in props), 'not all properties calculated '

    return results

def normaltermination(filepaths,debug=False):
    #-----
    def termination(filepath):
       with open(filepath,'r') as fid:
           text = fid.readlines()[-3:]
           if re.search('Normal termination',''.join(text)):
               fid.close()
               return 1
           #elif re.search('IGNORE',''.join(text)):
           #    fid.close()
           #    return 2
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
    mymol=Logfile(path) #read outputfile
    mymol.extract(coords=1) #extract file with also the coordinates
    from CINDES4.utils import utils
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
                            #newfile.extend([' '.join( map("{12.6f}".format, item))+'\n' for item in mymol.atomcoords[-1]])
                            xyz = mymol.atomcoords[-1]
                            xyz_f = [ item[0] + ' '.join( map( "{:12.6f}".format, item[1:])) + '\n' for item in xyz ]
                            #xyz_f = [ item[0] + ' '.join(
                            #                              map(
                            #                                   str, item[1:]
                            #                                 )
                            #                            ) + '\n' for item in xyz ]
                            print xyz_f
                            newfile.extend(xyz_f)
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
    try:
        index = sys.argv[2]
    except IndexError:
        index = 1
    mymol = Logfile(filename)
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
     
    if True:
        from CINDES4.cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        HOMO = myfile.myhomos[index]
        Ehomo= myfile.mymos[index]['alpha'][0][HOMO]
        Elumo= myfile.mymos[index]['alpha'][0][HOMO+1]

        Egap = Elumo - Ehomo
        print "E-HOMO :", Ehomo, Ehomo/27.2113838
        print "E-LUMO :", Elumo, Elumo/27.2113838
        print "BANDGAP:", Egap
        print Ehomo, Elumo, Egap
    if True:
        I = mymol.scfenergies[-1] - mymol.scfenergies[-2]
        #A = mymol.scfenergies[-1] - mymol.scfenergies[-3]
        print "I:", I, I*27.2113838
        #print "A:", A 
