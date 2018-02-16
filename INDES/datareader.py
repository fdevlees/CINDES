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

def extract_eahs(molecule, fileparameters):
    index = molecule.index
    EAHs=dict() #all EAHs from all different positions in here
    Npos=[] #positions with a nitrogen in here
    for pos in fileparameters['positions']: #extract al AH energies and take the lowest
        print "pos:", pos,
        file2= fileparameters['path'] + '/' + index + '/' + fileparameters['identify'] + index + '_' + str(pos) + '.log'

        datadict_file2 = read_file(file2, fileparameters['stabjobs'])
        # returns something like: '{'e':638.8, 'eAH':392.389 }

        #---- a bit tricky: get the pos-positions that correspond to a nitrogen-X (X=H,CH3) bond.
        confje = index.split('_')
        N=False #set initially to False
        try:
            siteindex = fileparameters['corresp'][pos] #find for each position the methyl index
        except KeyError:
            # there is theoretically a possibility that the position to add a A-X, (X=H,CH3) group is not an active or passive site
            # in this case the fileparameters['corresp'] does not contain the pos this is only possible when te possible reactive center
            # has no hydrogen for the case of phenenalenyl. for thiadiazinyl this is automatically an sp2 nitrogen position
            print "position of H atom is not a possible site. Therefore the program assumes H is attached to a nitrogen atom!"
            N=True
        else:
            # this is only executed when no Error is raised!
            if siteindex in fileparameters['sites']:# look if that index is used as a site
                if any(confje[ fileparameters['sites'].index(siteindex) ]==n for n in [['N'],'N']):
                    print "    there is a nitrogen on this position!    "
                    N=True

        if N:
            Npos.append(pos)
        EAHs[pos]={'eAH':datadict_file2['eAH'], 'N':N}

    print "EAHs:"
    pprint(EAHs)
    print "Npos:",Npos
    return EAHs

def calculate_stab(results, EAHs):
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
    #gasconstant = 8.3144621
    #----- end of parameters

    minpos = min(EAHs, key=lambda x:EAHs[x]['eAH'])
    E_ah = EAHs[minpos]['eAH']

    Domega = results['omega'] - 2.
    print "Domega:", Domega
    results['BDE_ah']  = ( results['eA'] + H_h - E_ah)*kJmol + avtc #avtc is AVerage Thermal Correction. 
    #if E_ah[1] in Npos:
    if EAHs[minpos]['N']:
        print "electronegativity correction for nitrogen is used"
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    else:
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h
    results['H_pos'] = minpos
    results['stab'] = stab
    return results

@log_io()
def datareader( mols_tocal, fileparameters):
    # 1. test normal termination
    program=fileparameters['program']
    if program=='gaussian':
        import gaussian as program
    elif program=='orca':
        import orca as program
    elif program=='nwchem':
        import nwchem as program
    else:
        raise SystemExit('not implemented')
    mols_tocal = program.normaltermination( mols_tocal, fileparameters)

    # 2. get a list of properties that need to be extracted for each molecule
    uni_props_set = fileparameters['props']

    # 3. obtain data for each molecule
    for molecule in mols_tocal:
        print "><"*15, molecule
        # make a copy of props_dict
        props_set = uni_props_set.copy()
        file1 = fileparameters['path'] + '/' + fileparameters['identify'] + molecule.index + '.log'

        # start by looking if stab is one of the crucial properties because it contains many others
        if 'stab' in props_set:
            # read EAHs for the stabfiles
            EAHs = extract_eahs(molecule, fileparameters)
            # expect to get something like: { 2:EAH2, 4:EAH4, 12:EAH12 }
        elif 'aromaticity' in props_set:
            pass
        else:
            EAHs = None

        # extract them
        readings = new_style_reader( file1, props_set, fileparameters, EAHs )
        molecule.props.update( readings )

        molecule.predicted = False

    return mols_tocal

def read_file(filename, jobs, program='gaussian'):
    if program=='gaussian':
        from CINDES4.cclib.parser.gaussianparser import Gaussian as Logfile
        key='termination'
        jobslines = open(filename).read().split(key)[:-1]
    elif program=='orca':
        from CINDES4.cclib.parser.orcaparser import ORCA as Logfile
        raise SystemExit('not implemented')
    elif program=='nwchem':
        from CINDES4.cclib.parser.nwchemparser import NWChem as Logfile
        key='NWChem Input Module'
        jobslines = open(filename).read().split(key)[1:-1]
    else:
        raise SystemExit('not implemented')
    # for every jobfile do a cclib extraction. faking the separate jobs as if it were single files
    print "njobs:", len(jobslines)
    from cStringIO import StringIO
    jobfiles = map(StringIO, jobslines)
    datadict = dict()
    for jobfile, job in zip(jobfiles, jobs):
        job_data = Logfile(jobfile).parse()
        #print "job_data:", job_data

        for inf in job['info']:
            if inf=='_':continue
            elif inf[0]=='e': # so it concerns an energy!:
                datadict[inf]=job_data.scfenergies[-1]/27.21138505 # this value is used in cclib
            elif inf in ['homo','lumo']:
                datadict['homo']=job_data.moenergies[-1][job_data.homos[0]]
                datadict['lumo']=job_data.moenergies[-1][job_data.homos[0]+1]
            elif inf=='dipole': datadict[inf]=job_data.dipole
            elif inf=='polar':
                datadict['polar'] = sum([job_data.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
            elif inf=='mw':
                datadict['mw'] = float( sum( job_data.atomnos) )
            elif inf=='rdv':
                print "job_data.atomcharges:", job_data.atomcharges
                spiden=job_data.atomcharges['natural']
                datadict['rdv'] = sum([ float(item[2])**2 for item in spiden if abs(item[2])>0.05 ])
                print datadict['rdv']
                print job_data.atomcharges
                try:
                    print job_data.atomspins
                except AttributeError:
                    pass
                raise SystemExit('rdv not tested yet')
            else:
                print "value not recognized:", inf
    #print "datadict:", datadict
    return datadict

def new_style_reader( file1, to_read_props, fileparameters, EAHs=None ):
    datadict = read_file(file1, fileparameters['jobs'], program=fileparameters['program'])

    # so now datadict should have all energy keys + homo/lumo + dipole
    # but not yet omega/solv/gap so:
    results=datadict.copy()
    if 'gap' in to_read_props:
        results['gap']=datadict['lumo']-datadict['homo']
    if 'solv' in to_read_props:
        results['solv']= (datadict['e1_solv']-datadict['e0_solv'])*627.5
    if any(i in to_read_props for i in ['ip','omega','stab']):
        results['ip']= (datadict['eIP']-datadict['e0'])*27.2113838
    if any(i in to_read_props for i in ['ea','omega','stab']):
        results['ea']= (datadict['e0']-datadict['eEA'])*27.2113838
    if any(i in to_read_props for i in ['omega', 'stab']):
        results['omega'] = ( ( results['ip'] + results['ea'] )**2 ) / ( 8 * ( results['ip'] - results['ea'] ))
    if 'stab' in to_read_props:
        results = calculate_stab(results, EAHs)

    print "results:", results
    return results

if __name__ == "__main__":
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
