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
    EAHs=[] #all EAHs from all different positions in here
    Npos=[] #positions with a nitrogen in here
    for pos in fileparameters['positions']: #extract al AH energies and take the lowest
        file2= fileparameters['path'] + '/' + index + '/' + fileparameters['identify'] + index + '_' + str(pos) + '.log'

        datadict_file2 = read_file(file2, fileparameters['stabjobs'])
        # returns something like: '{'e':638.8, 'eAH':392.389 }

        #---- HERE THE electronegativity part of the stab a bit tricky
        confje = index.split('_')
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
    E_ah =min(EAHs)
    I = results['IP']
    A = results['EA']
    #propsA['omega'] = ( ((I+A)**2 )/(8*(I-A)) )*eV #in eV
    Domega = results['omega'] - 2
    results['BDE_ah']  = ( results['Eopt'] + H_h - E_ah[0])*kJmol + avtc #avtc is AVerage Thermal Correction. 
    if E_ah[1] in Npos:
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    else:
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h
    results['H_pos'] = E_ah[1]
    results['stab'] = stab
    return results

@log_io()
def datareader( mols_tocal, fileparameters):
    # 1. get all paths

    # 2. test normal termination
    mols_tocal = normaltermination( mols_tocal, fileparameters)

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
            # read EAHs for the stabfiles
            EAHs = extract_eahs(molecule, fileparameters)
            # expect to get something like: { 2:EAH2, 4:EAH4, 12:EAH12 }
        else:
            EAHs = None

        # extract them
        if to_read_props:
            readings = new_style_reader( file1, to_read_props, fileparameters, EAHs )
            molecule.props.update( readings )

        molecule.predicted = False

    return mols_tocal

def read_file(filename, jobs):
    # for every jobfile do a cclib extraction. faking the separate jobs as if it were single files
    jobslines = open(filename).read().split('termination')[:-1]
    from cStringIO import StringIO
    jobfiles = map(StringIO, jobslines)
    datadict = dict()
    for jobfile, job in zip(jobfiles, jobs):
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
    return datadict

def new_style_reader( file1, to_read_props, fileparameters, EAHs=None ):
    datadict = read_file(file1, fileparameters['jobs'])

    # so now datadict should have all energy keys + homo/lumo + dipole
    # but not yet omega/solv/gap so:
    results=datadict.copy()
    if 'gap' in to_read_props:
        results['gap']=datadict['lumo']-datadict['homo']
    if 'solv' in to_read_props:
        results['solv']= (datadict['e1_solv']-datadict['e0_solv'])*627.5
    if any(i in to_read_props for i in ['IP','omega','stab']):
        results['IP']= (datadict['eIP']-datadict['e0'])
        print "results:IP", results['IP']
    if any(i in to_read_props for i in ['EA','omega','stab']):
        results['EA']= (datadict['e0']-datadict['eEA'])
    if any(i in to_read_props for i in ['omega', 'stab']):
        results['omega'] = ( ( results['IP'] + results['EA'] )**2 ) / ( 8 * ( results['IP'] - results['EA'] ))
    if 'stab' in to_read_props:
        # 
        results = calc_stab(results, EAHs)

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

def get_paths( mols, fileparameters):
    ''' get all paths that need to be examined later '''
    files=[]
    #path = fileparameters['path']
    for molecule in mols:
        #files.append(path + '/' + fileparameters['identify'] + molecule.index + '.log')
        #if fileparameters['stab']==1:
        #    for pos in fileparameters['positions']: #extract al AH energies and take the lowest
        #        files.append(path + '/' + molecule.index + '/' + fileparameters['identify'] + molecule.index + '_' + str(pos) + '.log')
        files.extend(get_molpaths(molecule, fileparameters))
    return files

def get_molpaths(mol, fileparameters):
    paths=[]
    paths.append(fileparameters['path'] + '/' + fileparameters['identify'] + mol.index + '.log')
    if fileparameters['stab']==1:
        for pos in fileparameters['positions']: #extract al AH energies and take the lowest
            paths.append(fileparameters['path'] + '/' + molecule.index + '/' + fileparameters['identify'] + molecule.index + '_' + str(pos) + '.log')
    return paths

def normaltermination(mols_tocal, fileparameters):
    #-----
    def termination(filepath):
       with open(filepath,'r') as fid:
           text = fid.readlines()[-3:]
           if re.search('Normal termination',''.join(text)):
               fid.close()
               return 1
           elif re.search('IGNORE',''.join(text)):
               fid.close()
               return 2
    #-----
    filepaths = get_paths( mols_tocal, fileparameters)
    debug=fileparameters['debug']
    copyfilepaths = filepaths[:] #copy to be able to append to it while looping over it
    for path in copyfilepaths: #test all for information which jobs crashed
        if not termination(path)==1:
            print "Error termination:",path
            bnewfile = errortermination(path,debug)
            #if bnewfile and debug:
            #    filepaths.append(path[:-4]+'zzz.com')
    extratime = 0
    once = 0

    #--- new:
    mols_toread=[]
    for mol in mols_tocal:
        # get path belonging to this particular mol
        molpaths=get_molpaths(mol, fileparameters)
        ignoremol=False
        for path in molpaths: #test one by one waiting for normal termination
            while True:
                if termination(path)==1:
                    break
                elif termination(path)==2:
                    print "\n\n{0}\n             INGORED: {1} IGNORED!\n{0}\n".format("    --oOo--"*10, path)
                    ignoremol=True
                    break
                else:
                    print "no normal termination for: ",path
                #time.sleep(300) # wait 5 minudtes
                time.sleep(300) # wait 5 minudtes
                extratime += 300
                print "extra waittime/h:", extratime/3600, "||",
            # when I'm here this path has normal termination
            if ignoremol: break # this ignores the other paths belonging to this mol
        # when I'm here every molpath of this mol should have normal termination
        if not ignoremol:
            mols_toread.append(mol)
    # when I'm here every mol should have normal termination


    #--- old:
    #for path in filepaths: #test one by one waiting for normal termination
    #    while True:
    #       if termination(path)==1:
    #           break
    #       elif termination(path)==2:
    #           print "\n    {} IGNORED!\n".format(path)
    #           break
    #       else:
    #           print "no normal termination for: ",path
    #       #time.sleep(300) # wait 5 minudtes
    #       time.sleep(300) # wait 5 minudtes
    #       extratime += 300
    #       print "extra waittime/h:", extratime/3600, "||",
    #mols_toread = mols_tocal
    # ---
    print "mols_toread:", mols_toread
    return mols_toread

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
