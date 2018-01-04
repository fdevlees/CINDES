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
    mols_tocal = normaltermination( mols_tocal, fileparameters)

    # 2. get a list of properties that need to be extracted for each molecule
    uni_props_set = fileparameters['props']

    # 3. obtain data for each molecule
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
    # for every jobfile do a cclib extraction. faking the separate jobs as if it were single files
    jobslines = open(filename).read().split('termination')[:-1]
    from cStringIO import StringIO
    jobfiles = map(StringIO, jobslines)
    datadict = dict()
    for jobfile, job in zip(jobfiles, jobs):

        if program=='gaussian':
            from CINDES4.cclib.parser.gaussianparser import Gaussian
            job_data = Gaussian(jobfile).parse()
        elif program=='orca':
            from CINDES4.cclib.parser.orcaparser import ORCA
            job_data = ORCA(jobfile).parse()
        else:
            raise SystemExit('not implemented')


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
    print "datadict:", datadict
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
    if any(i in to_read_props for i in ['IP','omega','stab']):
        results['IP']= (datadict['eIP']-datadict['e0'])*27.2113838
    if any(i in to_read_props for i in ['EA','omega','stab']):
        results['EA']= (datadict['e0']-datadict['eEA'])*27.2113838
    if any(i in to_read_props for i in ['omega', 'stab']):
        results['omega'] = ( ( results['IP'] + results['EA'] )**2 ) / ( 8 * ( results['IP'] - results['EA'] ))
    if 'stab' in to_read_props:
        results = calculate_stab(results, EAHs)

    print "results:", results
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
            paths.append(fileparameters['path'] + '/' + mol.index + '/' + fileparameters['identify'] + mol.index + '_' + str(pos) + '.log')
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
