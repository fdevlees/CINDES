from pprint import pprint
#from writings import log_io
from CINDES.utils.writings import log_io, print_title, sprint
import re
import numpy
import time
import logging
logging.basicConfig(level=logging.DEBUG)

import construction

corresp = {2:46,6:18,7:42,9:34,11:22,12:30}

def setEAHs(molecule):
    EAHs=dict()
    Npos=[] #positions with a nitrogen in here
    for job in molecule.jobs:
        if hasattr(job, 'pos'):
            EAHs[job.pos]={'eAH':molecule.props.pop('eAH_P{}'.format(str(job.pos))), 'N':job.N}
    if EAHs:
        molecule.props['EAHs']=EAHs
    return

def calculate_stab(results, molecule):
    #print "molecule:", molecule
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
    EAHs=results['EAHs']
    minpos = min(EAHs, key=lambda x:EAHs[x]['eAH'])
    E_ah = EAHs[minpos]['eAH']

    Domega = results['omega'] - 2.
    results['BDE_ah']  = ( results['eA'] + H_h - E_ah)*kJmol + avtc #avtc is AVerage Thermal Correction. 
    if EAHs[minpos]['N']:
        print "electronegativity correction for nitrogen is used"
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    else:
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h
    results['H_pos'] = minpos
    results['stab'] = stab
    return results


def normaltermination(mols, **kwargs):
    """ kwargs ignore / debug this new normal termination runs more parallel """
    # 1. test normal termination and submit errorjob
    notready=''
    for i, mol in enumerate(mols):
        for j, job in enumerate(mol.jobs):
            job.normaltermination(**kwargs)
            if not job.IsReady:
                notready += "{}.{}: {}\n".format(i,j, job.name)
    print "jobs not ready:\n", notready

    # 2. test normal termination and errorjob are ready or molecule is ignored
    extratime = 0
    timestep1 = 10
    timestep2 = 300
    while True:
        # CHECK READY:
        print "not ready:",
        for i, mol in enumerate(mols):
            if mol.ignoremol or mol.IsReady: continue
            for j, job in enumerate(mol.jobs):
                if job.IsReady or job.ignorejob: continue
                else:
                    print "{}.{}".format(i,j),
                job.ready(**kwargs)
                if job.ignorejob: mol.discard()

            if all(job.IsReady for job in mol.jobs):
                mol.IsReady=True
        print
        if all(mol.IsReady for mol in mols):
            break

        # WAIT:
        time.sleep(timestep1) # wait 5 minudtes
        extratime += timestep1
        # after some time use larger timesteps
        if extratime>=timestep2:timestep1=timestep2
        print "extra waittime/h:", "{:.2f}".format(round(extratime/3600.,2)), "||",

    # 3. return mols that are not ignored:
    mols_toread = filter(lambda mol:not mol.ignoremol, mols)
    return mols_toread

@log_io()
def datareader( mols, run):
    # 1. test normal termination
    mols_toread = normaltermination(mols, debug=run.debug, ignore=run.ignore)

    # 2. obtain data for each molecule
    for molecule in mols_toread:
        print "><"*15, molecule,
        for job in molecule.jobs:
            readings = read_file(job)

            if hasattr(job,'pos'):
                for old_key in readings.keys(): #the .keys is very important here. iterkeys or for just readings do not work!
                    readings["{}_P{}".format(old_key, str(job.pos))] = readings.pop(old_key)

            molecule.props.update(readings)
        molecule.predicted = False

        # only relevant for stab calculations
        setEAHs(molecule)

    return

def read_file(Job):
    # 1. look to which calc the logfile belongs when there were simultaneous calculations:
    program=Job.calc['program']
    jobs=Job.calc['jobs']
    filename=Job.logpath

    # 2. split logfile in different jobs
    if program=='gaussian':
        from CINDES.cclib.parser.gaussianparser import Gaussian as Logfile
        key='termination'
        # if keyword freq in line than there is an extra internal job!
        jobslines_v1 = open(filename).read().split(key)[:-1]
        jobslines = []
        for joblines_v1 in jobslines_v1:
            if 'roceeding to internal job step number' in joblines_v1.split('\n',2)[1]:
                print "freq job appended to main job"
                jobslines[-1]+= joblines_v1
            else:
                jobslines.append(joblines_v1)
    elif program=='orca':
        from CINDES.cclib.parser.orcaparser import ORCA as Logfile
        raise SystemExit('ORCA interface not implemented')
    elif program=='nwchem':
        from CINDES.cclib.parser.nwchemparser import NWChem as Logfile
        key='NWChem Input Module'
        splitted = open(filename).read().split(key)
        jobslines = splitted[1:-1]
    else:
        raise SystemExit('not implemented')
    print "njobs:", len(jobslines),
    if len(jobslines)==0:
        print "no jobs in logfile!"
        raise SystemExit('should not occur here')


    # 3. handle every subjob as a different logfile and read the needed job['info'] from it
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
            elif inf[0]=='g' and not inf=='gap':
                # note that scfenergies are given in eV by cclib but free energy in hartree
                datadict[inf]=job_data.freeenergy
                print "free energy found:", job_data.freeenergy
            elif inf in ['homo','lumo']:
                datadict['homo']=job_data.moenergies[-1][job_data.homos[0]]
                datadict['lumo']=job_data.moenergies[-1][job_data.homos[0]+1]
            elif inf=='dipole': datadict[inf]=job_data.dipole
            elif inf=='polar':
                datadict['polar'] = sum([job_data.polex[i]/3 for i in [0,2,5]]) # = 1/3*(axx+ayy+azz)
            elif inf=='mw':
                datadict['mw'] = float( sum( job_data.atomnos) )
            elif inf=='rdv':
                spiden = map(lambda x:x[0]-x[1], zip(job_data.npaa, job_data.npab))
                #spiden=job_data.atomcharges['natural']
                datadict['rdv'] = sum([ item**2 for item in spiden if abs(item)>0.05 ])
                #try:
                #    print job_data.atomspins
                #except AttributeError:
                #    pass
            elif inf=='spindensities':
                spiden = map(lambda x:round(x[0]-x[1], 8), zip(job_data.npaa, job_data.npab))
                datadict['spindensities']=spiden
            elif any(prop in inf for prop in ['pcharges','partialcharges']):
                print "partial charges", job_data.atomcharges
                round8 = lambda x:round(float(x), 8)
                try:
                    pcharges = zip(map(int, job_data.atomnos), map(round8, job_data.atomcharges['mulliken']))
                except KeyError:
                    pcharges = zip(map(int, job_data.atomnos), map(round8, job_data.atomcharges['natural']))
                datadict[inf]=pcharges
            else:
                print "value not recognized:", inf
    print


    return datadict

def set_combined_variables(mol, to_read_props):
    # so now datadict should have all energy keys + homo/lumo + dipole
    # but not yet omega/solv/gap so:
    results=mol.props
    if 'gap' in to_read_props:
        results['gap']   = results['lumo']-results['homo']
    if 'solv' in to_read_props:
        results['solv']  = (results['e1_solv']-results['e0_solv'])*627.5
    if any(i in to_read_props for i in ['ip','omega','stab']):
        results['ip']    = (results['eIP']-results['e0'])*27.2113838
    if any(i in to_read_props for i in ['ea','omega','stab']):
        results['ea']    = (results['e0']-results['eEA'])*27.2113838
    if any(i in to_read_props for i in ['omega', 'stab']):
        results['omega'] = (( results['ip'] + results['ea'] )**2 ) / ( 8 * ( results['ip'] - results['ea'] ))
    if 'stab' in to_read_props:
        results = calculate_stab(results, mol)
    if any(i in to_read_props for i in ['ipfukui','radfukui']):
        print "results:", results
        results['ipfukui'] = [ -( q_ip[1] - q_0[1] ) for q_ip, q_0 in zip(results['pchargesIP'], results['pcharges0'])]
    if any(i in to_read_props for i in ['eafukui','radfukui']):
        results['eafukui'] = [ -( q_0[1] - q_ea[1] ) for q_0, q_ea in zip(results['pcharges0'], results['pchargesEA'])]
    if 'radfukui' in to_read_props:
        results['radfukui'] = [ .5*(ipf + eaf) for ipf, eaf in zip(results['ipfukui'], results['eafukui'])]

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
        from CINDES.cclib.parser import ccopen
        myfile=ccopen(filename).parse()
        HOMO = myfile.myhomos[index]
        Ehomo= myfile.mymos[index]['alpha'][0][HOMO]
        Elumo= myfile.mymos[index]['alpha'][0][HOMO+1]

        Egap = Elumo - Ehomo
        print "E-HOMO :", Ehomo, Ehomo/27.2113838
        print "E-LUMO :", Elumo, Elumo/27.2113838
        print "BANDGAP:", Egap
        print Ehomo, Elumo, Egap
