#!/bin/env python

import re
import numpy
import time
import logging
import os#@@@david@@
from pprint import pprint

from CINDES.utils.writings import log_io, print_title, sprint
from .submitter import qsta
from . import construction

def round8(x): return round(float(x), 8)


def setEAHs(molecule):
    EAHs = dict()
    Npos = []  # positions with a nitrogen in here
    for job in molecule.jobs:
        if hasattr(job, 'pos'):
            eAH = molecule.props.pop('eAH_P{}'.format(str(job.pos)))
            EAHs[job.pos] = {
                'eAH': eAH,
                'Aatom': job.Aatom}

            tcAHkey = 'tchAH_P{}'.format(str(job.pos))
            if tcAHkey in molecule.props:
                tcAH = molecule.props.pop(tcAHkey)
                EAHs[job.pos]['tchAH'] = tcAH
                EAHs[job.pos]['enthalpy'] = eAH + tcAH
            else:
                print("no tch", end=' ')
                EAHs[job.pos]['enthalpy'] = eAH
    if EAHs:
        molecule.props['EAHs'] = EAHs
        print("EAHs:", EAHs)
    else:
        print('molecule has no EAHs!:', molecule)
    return


def calculate_stab(results, molecule):
    #print "molecule:", molecule
    # ---- some parameters needed
    bde_a = -12.68  # kJ/mol/eV^2
    bde_b = -218.1  # kJ/mol
    stab_h = 235.8  # kJ/mol
    Dw_h = 0.063  # eV

    # here the chi terms. The only relevant chi's are the ones higher than 3!
    # values come from Freija's radical stability scale paper!
    chi_h = 2.20
    chi_c = 2.60
    chi_n = 3.05
    chi_o = 3.50
    chi_f = 4.00
    chi_s = 2.60
    chi_cl= 3.15
    chi_br= 2.85
    H_h = -0.514457233  # a.u.
    E_h = -0.516817233  # a.u.
    kJmol = 2625.5
    eV = 27.2113838
    avtc = -28.1290706  # kJ/mol #average thermal correction for 5 random structures kJ/mol
    #gasconstant = 8.3144621
    # ----- end of parameters
    try:
        Domega = results['omega'] - 2.
    except KeyError:
        Domega = 0.0
        print("No electrophilicity term found so stab is calculated without omega term")

    EAHs = results['EAHs']

    once = False
    for i, EAH in EAHs.items():
        # 1. first set BDE for each AH molecule
        if 'tchA' in results:
            if not once:
                print("applying thermal corrections", results['tchA'], 'and', EAH['tchAH'])
                # note that if tch_AH is present it is already present in the enthalpy term
                once = True
            EAH['BDE_ah'] = (results['eA'] + results['tchA'] + H_h - EAH['enthalpy']) * kJmol
        else:
            EAH['BDE_ah'] = (results['eA'] + E_h - EAH['enthalpy']) * kJmol

        # 2. than set stab for each AH molecule
        if EAH['Aatom']==7:
            chi_term = bde_b * (chi_h - 3) * (chi_n - 3)
            print("+Domega term for N", end=' ')
            EAH['stab'] = EAH['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
        elif EAH['Aatom']==8:
            chi_term = bde_b * (chi_h - 3) * (chi_o - 3)
            print("+Domega term for O", end=' ')
            EAH['stab'] = EAH['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
        else:
            EAH['stab'] = EAH['BDE_ah'] - stab_h - bde_a * Domega * Dw_h

    # now finally set the final values having the highest BDE
    # get position of max EAH
    maxpos = max(EAHs, key=lambda x: EAHs[x]['BDE_ah'])
    #H_ah = EAHs[maxpos]['enthalpy']
    results['H_pos'] = maxpos
    results['stab'] = EAHs[maxpos]['stab']
    return results


def normaltermination(mols, run):
    """ kwargs ignore / debug this new normal termination runs more parallel """
    # 1. test normal termination and submit errorjob
    notready = ''
    for i, mol in enumerate(mols):
        for j, job in enumerate(mol.jobs):
            job.normaltermination(debug=run.debug)
            if not job.IsReady:
                notready += "{}.{}: {}\n".format(i, j, job.name)
                mol.IsReady = False
    if notready:
        print("jobs not ready:\n", notready)

    # 2. test normal termination and errorjob are ready or molecule is ignored
    extratime = 0
    normaltime= 0
    timestep1 = run.timestep
    timestep2 = max(300, run.timestep)
    print("normal waiting time=NT | extra waiting time=XT")
    while True:
        # CHECK READY:
        print("not ready:", end=' ')
        counter = 0
        for i, mol in enumerate(mols):
            # i. check if mol is ignored or ready
            if mol.ignoremol or mol.IsReady:
                continue

            # ii. check jobs of mol
            for j, job in enumerate(mol.jobs):
                # check if job is ignored or ready
                if job.IsReady or job.ignorejob:
                    continue
                else:
                    counter += 1
                    print("{}.{}".format(i, j), end=' ')

                # test ready. when extratime is too high job is ignored
                job.ready(extratime=extratime, ignore=run.ignore)

                # until now when one job is ignored automatically the whole molecule is discarded
                # this could be more complex if for example only a certain configuration is ignored!
                if job.ignorejob:
                    mol.discard()

            if all(job.IsReady for job in mol.jobs):
                mol.IsReady = True
        print()
        if all(mol.IsReady for mol in mols):
            break
        elif counter == 0:
            print("there are no zzz files anymore?!")
            break

        # test if there are still uncompleted zzz_files in queue. If not count extra time
        if not zzztester(mols):
            extratime += timestep1
        else:
            normaltime += timestep1
        # WAIT:
        time.sleep(timestep1)  # wait 5 minudtes
        # after some time use larger timesteps
        if extratime >= timestep2:
            timestep1 = timestep2
        if extratime:
            print("XT: {:.4f} ||".format(extratime / 3600.), end=' ')
        else:
            print("NT: {:.4f} ||".format(normaltime / 3600.), end=' ')


    # 3. return mols that are not ignored:
    mols_toread = [mol for mol in mols if not mol.ignoremol]
    return mols_toread


def zzztester(mols):
    """
    this function tests if the jobs are still queing or running based on the output of the 'qsta' command
    function needs:
    mols:
    -jobs
    myrun:
    -timestep
    """
    files = []
    for mol in mols:
        jobnames = [job.errorfile for job in mol.jobs if not job.errorfile is None]
        files.extend(jobnames)
    if files:
        #print "zzz-files:", files
        print("n zzz files:", len(files), end=' ')
    else:  # here return so we don't need the qsta
        return False

    # get qstat
    while True:
        qsta_raw = qsta()
        if qsta_raw == False:
            print("qsta not working! trying again after 1 minute")
            time.sleep(60)
        else:
            break

    # get list of jobs and list of their states
    qsta_out = [item.split() for item in qsta_raw.split('\n')]
    states = []
    jobs = []
    for item in qsta_out:
        if len(item) == 4:
            states.append(item[1])
            jobs.append(item[3])
        elif len(item) == 5:
            states.append(item[2])
            jobs.append(item[4])

    # inefficient loop
    #print "files:", files
    #print "jobs:", jobs
    #print "states:", states
    count = 0
    for filetje in files:
        for state, job in zip(states, jobs):
            if job in filetje:  # so there is a zzzjob in the queue!
                # filetje is whole path so job in filetje or filetje.split('/')[-1]==job
                if state in ['Q', 'R', 'H', 'E']:  # so if job still in queue and not has state=='C'
                    count += 1  # so count all the jobs still in queue
                else:
                    assert state == 'C'

    # return answer
    if count > 0:
        return True
    else:
        print("no zzzs (anymore) in queue", end=' ')
        return False


def wait_hasimagfreq(job):
    once = False
    while True:
        datadict = read_file(job)#@@@david: bugfix, can now have calculations without the hasimagfreq keyword in same CINDES run as those with
        #print('ddebug: datadict[hasimagfreq]=',datadict['hasimagfreq'],'type=',type(datadict['hasimagfreq']))
        if datadict.get('hasimagfreq',False):#used to fail when key didn't exist, should now not go into id when no hasimagfreq key@@
            print("still imaginary frequency for job:", job, end=' ')
            #print("!WIP!: will check for AutoMD keyword and when present attempt a manual displacement along imag freq")
            if datadict.get('AutoMD',False):
                import CINDES.scripts.AutoMD_G16
                print("\nAutoMD keyword detected for calculation with imaginary frequency")
                #print('pwd: ',os.getcwd())
                print("Will use following command: MakeG16MDGjf("+str(job).split()[-1]+","+datadict['CalcLine']+","+datadict['Mem']+","+datadict['Cores']+")")
                CINDES.scripts.AutoMD_G16.MakeG16MDGjf(str(job).split()[-1],datadict['CalcLine'],datadict['Mem'],datadict['Cores'] )
                #print("printing datadict next:")
                #print(datadict)
            if not once:
                try:
                    vibfreqs = datadict['vibfreqs']
                    print("vibfreqs:", vibfreqs)
                except KeyError:
                    pass
                once = True
            time.sleep(300)
        else:
            try:
                if datadict['hasimagfreq']==False:
                    print ("no imaginary frequency for job: ", job)
                else:
                    print("ERROR: unexpected value for hasimagfreq keyword job: ", job)
                    print("value hasimagfreq keyword: ", datadict['hasimagfreq'])
            except KeyError:
                print ("no hasimagfreq keyword for job: ", job)
            #print("no imaginary frequency or no hasimagfreq keyword for job:", job, end=' ')
            break
    return


@log_io()
def datareader(mols, run):
    # 1. test normal termination
    mols_toread = normaltermination(mols, run=run)

    # 2. obtain data for each molecule
    for molecule in mols_toread:
        print("><" * 15, molecule, end=' ')
        for job in molecule.jobs:

            # 2.0 optionally: wait until no imag freqs
            if run.hasimagfreq:
                wait_hasimagfreq(job)

            # 2.1 read the job
            readings = read_file(job)

            # 2.2 if there are multiple variants of the job, give each variant a different index _P#
            if hasattr(job, 'pos'):
                for old_key in list(readings.keys()):  # the .keys is very important here. iterkeys or for just readings do not work!
                    readings["{}_P{}".format(old_key, str(job.pos))] = readings.pop(old_key)

            # 2.3 assign Job data to Molecule
            molecule.props.update(readings)

            # 2.4 optionally set geometry data as Molecule attribute
            if job.assign_geom:
                set_geom_attribute(molecule, readings)
        molecule.predicted = False

        # only relevant for stab calculations
        if True:
            setEAHs(molecule)

    return

def set_geom_attribute(molecule, readings):
    from CINDES.utils import utils
    t = utils.PeriodicTable()
    formatstr = lambda x:"{:12.6f}".format(x)
    for key in readings:
        if key.startswith('xyz'):
            coords = readings[key]
            coords = [ list(map(formatstr, xyz)) for xyz in coords ]
            atomnos = readings['atomnos' + key.lstrip('xyz')]
            # add the elements
            for atomn, xyz in zip(atomnos, coords):
                xyz.insert(0, t.element[atomn])
            coordsstr = "\n".join([ " ".join(xyz) for xyz in coords ]) + "\n"
            setattr(molecule, key, coordsstr)
    return

def read_file(Job):
    # 1. look to which calc the logfile belongs when there were simultaneous calculations:
    program = Job.calc['program']
    jobs = Job.calc['jobs']
    filename = Job.logpath


    AutoMD=Job.calc.get('AutoMD',False)#@@@david: WIP added key to automatically MD when imagfreq, currently not used@@
    Cores=str(Job.calc.get('nprocs',4))
    Mem=str(int(Cores)*2)+'GB'#so far not getting specified mem just cores*2
    #CalcLine=Job.calc.jobs["hotline"]#@@

    # 2. split logfile in different jobs
    if program == 'gaussian':
        from CINDES.evaluation.cclib.parser.gaussianparser import Gaussian as Log
        key = 'termination'
        # if keyword freq in line than there is an extra internal job!
        jobslines_v1 = open(filename).read().split(key)[:-1]
        jobslines = []
        for joblines_v1 in jobslines_v1:
            if 'roceeding to internal job step number' in joblines_v1.split('\n', 2)[1]:
                print("freq job appended to main job")
                jobslines[-1] += joblines_v1
            else:
                jobslines.append(joblines_v1)
    elif program == 'orca':
        from CINDES.evaluation.cclib.parser.orcaparser import ORCA as Log
        raise SystemExit('ORCA interface not implemented')
    elif program == 'nwchem':
        from CINDES.evaluation.cclib.parser.nwchemparser import NWChem as Log
        key = 'NWChem Input Module'
        splitted = open(filename).read().split(key)
        jobslines = splitted[1:-1]
    elif program == 'vasp':
        from CINDES.evaluation.cclib.parser.vaspparser import VASP as Log
        jobslines = open(filename).read()
    else:
        raise SystemExit('not implemented')
    print("njobs:", len(jobslines), end=' ')
    if len(jobslines) == 0:
        print("no jobs in logfile!")
        raise SystemExit('should not occur here')

    # 3. handle every subjob as a different logfile and read the needed job['info'] from it
    from io import StringIO
    jobfiles = list(map(StringIO, jobslines))
    datadict = dict()

    datadict['AutoMD']=AutoMD#@@@david
    datadict['Cores']=Cores
    datadict['Mem']=Mem
    #datadict['CalcLine']=CalcLine#@@

    for jobfile, job in zip(jobfiles, jobs):
        datadict['CalcLine']=job["hotline"]
        try:
            job_data = Log(jobfile).parse()
        except Exception as e:
            print("parsing error with:", jobfile, job)
            raise e

        for inf in job['info']:
            if inf == '_':
                continue
            elif inf[0] == 'e':  # so it concerns an energy!:
                datadict[inf] = job_data.scfenergies[-1] / 27.21138505  # this value is used in cclib
            elif inf[0] == 'g' and not inf.startswith('gap'):
                # note that scfenergies are given in eV by cclib but free energy in hartree
                datadict[inf] = job_data.freeenergy
                print("free energy found:", job_data.freeenergy)
            elif inf.startswith('tch'):
                datadict[inf] = job_data.thermalcorrectionH
            elif inf.startswith('tc'):
                datadict[inf] = job_data.thermalcorrectionG
            elif inf.startswith('homo'):
                datadict[inf] = job_data.moenergies[-1][job_data.homos[0]]
            elif inf.startswith('lumo'):
                datadict[inf] = job_data.moenergies[-1][job_data.homos[0] + 1]
            elif inf.startswith('hono_occn'):
                nlumo = ( sum(job_data.atomnos) - job_data.charge ) / 2
                datadict[inf] = job_data.nooccnos[nlumo - 1]
            elif inf.startswith('luno_occn'):
                nlumo = ( sum(job_data.atomnos) - job_data.charge ) / 2
                datadict[inf] = job_data.nooccnos[nlumo]
            elif inf == 'dipole':
                datadict[inf] = job_data.dipole
            elif inf == 'polar':
                datadict['polar'] = sum([job_data.polex[i] / 3 for i in [0, 2, 5]])  # = 1/3*(axx+ayy+azz)
            elif inf == 'mw':
                datadict['mw'] = float(sum(job_data.atomnos))
            elif inf == 'rdv':
                spiden = [x[0] - x[1] for x in zip(job_data.npaa, job_data.npab)]
                datadict['rdv'] = sum([item**2 for item in spiden if abs(item) > 0.05])
            elif inf.startswith('spindensities'):
                # NB charge-beta - charge-alpha because spindensity is a positive value but electron charge is negative!
                spiden = [round(x[1] - x[0], 8) for x in zip(job_data.npaa, job_data.npab)]
                datadict[inf] = spiden
            elif any(prop in inf for prop in ['pcharges', 'partialcharges']):
                #print "partial charges", job_data.atomcharges

                try:
                    pcharges = list(zip(list(map(int, job_data.atomnos)), list(map(round8, job_data.atomcharges['natural']))))
                except KeyError:
                    print("partial charges not found for natural orbitals so Mulliken charges are used")
                    pcharges = list(zip(list(map(int, job_data.atomnos)), list(map(round8, job_data.atomcharges['mulliken']))))
                datadict[inf] = pcharges
            elif inf.startswith('xyz'):
                datadict[inf] = job_data.atomcoords[-1]
                datadict['atomnos' + inf.lstrip('xyz')] = job_data.atomnos
            elif inf.startswith('dipolealpha'):
                datadict[inf] = datadict[inf] = job_data.dipolealpha
            elif inf.startswith('dipolebeta'):
                datadict[inf] = datadict[inf] = job_data.dipolebeta
            elif inf.startswith('electricdipole'):
                datadict[inf] = datadict[inf] = job_data.electricdipole
            elif inf.startswith('vibfreqs'):
                datadict[inf] = job_data.vibfreqs
            elif inf.startswith('hasimagfreq'):
                datadict[inf] = any( freq<0.0 for freq in job_data.vibfreqs )
            elif inf.startswith('nics0'):
                datadict[inf] = job_data.shieldings[0]['Isotropic']
            elif inf.startswith('nicszz0'):
                datadict[inf] = job_data.shieldings[0]['matrix'][2,2]
            elif inf.startswith('nics1'):
                shieldings = job_data.shieldings
                datadict[inf] = 0.5 * ( shieldings[1]['Isotropic'] + shieldings[2]['Isotropic'] )
            elif inf.startswith('nicszz1'):
                shieldings = job_data.shieldings
                datadict[inf] = 0.5 * ( shieldings[1]['matrix'][2,2] + shieldings[2]['matrix'][2,2] )
            elif inf.startswith('chi_0'):
                datadict[inf] = job_data.chi_0
            elif inf.startswith('wbo'):#@@@david: added for the extraction of Wiberg Bond Order
                datadict[inf] = job_data.wbo#@@

            else:
                # try to see if there is an attribute from job_data matching inf. 
                # NB: the else statement of a for loop is executed when the for loop finishes without break statement!
                for attribute in dir(job_data):
                    if inf.startswith(attribute):
                        print("cclib attribute recognized:", attribute)
                        datadict[inf] = getattr(job_data, attribute)
                        break
                else:
                    print("value not recognized:", inf)            
    print()

    return datadict


def set_combined_variables(mol, to_read_props):

    # so now datadict should have all energy keys + homo/lumo + dipole
    # but not yet omega/solv/gap so:
    results = mol.props
    if 'gap' in to_read_props:
        results['gap'] = results['lumo'] - results['homo']
    if 'solv' in to_read_props:
        results['solv'] = (results['e1_solv'] - results['e0_solv']) * 627.5

    # the next properties are possible but not always calculated for radical stability values:
    try:
        if any(i in to_read_props for i in ['ip', 'omega', 'stab']):
            results['ip'] = (results['eIP'] - results['e0']) * 27.2113838
        if any(i in to_read_props for i in ['ea', 'omega', 'stab']):
            results['ea'] = (results['e0'] - results['eEA']) * 27.2113838
        if any(i in to_read_props for i in ['omega', 'stab']):
            results['omega'] = ((results['ip'] + results['ea'])**2) / (8 * (results['ip'] - results['ea']))
    except KeyError:
        if not 'stab' in to_read_props:
            raise

    if 'stab' in to_read_props:
        results = calculate_stab(results, mol)
    if any(i in to_read_props for i in ['ipfukui', 'radfukui']):
        #print "results:", results
        results['ipfukui'] = list(map(round8, [(q_ip[1] - q_0[1]) for q_ip, q_0 in zip(results['pchargesIP'], results['pcharges0'])]))
    if any(i in to_read_props for i in ['eafukui', 'radfukui']):
        results['eafukui'] = list(map(round8, [(q_0[1] - q_ea[1]) for q_0, q_ea in zip(results['pcharges0'], results['pchargesEA'])]))
    if 'radfukui' in to_read_props:
        results['radfukui'] = list(map(round8, [0.5 * (ipf + eaf) for ipf, eaf in zip(results['ipfukui'], results['eafukui'])]))

    if False: # this part can probably be removed
        if 'delta_hardness' in to_read_props:
            results['delta_hardness'] = ( results['lumo_CH3'] - results['homo_CH3'] ) - ( results['lumo_CH2'] - results['homo_CH2'] )
        if 'exaltation' in to_read_props:
            results['exaltation'] = results['chi_0_CH3'] - results['chi_0_CH2']

    return results
