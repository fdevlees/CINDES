#!/bin/env python

from pprint import pprint
#from writings import log_io
from CINDES.utils.writings import log_io, print_title, sprint
from submitter import qsta
import re
import numpy
import time
import logging
#logging.basicConfig(level=logging.DEBUG)

import construction

def round8(x): return round(float(x), 8)


def setEAHs(molecule):
    EAHs = dict()
    Npos = []  # positions with a nitrogen in here
    for job in molecule.jobs:
        if hasattr(job, 'pos'):
            EAHs[job.pos] = {
                'eAH': molecule.props.pop('eAH_P{}'.format(str(job.pos))),
                'Aatom': job.Aatom}
            try:
                EAHs[job.pos]['tchAH'] = molecule.props.pop('tchAH_P{}'.format(str(job.pos)))
            except KeyError:
                pass
    if EAHs:
        molecule.props['EAHs'] = EAHs
        print "EAHs:", EAHs
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
    EAHs = results['EAHs']
    minpos = min(EAHs, key=lambda x: EAHs[x]['eAH'])
    E_ah = EAHs[minpos]['eAH']
    try:
        tch_ah = EAHs[minpos]['tchAH']
    except KeyError:
        pass

    try:
        Domega = results['omega'] - 2.
    except KeyError:
        Domega = 0.0
        print "No electrophilicity term found so stab is calculated without omega term"

    if 'tchA' in results:
        print "applying thermal corrections", results['tchA'], 'and', tch_ah
        results['BDE_ah'] = (results['eA'] + results['tchA'] + H_h - ( E_ah + tch_ah )) * kJmol
    else:
        results['BDE_ah'] = (results['eA'] + E_h - E_ah) * kJmol  # avtc is AVerage Thermal Correction.

    #if EAHs[minpos]['N']:
    if EAHs[minpos]['Aatom']==7:
        chi_term = bde_b * (chi_h - 3) * (chi_n - 3)  # term is independent of the molecule itself. ongeveer 8.4 kJ/mol?
        print "electronegativity correction for nitrogen is used"
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    elif EAHs[minpos]['Aatom']==8:
        chi_term = bde_b * (chi_h - 3) * (chi_o - 3)  # term is independent of the molecule itself. ongeveer 8.4 kJ/mol?
        print "electronegativity correction for OXYGEN is used"
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h - chi_term
    else:
        stab = results['BDE_ah'] - stab_h - bde_a * Domega * Dw_h
    results['H_pos'] = minpos
    results['stab'] = stab
    return results


def calculate_EAHs(results, molecule):
    EAHs = results['EAHs']
    minpos = min(EAHs, key=lambda x: EAHs[x]['eAH'])
    E_ah = EAHs[minpos]['eAH']
    results['eAH'] = E_ah
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
    if notready:
        print "jobs not ready:\n", notready

    # 2. test normal termination and errorjob are ready or molecule is ignored
    extratime = 0
    normaltime= 0
    timestep1 = 10
    timestep2 = 300
    print "normal waiting time=NT | extra waiting time=XT"
    while True:
        # CHECK READY:
        print "not ready:",
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
                    print "{}.{}".format(i, j),

                # test ready. when extratime is too high job is ignored
                job.ready(extratime=extratime, ignore=run.ignore)

                # until now when one job is ignored automatically the whole molecule is discarded
                # this could be more complex if for example only a certain configuration is ignored!
                if job.ignorejob:
                    mol.discard()

            if all(job.IsReady for job in mol.jobs):
                mol.IsReady = True
        print
        if all(mol.IsReady for mol in mols):
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
            print "XT: {:.4f} ||".format(extratime / 3600.),
        else:
            print "NT: {:.4f} ||".format(normaltime / 3600.),


    # 3. return mols that are not ignored:
    mols_toread = filter(lambda mol: not mol.ignoremol, mols)
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
        print "n zzz files:", len(files),
    else:  # here return so we don't need the qsta
        return False

    # get qstat
    while True:
        qsta_raw = qsta()
        if qsta_raw == False:
            print "qsta not working! trying again after 1 minute"
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
        print "no zzzs (anymore) in queue",
        return False


@log_io()
def datareader(mols, run):
    # 1. test normal termination
    mols_toread = normaltermination(mols, run=run)

    # 2. obtain data for each molecule
    for molecule in mols_toread:
        print "><" * 15, molecule,
        for job in molecule.jobs:
            # 2.1 read the job
            readings = read_file(job)

            # 2.2 if there are multiple variants of the job, give each variant a different index _P#
            if hasattr(job, 'pos'):
                for old_key in readings.keys():  # the .keys is very important here. iterkeys or for just readings do not work!
                    readings["{}_P{}".format(old_key, str(job.pos))] = readings.pop(old_key)

            molecule.props.update(readings)
        molecule.predicted = False

        # only relevant for stab calculations
        if True:
            setEAHs(molecule)

    return


def read_file(Job):
    # 1. look to which calc the logfile belongs when there were simultaneous calculations:
    program = Job.calc['program']
    jobs = Job.calc['jobs']
    filename = Job.logpath

    # 2. split logfile in different jobs
    if program == 'gaussian':
        from CINDES.cclib.parser.gaussianparser import Gaussian as Logfile
        key = 'termination'
        # if keyword freq in line than there is an extra internal job!
        jobslines_v1 = open(filename).read().split(key)[:-1]
        jobslines = []
        for joblines_v1 in jobslines_v1:
            if 'roceeding to internal job step number' in joblines_v1.split('\n', 2)[1]:
                print "freq job appended to main job"
                jobslines[-1] += joblines_v1
            else:
                jobslines.append(joblines_v1)
    elif program == 'orca':
        from CINDES.cclib.parser.orcaparser import ORCA as Logfile
        raise SystemExit('ORCA interface not implemented')
    elif program == 'nwchem':
        from CINDES.cclib.parser.nwchemparser import NWChem as Logfile
        key = 'NWChem Input Module'
        splitted = open(filename).read().split(key)
        jobslines = splitted[1:-1]
    else:
        raise SystemExit('not implemented')
    print "njobs:", len(jobslines),
    if len(jobslines) == 0:
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
            if inf == '_':
                continue
            elif inf[0] == 'e':  # so it concerns an energy!:
                datadict[inf] = job_data.scfenergies[-1] / 27.21138505  # this value is used in cclib
            elif inf[0] == 'g' and not inf == 'gap':
                # note that scfenergies are given in eV by cclib but free energy in hartree
                datadict[inf] = job_data.freeenergy
                print "free energy found:", job_data.freeenergy
            elif inf.startswith('tch'):
                datadict[inf] = job_data.thermalcorrectionH
            elif inf.startswith('tc'):
                datadict[inf] = job_data.thermalcorrectionG
            elif inf.startswith('homo'):
                datadict[inf] = job_data.moenergies[-1][job_data.homos[0]]
            elif inf.startswith('lumo'):
                datadict[inf] = job_data.moenergies[-1][job_data.homos[0] + 1]
            elif inf == 'dipole':
                datadict[inf] = job_data.dipole
            elif inf == 'polar':
                datadict['polar'] = sum([job_data.polex[i] / 3 for i in [0, 2, 5]])  # = 1/3*(axx+ayy+azz)
            elif inf == 'mw':
                datadict['mw'] = float(sum(job_data.atomnos))
            elif inf == 'rdv':
                spiden = map(lambda x: x[0] - x[1], zip(job_data.npaa, job_data.npab))
                datadict['rdv'] = sum([item**2 for item in spiden if abs(item) > 0.05])
            elif inf.startswith('spindensities'):
                # NB charge-beta - charge-alpha because spindensity is a positive value but electron charge is negative!
                spiden = map(lambda x: round(x[1] - x[0], 8), zip(job_data.npaa, job_data.npab))
                datadict[inf] = spiden
            elif any(prop in inf for prop in ['pcharges', 'partialcharges']):
                #print "partial charges", job_data.atomcharges

                try:
                    pcharges = zip(map(int, job_data.atomnos), map(round8, job_data.atomcharges['natural']))
                except KeyError:
                    print "partial charges not found for natural orbitals so Mulliken charges are used"
                    pcharges = zip(map(int, job_data.atomnos), map(round8, job_data.atomcharges['mulliken']))
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
            else:
                print "value not recognized:", inf
    print

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
    #elif 'EAHs' in to_read_props:
    #    results = calculate_EAHs(results, mol)
    if any(i in to_read_props for i in ['ipfukui', 'radfukui']):
        #print "results:", results
        results['ipfukui'] = map(round8, [-(q_ip[1] - q_0[1]) for q_ip, q_0 in zip(results['pchargesIP'], results['pcharges0'])])
    if any(i in to_read_props for i in ['eafukui', 'radfukui']):
        results['eafukui'] = map(round8, [-(q_0[1] - q_ea[1]) for q_0, q_ea in zip(results['pcharges0'], results['pchargesEA'])])
    if 'radfukui' in to_read_props:
        results['radfukui'] = map(round8, [.5 * (ipf + eaf) for ipf, eaf in zip(results['ipfukui'], results['eafukui'])])

    if False: # this part can probably be removed
        if 'delta_hardness' in to_read_props:
            results['delta_hardness'] = ( results['lumo_CH3'] - results['homo_CH3'] ) - ( results['lumo_CH2'] - results['homo_CH2'] )
        if 'exaltation' in to_read_props:
            results['exaltation'] = results['chi_0_CH3'] - results['chi_0_CH2']

    return results

