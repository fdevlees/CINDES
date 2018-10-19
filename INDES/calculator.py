"""gaussian functions
The main function of this module is to calculate the values of all the given indices

The input:
    - the indices of the molecules to calculate
    - the configurations of the molecules to calculate (this may be redundant)
    - the data list. The results are appended to this list and returned.
    - the molecular framework given as a dictionary with keys: 'core', 'active' and 'passive'
    - the instance of the Run class. This instance contains all keywords/options that control the workflow

The output:
    - the data that is calculated. It is appended to the given data list!

J.L. Teunissen, 20th June 2016

"""
debug = False
safe = False

import submitter as subm
from copy import deepcopy
import glob
from pprint import pprint
import pprint
import time
#from writings import log_io, sprint
from CINDES.utils.writings import log_io, print_title, sprint
import logging
import construction as zcon
import datareader

# for stab:
import re
import os
import shutil
once = 0

# this function to redirect the output of the optga keyword
import sys
from contextlib import contextmanager


@contextmanager
def custom_redirection(fileobj):
    old = sys.stdout
    sys.stdout = fileobj
    try:
        yield fileobj
    finally:
        sys.stdout = old


def invoke_script(calc, namespace, ID=0):
    if not isinstance(calc, dict):
        return
    if not 'script' in calc:
        print "no script invocation"
        return
    print "in invoke script:):"
    module_obj = __import__(calc['script'])
    module_obj.main(namespace, ID)
    return


def get_secret_data(tablefilename, mols_tocal, mols_nocal, myrun):
    '''checks for confs already calculated:
        uses myrun.~
        -props (set)
    '''
    import json
    with open(tablefilename, 'rb') as f:
        db = json.load(f)
    if debug:
        print "secret_table:"
        sprint(10, secret_table)
    table = {key: value for key, value in db.iteritems() if myrun.props <= value.viewkeys()}
    for mol in mols_tocal[:]:
        try:
            mol.props = table[mol.index]
        except KeyError:
            continue
        mol.predicted = False
        mols_tocal.remove(mol)
        mols_nocal.append(mol)
    return mols_tocal, mols_nocal


def runjobs(mols_tocal, myrun, i):
    calc = myrun.calcs[i]

    def call(function, calc, *args, **kwargs):
        if isinstance(calc, list) or isinstance(calc, tuple):
            for j, cal in enumerate(calc):
                invoke_script(cal, locals(), (i + 1) * 100 + (j + 1))
                function(calc=cal, *args, **kwargs)
        else:
            invoke_script(calc, locals(), i + 1)
            function(calc=calc, *args, **kwargs)
    # 0. filter off ignored molecules
    mols_calc = filter(lambda x: not x.ignoremol, mols_tocal)

    # 1. Make the jobs and add them to the molecules:
    call(jobmaker, mols=mols_calc, myrun=myrun, calc=calc)
    # 2. now the jobs have to be submitted (this function contains a try_ready test)
    jobids = submission(mols_calc, myrun)

    # 3. test of all jobs are ready
    if jobids:
        jobtester(mols_calc, myrun, jobids)
    else:
        print "no jobids so assume no jobs submitted"

    # 4. test normal termination and read jobs (only myrun variable used is actually debug)
    datareader.datareader(mols_calc, myrun)

    return mols_tocal


def do_calcs(mols_tocal, myrun):
    for i, calc in enumerate(myrun.calcs):
        runjobs(mols_tocal, myrun, i)
        # if there need to be set some new geometries for new calculation.
        invoke_script(calc, locals(), (i+1)*100)
        # delete jobs such that new jobs can be set up.
        for mol in mols_tocal:
            mol.deletejobs()

    # after every calculation is performed:
    for mol in mols_tocal:
        if not mol.ignoremol:
            datareader.set_combined_variables(mol, myrun.props)
    return mols_tocal

# PROCEDURE


def procedure(myrun, mols_tocal, mols_nocal):
    global once
    #print "nconfs:", len(population)
    print "| n_indices_tocal:", len(mols_tocal)
    print "|    n_data_nocal:", len(mols_nocal)
    if myrun.no1sub == 1 and once == 0:
        once = 1
        print " "
    elif myrun.secret_file:
        print "SECRET DATA activated:", myrun.secret_file
        tablefilename = myrun.secret_file
        mols_tocal, mols_nocal = get_secret_data(tablefilename, mols_tocal, mols_nocal, myrun)

    if not mols_tocal == []:
        # 0. Set the molecular geometries
        geommaker(mols_tocal, myrun)
        # 1. And perform the calculations
        do_calcs(mols_tocal, myrun)

    # 5. merge data_calc and data_nocal to data_all
    mols_all = mols_tocal + mols_nocal

    # 6. set target property i.e. mol.Pvalue and mol.boundaries
    set_target_properties(mols_all, myrun)

    return mols_all

# 0. geom making


@log_io()
def geommaker(mols_tocal, myrun):

    if myrun.symlinks:
        print "symmetry will be applied |",
    if myrun.optga:
        print "Output of Dihedral GA Optimizer is redirected to optga.out. optimizing... |",
    print
    e = None
    for molecule in mols_tocal:

        # set zmatrices:
        # for multiple zmats:
        if isinstance(myrun.zmatrixfile, list):
            print "myrun.zmatrixfile:", myrun.zmatrixfile
            for zmatfile in myrun.zmatrixfile:
                tzmat = myrun.TZmatrices[zmatfile]
                c, a, p = map(deepcopy, (tzmat['core'], tzmat['active'], tzmat['passive']))
                zmat = zcon.constructor2(molecule.conf, c, a, p, links=myrun.symlinks, defaultgroups=myrun.defaultgroups)
                setattr(molecule, zmatfile, zmat)
        # for a single zmat
        else:
            c = deepcopy(myrun.TZmat['core'])
            a = deepcopy(myrun.TZmat['active'])
            p = deepcopy(myrun.TZmat['passive'])
            zmat = zcon.constructor2(molecule.conf, c, a, p, links=myrun.symlinks, defaultgroups=myrun.defaultgroups)
            setattr(molecule, myrun.zmatrixfile, zmat)
        # always give a default zmat
        setattr(molecule, 'zmat', zmat)
        # for an extra zmat
        if myrun.extrajobs:
            if True:  # i.e. give second geom similar geometry as default geom
                molecule.zmat2 = molecule.zmat

        # Try to print SMILES
        try:
            smiles = molecule.get_format()
            print "smiles:", smiles,
        except (IndexError, NameError, KeyError) as e:
            pass  # only the last error is printed after the whole molecule loop

        if myrun.optga:
            import sys
            from CINDES.utils.ga_dihedrals import reduce_conflicts
            # this function sets molecule.conf with optimized dihedrals in the conf attribute
            c = deepcopy(myrun.TZmat['core'])
            a = deepcopy(myrun.TZmat['active'])
            p = deepcopy(myrun.TZmat['passive'])

            with open('optga.out', 'a') as out:
                with custom_redirection(out):
                    reduce_conflicts(molecule, c, a, p)
        elif False:
            import fafoom
            # conformational analysis could be implemented here
            # the program is not ready for multiple conformations at the moment!
            # possibly the molecule has only 1 index attribute without dihedral angles
            # but multiple conf lists with different dihedral angles.
            # or when we also want to consider more than one dihedral per group
            # for example not only the core-COOH but also the coreCOO-H dihedral
            # then the zma or xyz attribute is a list of multiple zmatrices or cartesian coordinates
            # the the filenames should not be similar so also an index should be included in the filename
            # subsequently the program has to check readyness of all structures
            # and read them
            # and take the props of the lowest energy structure.

            # on the other hand: fafoom could be used to find the configuration with the smallest rdkit-ff configuration
            # nevertheless probably also a tweeked ga_dihedrals with a different eval_function could do this!
        else:
            # no special action. the first assigned geometry is used as a start
            pass
    if e:
        print "OpenBabel Smiles error:", e
    return

# 1. file making


@log_io()
def jobmaker(mols, myrun, calc):  # ----- dict with info for filewriter has to pass here)
    '''jkl'''

    if debug:
        print "DEBUG: calc", calc

    # 1. Decide program
    if calc['program'] == 'gaussian':
        import gaussian as program
    elif calc['program'] == 'nwchem':
        import nwchem as program
    else:
        raise SystemExit('program not recognized')

    # 2. Write inputfile(s)
    for molecule in mols:
        if 'positions' in calc:  # so multiple jobs
            assert 'geom' in calc and (calc['geom'] in ['H', 'AH'])
            # 1. make a folder with the indexname in /data/indices[i]
            if not os.path.exists(calc['path'] + '/' + molecule.index):  # path is $WORKDIR/data
                os.makedirs(calc['path'] + '/' + molecule.index)
                # and make sure ID_gauss is in the folder!
                shutil.copy(calc['path'] + '/' + myrun.script, calc['path'] + '/' + molecule.index)
            # 2. use zmat to make the AH files with the positions stored in fileparameters['positions']
            for pos in calc['positions']:
                zmat2 = deepcopy(molecule.zmat)
                zmat2, h, N = add_hydrogen(zmat2, pos, myrun.ncore)
                attr = 'zmat{}_{:d}'.format(calc['geom'], pos)
                setattr(molecule, attr, zmat2)
                job = program.filewriter(molecule, calc, pos)

                # job.N will be deprecated and replaced by job.Aatom which is the atomic number of 
                # the attached atom. 
                job.N = N
                if job.N:
                    job.Aatom=7
        elif 'fafoom' in calc:  # so first find a lower xyz
            from CINDES.utils import fafoom_utils
            print "trying fafoom..."
            if calc['fafoom'] == 1:
                # only find the lowest conformer
                molecule.xyz = fafoom_utils.GetLowestXYZ(molecule)
                program.filewriter(molecule, calc)
            elif calc['fafoom'] == 2:
                if not os.path.exists(calc['path'] + '/' + molecule.index):  # path is $WORKDIR/data
                    os.makedirs(calc['path'] + '/' + molecule.index)
                    # and make sure ID_gauss is in the folder!
                    shutil.copy(calc['path'] + '/' + myrun.script, calc['path'] + '/' + molecule.index)
                # find a set of conformers
                conformers = fafoom_utils.GetConformers(molecule)
                for i, conformer in enumerate(conformers, 1):  # enumerate starts at 1!
                    setattr(molecule, 'xyz{}'.format(i), conformer.GetProp('xyz'))
                    program.filewriter(molecule, calc, i)
        elif 'geom' in calc and isinstance(getattr(molecule, calc['geom']), list) and calc['geom'][-1]=='s':
            # it is assumed that there are multiple geometries when the geometries attribute ends with an s!
            # 1. make a folder with the indexname in /data/indices[i]
            if not os.path.exists(calc['path'] + '/' + molecule.index):  # path is $WORKDIR/data
                os.makedirs(calc['path'] + '/' + molecule.index)
                # and make sure ID_gauss is in the folder!
                shutil.copy(calc['path'] + '/' + myrun.script, calc['path'] + '/' + molecule.index)
            # 2. 
            for i, geom in enumerate(getattr(molecule, calc['geom'])):
                geomattr = 'xyz{}_{:d}'.format(calc['geom'], i)
                print "geomattr:", geomattr
                setattr(molecule, geomattr, geom)
                job = program.filewriter(molecule, calc, i)
                job.Aatom = molecule.Aatoms[i]
                # this line has to be removed later if job.N is fully deprecated!
                job.N = None
        else:  # so single job
            program.filewriter(molecule, calc)
    return


def add_hydrogen(zmat, pos, ncore):
    """makes new file with hydrogen attached on first dihedral"""
    zmatnew = deepcopy(zmat)
    spos = str(pos)  # spos is string of pos. pos = position
    h = 0
    logging.debug(pprint.pformat(zmat))
    if zmat[pos - 1][0] == 'N':
        N = True
        item = zmat[pos - 1]
        h = 1
        if len(item) == 1:  # when pos is 1 so first index of a zmat
            hline = ['H', 1, 0.9, 2, 109.5, 3, 176.0]
        elif len(item) == 3:  # when pos is 2 so second index of a zmat
            hline = ['H', 2, 0.9, 3, 109.5, 4, 176.0]
        else:
            bondindex = item[1]  # or if item only has length 1
            dihedralindex = item[3]
            hline = ['H', spos, 0.9, bondindex, 109.5, item[3], 176.0]  # LOOK AT THIS
    else:
        N = False
        for item in zmatnew[ncore:]:
            if str(item[1]) == spos:
                if item[0] == 'H':
                    h = 1
                hline = item.copy()
                item[6] = '126.0'
                hline[6] = '234.0'
                hline[2] = 0.9
                hline[0] = 'H'
                #print "hline:",hline
    zmatnew.append(hline)
    return zmatnew, h, N

#    IF I ever want to make the structures with H on the other side attached I need something like this:
# def maker2(zmat,pos,index,**fileparameters):
#    '''makes new file with hydrogen attached on second dihedral.
#    this is only necessary when there is not already another hydrogen on the compound
#    or that the site is nitrogen or possibly sulfur doped. '''
#    zmatnew = zmat[:]
#    spos = str(pos)
#    h=0
#    #print "pos:",spos
#    for item in zmatnew[fileparameters['ncore']:]:
#        if str(item[1]) == spos:
#            hline=item[:]
#            item[6] = '234.0'
#            hline[6] = '126.0'
#            hline[0] = 'H'
#    zmatnew.append(hline)
#    filewriterAH(zmatnew,spos + '_2',index,**fileparameters)
#    return


# 2. submission
@log_io()
def submission(mols_tocal, myrun):
    global once
    if once == 1 and myrun.no1sub == 1:
        print "submit skipped"
        once = 2
        jobids = None
    else:
        jobids = submit_normal(mols_tocal, myrun)  # In here is decided to run on shell or to really submit!
    logging.info("----- END all jobs are submitted ----------")
    if safe:
        time.sleep(15)  # wait 15 seconds. to be sure that the jobs appear in the qstat command
    return jobids


def submit_normal(mols_tocal, myrun):
    ''' this function submits all the jobs of every mol.jobs list
    when worker=True the worker framework will be used:
        first a list of jobs is gathered in jobids and this list is used
        to submit the jobs simultaneously

        worker/1.6.8-intel-2018a
        '''

    jobids = []
    arrayjob = False
    worker = myrun.worker
    for molecule in mols_tocal:
        for job in molecule.jobs:

            # 1. try ready part
            if myrun.try_ready:
                if arrayjob:
                    name = job.logpath
                    if glob.glob(name):
                        print "already calculated:", name
                        continue
                else:
                    name1 = job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                    name2 = job.filepath + '.o[0-9][0-9][0-9][0-9]*'
                    if glob.glob(name1) or glob.glob(name2):
                        print "already calculated:", name2
                        continue

            # 2. submit part
            if worker:
                jobids.append(job)
            else:
                jobid = job.submit()
                jobids.append(jobid)
                time.sleep(1)

    if worker:
        with open('loglist.csv', 'w') as f:
            f.write('job,log\n')
            for job in jobids:
                f.write('{},{}\n'.format(job.filepath,job.logpath))

        jobids = subm.submitworker()
    return jobids

# 3. testing


@log_io(signator='=')
def jobtester(mols_tocal, myrun, jobids=None):
    """ this tester tests if the jobs are ready by looking for a file <name><.extension>.o<6digits>.

        - even if try_ready is activated all indices are used. And the already ready ones are immediately recognized as ready.
        - they are just not submitted again.
        - note that jobids are not used!
        - function returns nothing but returns when all jobs are ready! this function therefore can take very long!
    """
    if jobids is None:
        jobids = []  # this to avoid the mutable default gotcha
    if myrun.nosub == 2 or len(mols_tocal) == 0:
        print "Job tester skipped because jobs are evaluated on login node or no jobs to be calculated"
        return
    print "len(mols_tocal):", len(mols_tocal)
    test_ready = myrun.test_ready
    path = myrun.path
    fileparameters = myrun.__dict__
    if test_ready == 1:
        test_ready1(mols_tocal, myrun)
    elif test_ready == 2:  # default
        test_ready2(mols_tocal, myrun)
    elif test_ready == 3:
        test_ready3(mols_tocal, myrun)
    else:
        raise SystemExit('no valid test_ready value')
    return


def test_ready1(indices, myrun):
    ''' test ready based on the presence of inputfile.o$$ file '''
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = []  # here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + \
            myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab'] == 1:  # property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + \
                    '_' + str(pos) + myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
                paths.append(path2)
    while True:  # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths == []:
            break
        print "time/h:", tijdje / 3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return


def test_ready2(mols_tocal, myrun):
    """
    this function tests if the jobs are still queing or running based on the output of the 'qsta' command
    function needs:
    mols:
    -jobs
    myrun:
    -timelimit
    -timestep
    -extrawaittime
    """
    completedjobs = []
    files = []
    fileparameters = myrun.__dict__

    # 1. get the list of entries that are in the queue
    if myrun.worker:
        files = ['my-gaussian-worker-job']
    else:
        for mol in mols_tocal:
            jobnames = [job.filename for job in mol.jobs]
            files.extend(jobnames)
    print "files:", files

    # 2. and wait until all are completed or not anymore in queue
    tijdje = 0
    while True:
        count = 0
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        filescopy = files[:]
        njobs = len(filescopy)
        qsta_raw = subm.qsta()
        if qsta_raw == False:
            print "qsta not working!"
            time.sleep(fileparameters['timestep'])
            continue
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
        if debug:
            print "states:", states
            print "jobs:", jobs
            print "files:", files
        for filetje in filescopy:
            for state, job in zip(states, jobs):
                if filetje == job:
                    if state in ['Q', 'R']:  # so if job still in queue and not has state=='C'
                        count += 1  # so count all the jobs still in queue
                    elif state in ['H', 'E']:
                        print "ERROR jobs on hold or Error"
                        count += 1
                    else:
                        assert state == 'C'
                        if job not in completedjobs:
                            print 'JOB completed:', job
                            completedjobs.append(job)
                    if debug:
                        print "found a job: ", state, job, filetje
        t = "{:.2f}".format(tijdje / 3600.)
        print "njobs -running: {:d} -ready: {:d} | waittime={} hrs".format(count, njobs - count, t)
        if count == 0:  # so no jobs anymore in queue
            break
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return


def test_ready3(indices, myrun):
    '''this could be something using a line as :
        touch ${PBS_JOBID}.completed
    in the jobscript ID_gauss
    '''
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = []  # here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.completed'
        paths.append(path1)
        if fileparameters['stab'] == 1:  # property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + \
                    indices[i] + '_' + str(pos) + '.completed'
                paths.append(path2)
    while True:  # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje > fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy = paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths == []:
            break
        print "time/h:", tijdje / 3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje += fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime'])  # just wait for the files to write back before opening them
    return
    pass

# 6. set molecular property attributes


def set_target_properties(molecules, myrun):
    ''' set mol.Pvalue and if boundary conditions mol.boundaries
    uses myrun attributes:
        -property
        -function
        -func_args
        -bcprop
    and molecule attributes:
        -props
    and sets molecule attributes:
        -Pvalue
        -boundaries
    '''
    if myrun.property == '_':
        return
    for mol in molecules:
        if mol.ignoremol:
            print mol, 'ignored'
            if myrun.optimum == 'maximum':
                mol.Pvalue = -float("inf")
            else:
                mol.Pvalue = float("inf")
            continue
        if mol.Pvalue:
            print "molecular target property already set. Predicted?", mol
            if myrun.bc:
                print "boundary condition cannot be set"
            print "molecule has probably no .props attribute"
            continue
        if myrun.property == 'func':
            try:
                kwargs = {prop: mol.props[prop] for prop in myrun.func_args}
            except KeyError:
                print "error mol:", mol
                print "props:", mol.props
            mol.Pvalue = myrun.function(**kwargs)
        else:
            try:
                mol.Pvalue = mol.props[myrun.property]
            except KeyError:
                print "error mol:", mol, mol.index
                print "mol.ignoremol", mol.ignoremol
                print "props:", mol.props
                print "molecule doesnt have the required property! so is ignored!"
                mol.discard()

        if myrun.bc:
            try:
                mol.boundaries = [mol.props[bcp] for bcp in [myrun.bcprop]]
            except KeyError as e:
                print e
                pass
    return
