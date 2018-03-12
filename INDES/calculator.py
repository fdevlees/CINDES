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
debug=False
safe=False

import submitter as subm
from copy import  deepcopy
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
once=0

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

def invoke_script(calc, namespace):
    if not isinstance(calc, dict): return
    if not 'script' in calc: 
        print "no script invocation"
        return
    print "in invoke script:):"
    module_obj=__import__(calc['script'])
    module_obj.main(namespace)
    return

def get_secret_data(tablefilename,mols_tocal, mols_nocal, myrun):
    '''checks for confs already calculated:
        uses myrun.~
        -props (set)
    '''
    import json
    with open(tablefilename,'rb') as f:
        db = json.load(f)
    if debug:
        print "secret_table:"
        sprint(10,secret_table)
    table = { key:value for key,value in db.iteritems() if myrun.props <= value.viewkeys() }
    for mol in mols_tocal[:]:
        try:
            mol.props=table[mol.index]
        except KeyError:
            continue
        mol.predicted=False
        mols_tocal.remove(mol)
        mols_nocal.append(mol)
    return mols_tocal, mols_nocal

def runjobs(mols_tocal, myrun, calc):
    def call( function, calc, *args, **kwargs):
        if isinstance(calc, list) or isinstance(calc, tuple):
            for cal in calc: function(calc=cal, *args, **kwargs)
        else: function(calc=calc, *args, **kwargs)
    # 0. filter off ignored molecules
    mols_calc = filter(lambda x:not x.ignore, mols_tocal)

    # 1. Make the jobs and add them to the molecules:
    #call(jobmaker, mols=mols_tocal, myrun=myrun, calc=calc)
    call(jobmaker, mols=mols_calc, myrun=myrun, calc=calc)
    # 2. now the jobs have to be submitted (this function contains a try_ready test)
    jobids = submission(mols_calc, myrun)
 
    # 3. test of all jobs are ready
    if jobids: jobtester(mols_calc,myrun,jobids)
    else: print "no jobids so assume no jobs submitted"
 
    # 4. test normal termination and read jobs (only myrun variable used is actually debug)
    datareader.datareader(mols_calc, myrun)

    return mols_tocal

def do_calcs(mols_tocal, myrun):
    for i, calc in enumerate(myrun.calcs):
        runjobs(mols_tocal, myrun, calc=myrun.calcs[i])

        # if there need to be set some new geometries for new calculation.
        invoke_script(calc, locals())
        
        # 5. delete jobs such that new jobs can be set up.
        for mol in mols_tocal: mol.deletejobs()

    # after every calculation is performed:
    for mol in mols_tocal: 
        if not mol.ignore: datareader.set_combined_variables(mol, myrun.props)
        
    return mols_tocal

# PROCEDURE
def procedure(myrun, mols_tocal, mols_nocal, TZMat):
    global once
    #print "nconfs:", len(population)
    print "| n_indices_tocal:", len(mols_tocal)
    print "|    n_data_nocal:", len(mols_nocal)
    if myrun.no1sub==1 and once==0:
        once = 1
        print " "
    elif myrun.secret_file:
        print "SECRET DATA activated:", myrun.secret_file
        tablefilename = myrun.secret_file
        mols_tocal , mols_nocal = get_secret_data(tablefilename, mols_tocal, mols_nocal, myrun)

    if not mols_tocal==[]:
        # 0. Set the molecular geometries
        geommaker(mols_tocal,myrun,**TZMat)
        # 1. And perform the calculations
        do_calcs(mols_tocal, myrun)

    # 5. merge data_calc and data_nocal to data_all
    mols_all = mols_tocal + mols_nocal

    # 6. set target property i.e. mol.Pvalue and mol.boundaries
    set_target_properties( mols_all, myrun)

    return mols_all

#0. geom making
@log_io()
def geommaker(mols_tocal,myrun,passive, active, core):
    if myrun.symlinks: print "symmetry will be applied |",
    if myrun.optga: print "Output of Dihedral GA Optimizer is redirected to optga.out. optimizing... |",
    print
    e=None
    for molecule in mols_tocal:
        c = deepcopy(core)
        a = deepcopy(active)
        p = deepcopy(passive)
        molecule.set_zmat( zcon.constructor2(molecule.conf,c,a,p, links=myrun.symlinks) )
        if myrun.extrajobs:
            if True: # i.e. give second geom similar geometry as default geom
                molecule.zmat2=molecule.zmat

        # Try to print SMILES
        try:
            smiles= molecule.get_format()
            print "smiles:", smiles,
        except (IndexError, NameError, KeyError) as e:
            pass # only the last error is printed after the whole molecule loop

        if myrun.optga:
            import sys
            from CINDES.utils.ga_dihedrals import reduce_conflicts
            # this function sets molecule.conf with optimized dihedrals in the conf attribute
            c = deepcopy(core)
            a = deepcopy(active)
            p = deepcopy(passive)

            with open('optga.out','a') as out:
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
    if e: print "OpenBabel Smiles error:", e
    return

# 1. file making
@log_io()
def jobmaker(mols,myrun, calc): #----- dict with info for filewriter has to pass here)
    '''jkl'''

    # 1. Decide program
    if calc['program']=='gaussian':
        import gaussian as program
    elif calc['program']=='nwchem':
        import nwchem as program
    else:
        raise SystemExit('program not recognized')

    # 2. Write inputfile(s)
    for molecule in mols:
        if 'positions' in calc: # so multiple jobs
            if 'geom' in calc and (calc['geom'] in ['H', 'AH']):
                # 1. make a folder with the indexname in /data/indices[i]
                if not os.path.exists(calc['path'] + '/' + molecule.index): #path is $WORKDIR/data
                    os.makedirs(calc['path'] + '/' + molecule.index)
                    # and make sure ID_gauss is in the folder!
                    shutil.copy(calc['path'] + '/' + myrun.script, calc['path']+'/'+molecule.index)
                # 2. use zmat to make the AH files with the positions stored in fileparameters['positions']
                for pos in calc['positions']:
                    zmat2 = deepcopy(molecule.zmat)
                    zmat2, h, N = add_hydrogen(zmat2, pos, myrun.ncore)
                    setattr(molecule, 'zmat{}'.format(pos), zmat2)
                    job = program.filewriter(molecule, calc, pos)
                    job.N=N
                    # here set somehow if the pos belongs to nitrogen
                    #raise NotImplementedError('here implement Npos')
        else: # so single job
            program.filewriter(molecule, calc) #------------------------------------------------HERE IS THE FILEWRITER CALL
    return

def add_hydrogen(zmat, pos, ncore):
    """makes new file with hydrogen attached on first dihedral"""
    zmatnew = deepcopy(zmat)
    spos = str(pos) #spos is string of pos. pos = position
    h=0
    logging.debug(pprint.pformat(zmat))
    if zmat[pos-1][0] == 'N':
        N=True
        item = zmat[pos-1]
        h=1
        if len(item)==1: #when pos is 1 so first index of a zmat
            hline = ['H',1,0.9,2,109.5,3,176.0]
        elif len(item)==3: #when pos is 2 so second index of a zmat
            hline = ['H',2,0.9,3,109.5,4,176.0]
        else:
            bondindex=item[1] #or if item only has length 1
            dihedralindex=item[3]
            hline= ['H',spos,0.9,bondindex,109.5,item[3],176.0]#LOOK AT THIS
    else:
        N=False
        for item in zmatnew[ncore:]:
            if str(item[1]) == spos:
                if item[0] == 'H': h=1
                hline=item.copy()
                item[6] = '126.0'
                hline[6] = '234.0'
                hline[2] = 0.9
                hline[0] = 'H'
                #print "hline:",hline
    zmatnew.append(hline)
    return zmatnew, h, N

#    IF I ever want to make the structures with H on the other side attached I need something like this:
#def maker2(zmat,pos,index,**fileparameters):
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
def submission(mol_tocal,myrun):
    global once
    if once==1 and myrun.no1sub==1:
        print "submit skipped"
        once = 2
        jobids = None
    else:
        #MOST IMPORTANT PART
        #if myrun.try_ready==1:
        #    print "try_ready activated"
        #    mol_tosubmit = try_ready_test(mol_tocal, myrun.__dict__)
        #else: mol_tosubmit=mol_tocal
        #jobids = submit_normal(mol_tosubmit, myrun) #In here is decided to run on shell or to really submit!
        jobids = submit_normal(mol_tocal, myrun) #In here is decided to run on shell or to really submit!
    logging.info("----- END all jobs are submitted ----------")
    if safe: time.sleep(15) # wait 15 seconds. to be sure that the jobs appear in the qstat command
    return jobids

'''
def try_ready_test(mol_tocal, fileparameters):
    """ Jobtester 3 looks which files shouldn't be submitted anymore. These are removed from the indices list and this list is returned

        - It tested if the .com.o123899 file already exists. Actually it should test if the logfile ends in normal termination.?
        - Note that this function does return new indices and no jobids
    """
    #extension=fileparameters['extension']
    #mol_submit = [ mol.copy() for mol in mol_tocal ] 
    arrayjob=False
    for mol in mol_tocal:
        for job in mol.jobs:
            if arrayjob:name=job.logpath
            else:name=job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
            print "name=:", name
            if glob.glob(name): # test A
                print "already calculated:", name
                mol.jobs.remove(job)
        #if not mol.jobs: #so if it is an empty list
        #    mol_submit.remove(mol)
    return mol_tocal
'''

def submit_normal(mols_tocal, myrun):
    jobids = []
    arrayjob=False
    for molecule in mols_tocal:
        for job in molecule.jobs:

            # 1. try ready part
            if myrun.try_ready:
                if arrayjob:name=job.logpath
                else:name=job.filepath[:-4] + '.o[0-9][0-9][0-9][0-9]*'
                #print "name=:", name
                if glob.glob(name): # test A
                    print "already calculated:", name
                    continue

            # 2. submit part
            if job.calc['nosub'] ==2:
                time.sleep(1)
                subm.nosubmit(job)
                print '{} submitted'.format(job)
            else:
                jobid = subm.submit(job, myrun.script).strip()
                jobids.append(jobid)
    return jobids

# 3. testing
@log_io(signator='=')
def jobtester(mols_tocal,myrun,jobids=[]):
    """ this tester tests if the jobs are ready by looking for a file <name><.extension>.o<6digits>.

        - even if try_ready is activated all indices are used. And the already ready ones are immediately recognized as ready. 
        - they are just not submitted again.
        - note that jobids are not used!
        - function returns nothing but returns when all jobs are ready! this function therefore can take very long!
    """
    if myrun.nosub==2 or len(mols_tocal)==0:
        print "Job tester skipped because jobs are evaluated on login node or no jobs to be calculated"
        return
    print "len(mols_tocal):", len(mols_tocal)
    test_ready = myrun.test_ready
    path = myrun.path
    fileparameters = myrun.__dict__
    if test_ready==1:
        test_ready1(mols_tocal,myrun)
    elif test_ready==2: # default
        test_ready2(mols_tocal,myrun)
    elif test_ready==3:
        test_ready3(mols_tocal,myrun)
    else:
        raise SystemExit('no valid test_ready value')
    return

def test_ready1(indices,myrun):
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = [] #here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + myrun.extension + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
                paths.append(path2)
    while True: # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje>fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy= paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths==[]:
            break
        print "time/h:", tijdje/3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje+=fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime']) #just wait for the files to write back before opening them
    return

def test_ready2(mols_tocal,myrun):
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
    files=[]
    fileparameters=myrun.__dict__
    for mol in mols_tocal:
        jobnames = [ job.filename for job in mol.jobs ]
        files.extend(jobnames)
    print "files:", files
    tijdje = 0
    while True:
        count=0
        if tijdje>fileparameters['timelimit']:
            print "time is up"
            break
        filescopy = files[:]
        njobs = len(filescopy)
        qsta_raw = subm.qsta()
        if qsta_raw==False:
            print "qsta not working!"
            time.sleep(fileparameters['timestep'])
            continue
        qsta_out = [ item.split() for item in subm.qsta().split('\n') ]
        #states,jobs = zip(*[ (item[2],item[4]) for item in qsta_out if len(item)>4 ])
        states = []
        jobs = []
        for item in qsta_out:
            if len(item)==4:
                states.append(item[1])
                jobs.append(item[3])
            elif len(item)==5:
                states.append(item[2])
                jobs.append(item[4])
        if debug:
            print "states:",states
            print "jobs:",jobs
            print "files:", files
        for filetje in filescopy:
            for state,job in zip(states,jobs):
                if filetje==job:
                    if state in ['Q','R']: #so if job still in queue and not has state=='C'
                        count += 1 #so count all the jobs still in queue
                    elif state in ['H','E']:
                        print "ERROR jobs on hold or Error"
                        count +=1
                    else:
                        assert state=='C'
                        if job not in completedjobs:
                            print 'JOB completed:',job
                            completedjobs.append(job)
                    if debug:
                        print "found a job: ", state, job, filetje
        print "there are still %d jobs in queue and %d jobs are ready | waittime=%f uur" % (count, njobs-count,float(tijdje)/3600.)
        if count==0: #so no jobs anymore in queue
            break
        time.sleep(fileparameters['timestep'])
        tijdje+=fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime']) #just wait for the files to write back before opening them
    return

def test_ready3(indices,myrun):
    '''this could be something using a line as :
        touch ${PBS_JOBID}.completed
    in the jobscript ID_gauss
    '''
    path = myrun.path
    fileparameters = myrun.__dict__
    tijdje = 0
    paths = [] #here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.completed'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.completed'
                paths.append(path2)
    while True: # then we remove each item of the paths that exists. If every path exists, all jobs are ready
        if tijdje>fileparameters['timelimit']:
            print "time is up"
            break
        pathscopy= paths[:]
        for pathje in pathscopy:
            if glob.glob(pathje):
                paths.remove(pathje)
                print "ready: ", pathje[:-25]
        if paths==[]:
            break
        print "time/h:", tijdje/3600, "len paths:", len(paths),
        time.sleep(fileparameters['timestep'])
        tijdje+=fileparameters['timestep']
    logging.info("All jobs are READY")
    time.sleep(fileparameters['extrawaittime']) #just wait for the files to write back before opening them
    return
    pass

#6. set molecular property attributes
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
    #print "I'm here: molecules:", molecules,
    for mol in molecules:
        if mol.ignore:
            print mol, 'ignored'
            if myrun.optimum=='maximum':mol.Pvalue=float("inf")
            else: mol.Pvalue=-float("inf")
            continue
        if mol.Pvalue:
            print "molecular target property already set. Predicted?", mol
            if myrun.bc: print "boundary condition cannot be set"
            print "molecule has probably no .props attribute"
            continue
        if myrun.property=='func':
            kwargs = { prop:mol.props[prop] for prop in myrun.func_args }
            mol.Pvalue = myrun.function(**kwargs)
            #print "function value:", mol.Pvalue
        else:
            #print "I'm here too:", mol.props
            #print "myrun.property:", myrun.property
            mol.Pvalue = mol.props[ myrun.property ]
            #print "myrun.Pvalue:", mol.Pvalue

        if myrun.bc:
            try:
                #print "I'm here: myrun.bcprop", myrun.bcprop, "mol.props?:", mol.props
                mol.boundaries = [ mol.props[bcp] for bcp in [myrun.bcprop] ]
            except KeyError as e:
                print e
                pass
    #print "test json attributes:"
    #for mol in molecules:
    #    print mol.index, mol.boundaries, mol.Pvalue
    return


