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
debug=0
safe=False

import submitter as subm
from copy import  deepcopy
import glob
from pprint import pprint
import pprint
import time
#from writings import log_io, sprint
from CINDES4.utils.writings import log_io, print_title, sprint
import logging
import construction as zcon
import datareader

# for stab:
import re
import os
import shutil
once=0

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

# PROCEDURE
def procedure(myrun, mols_tocal, mols_nocal, TZmat):
    global once
    #print "nconfs:", len(population)
    print "| n_indices_tocal:", len(mols_tocal)
    print "|    n_data_nocal:", len(mols_nocal)
    if myrun.no1sub==1 and once==0:
        once = 1
        print " "
    elif myrun.nosub==3:
        print "SECRET DATA activated:", myrun.nosub_file
        tablefilename = myrun.nosub_file
        mols_tocal , mols_nocal = get_secret_data(tablefilename, mols_tocal, mols_nocal, myrun)

    if not mols_tocal==[]:
        # 1. Make the files
        filemaker(mols_tocal,myrun,**TZmat) #----------------------------------HERE IS THE FILEWRITER CALL

        # 2. now the jobs have to be submitted 
        jobids = submission(mols_tocal,myrun)

        # 3. test of all jobs are ready
        jobtester(mols_tocal,myrun,jobids)

        # 4. test normal termination and read jobs 
        mols_calc = datareader.datareader(mols_tocal,myrun.__dict__)

    else: mols_calc = []

    # 5. merge data_calc and data_nocal to data_all
    mols_all = mols_calc + mols_nocal

    # 6. set target property i.e. mol.Pvalue and mol.boundaries
    set_target_properties( mols_all, myrun)

    return mols_all

# 1. file making
@log_io()
def filemaker(mols_tocal,myrun,passive,active,core): #----- dict with info for filewriter has to pass here)
    ''' jkl'''
    path = myrun.path
    fileparameters = myrun.__dict__
    for molecule in mols_tocal:
        c = deepcopy(core)
        a = deepcopy(active)
        p = deepcopy(passive)
        molecule.set_zmat( zcon.constructor2(molecule.conf,c,a,p, links=myrun.symlinks) )
        if not myrun.stab==1:
            zcon.filewriter2(molecule.zmat,molecule.index,**fileparameters) #------------------------------------------------HERE IS THE FILEWRITER CALL
        else:
            # 1. WRITE radical input with filewriterA
            zcon.filewriterA(molecule.zmat,molecule.index,**fileparameters) #here we have to use makers to construct the AH files

            # 2. make a folder with the indexname in /data/indices[i]
            if not os.path.exists(path + '/' + molecule.index): #path is $WORKDIR/data
                os.makedirs(path + '/' + molecule.index)
                # and make sure ID_gauss is in the folder!
                shutil.copy(path +'/ID_gauss',path+'/'+molecule.index)

            # 3. reopen written A-file to extract Z-matrix to make the AH files
            filename = fileparameters['path'] + '/' + fileparameters['identify'] + molecule.index + ".com" #same line as in filewriter. open it again.
            zmat = extract_zmat(filename)

            # 4. use zmat to make the AH files with the positions stored in fileparameters['positions']
            for pos in fileparameters['positions']:
                zmat2 = deepcopy(zmat)
                hornot = maker1(zmat2,pos,molecule.index,**fileparameters) #returns a value indicating if there is already a hydrogen (or a nitrogen)
                # FOR NOW ONLY DO ONE POSSIBILITY THIS IS EASIER BECAUSE WE KNOW EXACTLY HOW MANY JOBS THERE HAVE TO BE SUBMITTED
                #if not hornot == 1: #if not there are two ways to place the hydrogen.
                    #maker2(zmat,pos,indices[i],**fileparameters)

        # Try to print SMILES
        try:
            smiles= molecule.get_format()
            print "smiles:", smiles,
        except IndexError:
            print "IndexError while trying to make smiles for molecule"
        except NameError:
            print "NameError while trying to make smiles for molecule"
        except KeyError:
            print "KeyError while trying to make smiles for molecule"
    return

def extract_zmat(filename):
    """This function is used to make the A-H files for the stab property"""

    logging.debug("filename: " + filename)
    multcharge = re.compile('^\-?[0-9]\s[0-9]') #only set the compiler
    with open(filename) as fid: #again open as fid
        for line in fid:
            if multcharge.match(line): #from where there is a match it reads the subsequant lines as the zmat
                zmat=[]
                line=next(fid)
                while not line == '\n': #until empty line
                    zmat.append(line.split())
                    line=next(fid)
                break #so that only the first match is used. after the other matches there is no zmat
    return zmat

def maker1(zmat,pos,index,**fileparameters):
    """makes new file with hydrogen attached on first dihedral"""
    zmatnew = deepcopy(zmat)
    spos = str(pos) #spos is string of pos. pos = position
    h=0
    #print "zmat[spos-1]:",zmat[pos-1]
    #print "pos:",spos
    #print "fileparamters ncore:", fileparameters['ncore']
    logging.debug(pprint.pformat(zmat))
    #print "zmat[fileparameters['ncore']:]:"
    #pp.pprint(zmat[fileparameters['ncore']:])
    if zmat[pos-1][0] == 'N':
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
        for item in zmatnew[fileparameters['ncore']:]:
            #print "in loop", "spos:",spos,"str(item[1]",str(item[1])
            if str(item[1]) == spos:
                if item[0] == 'H': h=1
                hline=item[:]
                item[6] = '126.0'
                hline[6] = '234.0'
                hline[0] = 'H'
                #print "hline:",hline
    zmatnew.append(hline)
    zcon.filewriterAH(zmatnew,spos,index,**fileparameters) #now it is important where this will be written.
    return h

def maker2(zmat,pos,index,**fileparameters):
    '''makes new file with hydrogen attached on second dihedral.
    this is only necessary when there is not already another hydrogen on the compound
    or that the site is nitrogen or possibly sulfur doped. '''
    zmatnew = zmat[:]
    spos = str(pos)
    h=0
    #print "pos:",spos
    for item in zmatnew[fileparameters['ncore']:]:
        if str(item[1]) == spos:
            hline=item[:]
            item[6] = '234.0'
            hline[6] = '126.0'
            hline[0] = 'H'
    zmatnew.append(hline)
    zcon.filewriterAH(zmatnew,spos + '_2',index,**fileparameters)
    return

# 2. submission
@log_io()
def submission(mol_tocal,myrun):
    global once
    fileparameters = myrun.__dict__
    #mol_tocal_all = deepcopy(mol_tocal) #here i copy the indices. The indices are submitted. The indicesall are not all submitted but are all read out.
    if once==1 and fileparameters['no1sub']==1:
        print "submit skipped"
        once = 2
        jobids = None
    elif myrun.stab==1: #then submit also the jobs in folders
        if fileparameters['try_ready']==1:
            print "try_ready activated"
            mol_submit = try_ready_test(mol_tocal,myrun.path,fileparameters,returnpath=False)
            print "mol_submit:", mol_submit
            jobids = submit_stab(mol_submit,myrun)
        else:
            jobids = submit_stab(mol_tocal,myrun)
    else:
        #MOST IMPORTANT PART
        if myrun.try_ready==1:
            print "try_ready activated"
            mol_submit = try_ready_test(mol_tocal,myrun.path,fileparameters)
            jobids = submit_normal(mol_submit, myrun)
        else:
            jobids = submit_normal(mol_tocal,myrun) #In here is decided to run on shell or to really submit!
    logging.info("----- END all jobs are submitted ----------")
    if safe: time.sleep(15) # wait 15 seconds. to be sure that the jobs appear in the qstat command
    return jobids

def try_ready_test(mol_tocal,path,fileparameters,returnpath=False):
    """ Jobtester 3 looks which files shouldn't be submitted anymore. These are removed from the indices list and this list is returned

        - It tested if the .com.o123899 file already exists. Actually it should test if the logfile ends in normal termination.?
        - Note that this function does return new indices and no jobids
    """

    arrayjob=True

    # 1. make a list of paths that need to exist when job is ready
    paths = [] #here we are going to make a list of paths of the jobs
    if 'positions' in fileparameters: positions = fileparameters['positions']
    #for i in range(len(indices)):
    for mol in mol_tocal:
        if arrayjob:
            path1 = path + '/' + fileparameters['identify'][:-1] + '*_' + mol.index + '.log'
        else:
            path1 = path + '/' + fileparameters['identify'][:-1] + '*_' + mol.index + '.com.o[0-9][0-9][0-9][0-9][0-9]*'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + mol.index + '/' + fileparameters['identify'] + mol.index + '_' + str(pos) + '.com.o[0-9][0-9][0-9][0-9][0-9]*'
                paths.append(path2)

    newpaths = paths[:]
    mol_submit = mol_tocal[:]
    if fileparameters['stab']==1: #test if all necessary A and AH calculations are performed
        k=0
        #for i in range(len(indices)): # all indices
        for mol in mol_tocal:
            l=0
            if glob.glob(paths[k]): # test A
                print "already calculated:", paths[k]
                newpaths.remove(paths[k])
                l+=1
            for j in range(len(positions)): # test all AH
                k +=1
                if glob.glob(paths[k]):
                    print "already calculated:", paths[k]
                    newpaths.remove(paths[k])
                    l+=1
            print "len(positions):", len(positions)
            print "l:", l
            if l == len(positions) + 1: #if all AH and A then remove from indices
                mol_submit.remove(mol)
            k+=1
    else:
        #for i in range(len(paths)):
        for mol, path in zip(mol_tocal, paths):
            if glob.glob(path):
                print "already calculated:", mol
                mol_submit.remove(mol)
    if returnpath:
        return mol_submit,newpaths
    else:
        return mol_submit

def submit_normal(mols_tocal,myrun):
    jobids = []
    for molecule in mols_tocal:
        name = molecule.index + '.com'
        if myrun.nosub ==2:
            time.sleep(1)
            jobid = subm.nosubmit(myrun.path,molecule.index ,myrun.identify)
            print molecule.index + 'submitted'
        else:
            #print "name:", name
            #print "myrun.path:", myrun.path
            #print "myrun.identify:", myrun.identify
            jobid = subm.submit(myrun.path,name,myrun.identify).strip()
        jobids.append(jobid)
    return jobids

def submit_stab(mol_submit,myrun,jobids=[]):
    path = myrun.path
    for molecule in mol_submit:
        name1 = molecule.index + '.com'
        jobid = subm.submit(path,name1,myrun.identify).strip()
        jobids.append(jobid)
        for pos in myrun.positions:
            path2 = path + '/' + molecule.index
            name2 = molecule.index + '_' + str(pos) + '.com'
            jobid = subm.submit(path2,name2,myrun.identify).strip()
            jobids.append(jobid)
    return jobids

# 3. testing
@log_io(signator='=')
def jobtester(mols_tocal,myrun,jobids=[]):
    """ this tester tests if the jobs are ready by looking for a file <name>.com.o<6digits>.

        - even if try_ready is activated all indices are used. And the already ready ones are immediately recognized as ready. 
        - they are just not submitted again.
        - note that jobids are not used!
        - function returns nothing but returns when all jobs are ready! this function therefore can take very long!
    """
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
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.com.o[0-9][0-9][0-9][0-9][0-9][0-9]'
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
    indices = [ mol.index for mol in mols_tocal ]
    completedjobs = []
    fileparameters = myrun.__dict__
    tijdje = 0
    files = [] #here we are going to make a list of filenames of the jobs
    for i in range(len(indices)):
        file1 = fileparameters['identify'] + indices[i] + '.com'
        files.append(file1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                file2 = fileparameters['identify'] + indices[i] + '_' + str(pos) + '.com'
                files.append(file2)
    while True:
        count=0
        if tijdje>fileparameters['timelimit']:
            print "time is up"
            break
        filescopy = files[:]
        njobs = len(filescopy)
        qsta_raw = subm.qsta()
        if qsta_raw==False:
            print "No jobs!"
            break
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
    for mol in molecules:
        if mol.Pvalue:
            print "molecular target property already set. Predicted?", mol
            if myrun.bc: print "boundary condition cannot be set"
            print "molecule has probably no .props attribute"
            continue
        if myrun.property=='func':
            kwargs = { prop:mol.props[prop] for prop in myrun.func_args }
            mol.Pvalue = myrun.function(**kwargs)
            print "function value:", mol.Pvalue
        else:
            mol.Pvalue = mol.props[ myrun.property ]
        if myrun.bc:
            try:
                print "I'm here: myrun.bcprop", myrun.bcprop, "mol.props?:", mol.props
                mol.boundaries = [ mol.props[bcp] for bcp in [myrun.bcprop] ]
            except KeyError as e:
                print e
                pass
    #print "test json attributes:"
    #for mol in molecules:
    #    print mol.index, mol.boundaries, mol.Pvalue
    return


