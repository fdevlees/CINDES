import logging
import glob
import time
from copy import deepcopy
import shutil
import os
import re
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
#logging.basicConfig(format='%(levelname)s:%(message)s',level=logging.DEBUG)
import pprint
import construction as zcon
import submitter as subm
import orcareader
# here orca functions

def submittingprocedure(core,active,passive,confs,indices,data,fileparameters):
    # here submitting thing knows at least the path
    path = fileparameters['path']
    filemaker(core,active,passive,confs,indices,**fileparameters) #----------------------------------HERE IS THE FILEWRITER CALL
    print "----- END making of the files -------------"
    #---           -------------           ---------------           --------------           ---#
    #   SYSTEM EXIT             SYSTEM EXIT               SYSTEM EXIT              SYSTEM EXIT
    #---           -------------           ---------------           --------------           ---#
#    raise SytemExit('Exit')
    # now the jobs have to be submitted 
    jobids = []
    if fileparameters['stab']==1:#then submit also the jobs in folders
        for item in indices:
            jobid = subm.submit(path,item,fileparameters['identify'],script='ID_orca').strip()
            jobids.append(jobid)
            for pos in fileparameters['positions']:
                path2 = path + '/' + item
                item2 = item + '_' + str(pos)
                jobid = subm.submit(path2,item2,fileparameters['identify'],script='ID_orca').strip()
                jobids.append(jobid)
    else:
        for item in indices:
            if fileparameters['nosub']==2:
                jobid = subm.nosubmit_orca(path,item,fileparameters['identify'])
                print item + 'submitted'
            else:
                jobid = subm.submit(path,item,fileparameters['identify'],script='ID_orca').strip()
            jobids.append(jobid)
    print "----- all jobs are submitted ----------"
    #status = subm.jobstatus(jobid) #CANNOT WORK WITH STAB
    # test of all jobs are ready | later change to two minutes or so. 
    jobtester(indices,jobids,path,fileparameters)
    data = orcareader.datareader(indices,jobids,path,data,fileparameters)
#    raise SystemExit(0)
    return data

def jobtester(indices,jobids,path,fileparameters):
    tijdje = 0
    paths = [] #here we are going to make a list of paths of the jobs
    for i in range(len(indices)):
        path1 = path + '/' + fileparameters['identify'] + indices[i] + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
        paths.append(path1)
        if fileparameters['stab']==1: #property is global variable
            for pos in fileparameters['positions']:
                path2 = path + '/' + indices[i] + '/' + fileparameters['identify'] + indices[i] + '_' + str(pos) + '.o[0-9][0-9][0-9][0-9][0-9][0-9]'
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
    print "All jobs are READY"
    if fileparameters['nosub']==0:
        time.sleep(180) #just wait for the files to write back before opening them
    return           

def filemaker(core,active,passive,confs,indices,**fileparameters): #----- dict with info for filewriter has to pass here
    from copy import deepcopy
    for i in range(len(confs)):
	c = deepcopy(core)
	a = deepcopy(active)
	p = deepcopy(passive)
	logging.debug("i=" + str(i))
	mat = zcon.constructor2(confs[i],c,a,p)
        logging.debug("in filemaker kwargs:")
        logging.debug(pprint.pformat(fileparameters))
        #now use this variable to test if we have to make all kind of extra files or not

        if not fileparameters['stab']==1:
	    orca_writer(mat,indices[i],**fileparameters) #------------------------------------------------HERE IS THE FILEWRITER CALL
        else:
            orca_writerA(mat,indices[i],**fileparameters) #here we have to use makers to construct the AH files
            #we make a folder with the indexname in /data/indeces[i]
            if not os.path.exists(fileparameters['path'] + '/' + indices[i]): #path is $WORKDIR/data
                os.makedirs(fileparameters['path'] + '/' + indices[i])
                # and make sure ID_gauss is in the folder!
                shutil.copy(fileparameters['path'] +'/ID_orca',fileparameters['path']+'/'+indices[i])
            filename = fileparameters['path'] + '/' + fileparameters['identify'] + str(indices[i]) #same line as in filewriter. open it again.
# ----part of changer
            logging.debug("filename: " + filename)
            multcharge = re.compile('^\*\sint\s\-?[0-9]\s[0-9]') #only set the compiler
            with open(filename) as fid: #again open as fid
                for line in fid:
                    if multcharge.match(line): #from where there is a match it reads the subsequant lines as the zmat
                        print "match!"
                        print multcharge.match(line).group()
                        zmat=[]
                        line=next(fid)
                        while not line == '*\n': #until empty line
                            zmat.append(line.split())
                            line=next(fid)
                        break #so that only the first match is used. after the other matches there is no zmat
            for pos in fileparameters['positions']:
                zmat2 = deepcopy(zmat)
                hornot = orca_maker1(zmat2,pos,indices[i],**fileparameters) #returns a value indicating if there is already a hydrogen (or a nitrogen)
                # FOR NOW ONLY DO ONE POSSIBILITY THIS IS EASIER BECAUSE WE KNOW EXACTLY HOW MANY JOBS THERE HAVE TO BE SUBMITTED
                #if not hornot == 1: #if not there are two ways to place the hydrogen.
                    #orca_maker2(zmat,pos,indices[i],**fileparameters)
            print "DONE"
    return
# ---- end part of changer

def orca_maker1(zmat,pos,index,**fileparameters):
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
    if zmat[pos-1][0] == 'N': #there is nothing attached yet so cannot copy 
        item = zmat[pos-1] #line of N
        h=1
        if pos==1: #when pos is 1 so first index of a zmat
            #hline = ['H',1,0.9,2,109.5,3,126.0]
            hline = ['H',1,2,5,0.9,109.5,126.0]
        elif pos==2: #when pos is 2 so second index of a zmat
            #hline = ['H',2,0.9,3,109.5,4,126.0]
            hline = ['H',2,3,4,0.9,109.5,126.0]
        else:
            bondindex=item[1] #or if item only has length 1
            #dihedralindex=item[3] 
            dihedralindex=item[2] 
            #hline= ['H',spos,0.9,bondindex,109.5,item[3],126.0]#LOOK AT THIS
            hline= ['H',spos,bondindex,dihedralindex,0.9,109.5,126.0]#LOOK AT THIS
    else:
        for item in zmatnew[fileparameters['ncore']:]:
            print "in loop", "spos:",spos,"str(item[1]",str(item[1])
            if str(item[1]) == spos: # so when at the atom that is connected to that position
                if item[0] == 'H': h=1 # check if that atom is an hydrogen atom
                hline=item[:] # copy that line
                item[6] = '126.0' #last position is dihedral. change to 126 for the atom 
                hline[6] = '234.0' #last position is dihedral. change to 234 for new H atom
                hline[0] = 'H' #and make that a hydrogen
                print "hline:",hline
    zmatnew.append(hline)      
    orca_writerAH(zmatnew,spos,index,**fileparameters) #now it is important where this will be written.
    return h

def orca_writer(zmat,index,**paras):
    filename= paras['identify'] + str(index)
    pol=0
    diff=0
    # BASISSET refactoring
    basis = paras['basisset']
    logging.info('basis std:' + basis)
    if not paras['nosub']==2:
        if paras['functional']=='B3P86':
            paras['functional']='B3P'
    if basis[0]=='6': #so if basis in form starting with 6 as 6-31G(d,p) etc
        basis = '_' + basis.replace('-','_') #change to _6_31G(d,p)
        if '(' in basis:
            b= re.split('\(',basis) #changes to ['_6_31G','(','d,p)']
            basis= b[0]
            pol = '_'+b[1].replace(',','').strip(')')
            logger.debug('basis: '+basis)
            logger.debug('polarization: '+pol)
        if '++' in basis:
            diff = '_pp'
            basis.replace('++','')
            logging.info('diffuse functions: '+diff)
        elif '+' in basis:
            diff = '_p'
            basis.replace('+','')
            logging.info('diffuse functions: '+diff)
        else:
           diff = None
    #HERE WRITING
    with open(paras['path'] + '/' + filename,'w') as fid:
        fid.write('# '+ filename + '\n')
        if paras['nosub']==2:
            fid.write('! 6-31G\n')
        else:
            fid.write(paras['orcaline'] + '\n')
        #base
        fid.write('%base \"'+ filename +'\"\n')
        #mulliken
        fid.write('%output\n  Print[ P_Mulliken ] 1\nend\n')
        #nosub part set maxiter 0
        if paras['nosub']==2:
            fid.write('%scf\nmaxiter 1\nguess patom\nend\n')
        else:
            #basisset
            fid.write('%basis\n  Basis ' + basis+'\n')
            if pol: fid.write('  Pol ' + pol + '\n')
            if diff: fid.write('  Diff ' + diff + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n')
        # charge multiplicity
        fid.write('* int '+str(paras['charge'])+ ' ' + str(paras['mult']) + '\n')
        #------THIS is still Gaussian: -------
        for i in range(len(zmat)):
            for j in [0,1,3,5,2,4,6]:
                try:
                    fid.writelines("%s " % zmat[i][j])
                except IndexError:
                    if j in [1,3,5]:
                        fid.writelines("0 ")
                    else:
                        fid.writelines("0.0 ")
            fid.write("\n")
        #-----------
        fid.write("*\n\n")
        if not paras['multiplejobs']==0 and not paras['nosub'] ==2:
            if paras['polar'] == 1:
                fid.write('$new_job\n')
                fid.write('! dft\n') #orcaline2
                fid.write('%base \"'+ filename +'\"\n')
                #polarizability
                fid.write('%elprop Polar 1\n end\n')
                #basisset
                fid.write('%basis\n  Basis ' + basis+'\n')
                if pol: fid.write('  Pol ' + pol + '\n')
                if diff: fid.write('  Diff ' + diff + '\n')
                fid.write('end\n')
                #dft functional
                fid.write('%method\n  functional '+paras['functional'] + '\n')
                if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                    fid.write('  LDAOpt C_VWN3\n')
                fid.write('end\n')
                # charge multiplicity
                fid.write('* xyzfile 1 2\n\n')
            if paras['ip'] == 1:
                fid.write('$new_job\n')
                fid.write('! dft\n') #orcaline2
                fid.write('%base \"'+ filename +'\"\n')
                #basisset
                fid.write('%basis\n  Basis ' + basis+'\n')
                if pol: fid.write('  Pol ' + pol + '\n')
                if diff: fid.write('  Diff ' + diff + '\n')
                fid.write('end\n')
                #dft functional
                fid.write('%method\n  functional '+paras['functional'] + '\n')
                if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                    fid.write('  LDAOpt C_VWN3\n')
                fid.write('end\n')
                # charge multiplicity
                fid.write('* xyzfile 1 2\n\n')
            if paras['ea'] == 1:
                fid.write('$new_job\n')
                fid.write('! dft\n') #orcaline2
                fid.write('%base \"'+ filename +'\"\n')
                #basisset
                fid.write('%basis\n  Basis ' + basis+'\n')
                if pol: fid.write('  Pol ' + pol + '\n')
                if diff: fid.write('  Diff ' + diff + '\n')
                fid.write('end\n')
                #dft functional
                fid.write('%method\n  functional '+paras['functional'] + '\n')
                if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                    fid.write('  LDAOpt C_VWN3\n')
                fid.write('end\n')
                # charge multiplicity
                fid.write('* xyzfile -1 2\n\n')
    print "---- FILE PRINTED SUCCESFULLY -----"
    return
                
def basissplit(basis):
    if basis[0]=='6': #so if basis in form starting with 6 as 6-31G(d,p) etc
        basis = '_' + basis.replace('-','_') #change to _6_31G(d,p)
        if '(' in basis:
            b= re.split('\(',basis) #changes to ['_6_31G','(','d,p)']
            basis= b[0]
            pol = '_'+b[1].replace(',','').strip(')')
            logger.debug('basis: '+basis)
            logger.debug('polarization: '+pol)
        else:
            pol = None
        if '++' in basis:
            diff = '_pp'
            basis = basis.replace('++','')
            logging.info('diffuse functions: '+diff)
        elif '+' in basis:
            diff = '_p'
            basis = basis.replace('+','')
            logging.info('diffuse functions: '+diff)
        else:
           diff = None
    logger.debug('basis:'+basis+'\npol:'+str(pol)+'\ndiff:'+str(diff))
    return basis,pol,diff

def orca_writerA(zmat,index,**paras):
    import re
    filename= paras['identify'] + str(index)
    if paras['functional']=='B3P86':
        paras['functional']='B3P'
    # BASISSET refactoring
    basis1,pol1,diff1 = basissplit(paras['basisset1'])
    basis2,pol2,diff2 = basissplit(paras['basisset2'])
    #HERE WRITING
    with open(paras['path'] + '/' + filename,'w') as fid:
        fid.write('# '+ filename + '\n')
        fid.write(paras['orcaline1'] + '\n')
        fid.write('%base \"'+ filename +'\"\n')
        fid.write('%output\n  Print[ P_Mulliken ] 1\nend\n')
        #basisset
        fid.write('%basis\n  Basis ' + basis1+'\n')
        if pol1: fid.write('  Pol ' + pol1 + '\n')
        if diff1: fid.write('  Diff ' + diff1 + '\n')
        fid.write('end\n')
        #dft functional
        fid.write('%method\n  functional '+paras['functional'] + '\n')
        if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
            fid.write('  LDAOpt C_VWN3\n')
        fid.write('end\n')
        # charge multiplicity
        fid.write('* int 0 2\n')
        #------THIS is not still Gaussian: -------
        for i in range(len(zmat)):
            for j in [0,1,3,5,2,4,6]:
                try:
                    fid.writelines("%s " % zmat[i][j])
                except IndexError:
                    if j in [1,3,5]:
                        fid.writelines("0 ")
                    else:
                        fid.writelines("0.0 ")
            fid.write("\n")
                

        #-----------
        fid.write("*\n\n")
        if True:
        # JOB2 B3P86 6-311+G(d,p) 0 2
            fid.write('$new_job\n')
            fid.write('! dft\n') #orcaline2
            fid.write('%base \"'+ filename +'\"\n')
            #basisset
            fid.write('%basis\n  Basis ' + basis2 +'\n')
            if pol2: fid.write('  Pol ' + pol2 + '\n')
            if diff2: fid.write('  Diff ' + diff2 + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional2'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n\n')
            # charge multiplicity
            fid.write('* xyzfile 0 2\n\n')
        # JOB3 B3LYP 6-311+G(d,p) 0 2
            fid.write('$new_job\n')
            fid.write('! dft\n') #orcaline2
            fid.write('%base \"'+ filename +'\"\n')
            #basisset
            fid.write('%basis\n  Basis ' + basis2 +'\n')
            if pol2: fid.write('  Pol ' + pol2 + '\n')
            if diff2: fid.write('  Diff ' + diff2 + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional1'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n')
            # charge multiplicity
            fid.write('* xyzfile 0 2\n\n')
        # JOB4
            fid.write('$new_job\n')
            fid.write('! dft\n') #orcaline2
            fid.write('%base \"'+ filename +'\"\n')
            #basisset
            fid.write('%basis\n  Basis ' + basis2 +'\n')
            if pol2: fid.write('  Pol ' + pol2 + '\n')
            if diff2: fid.write('  Diff ' + diff2 + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional1'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n')
            # charge multiplicity
            fid.write('* xyzfile 1 1\n\n')
        # JOB5 -1 1 B3LYP 6-311+G(d,p)
            fid.write('$new_job\n')
            fid.write('! dft\n') #orcaline2
            fid.write('%base \"'+ filename +'\"\n')
            #basisset
            fid.write('%basis\n  Basis ' + basis2 +'\n')
            if pol2: fid.write('  Pol ' + pol2 + '\n')
            if diff2: fid.write('  Diff ' + diff2 + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional1'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n')
            # charge multiplicity
            fid.write('* xyzfile -1 1\n\n')
    print "---- FILE PRINTED SUCCESFULLY -----"
    return

def orca_writerAH(zmat,spos,index,**paras):
    import re
    filename= paras['identify'] + str(index)+'_'+str(spos)
    for key in ['functional1','functional2']:
        if paras[key]=='B3P86':
            paras[key]='B3P'
    # BASISSET refactoring
    basis1,pol1,diff1 = basissplit(paras['basisset1'])
    basis2,pol2,diff2 = basissplit(paras['basisset2'])
    #HERE WRITING
    with open(paras['path'] + '/' +index + '/' + filename,'w') as fid:
        fid.write('# '+ filename + '\n')
        fid.write(paras['orcaline1'] + '\n')
        fid.write('%base \"'+ filename +'\"\n')
        fid.write('%output\n  Print[ P_Mulliken ] 1\nend\n')
        #basisset
        fid.write('%basis\n  Basis ' + basis1+'\n')
        if pol1: fid.write('  Pol ' + pol1 + '\n')
        if diff1: fid.write('  Diff ' + diff1 + '\n')
        fid.write('end\n')
        #dft functional
        fid.write('%method\n  functional '+paras['functional1'] + '\n')
        if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
            fid.write('  LDAOpt C_VWN3\n')
        fid.write('end\n')
        # charge multiplicity
        fid.write('* int 0 1\n')
        #------THIS is still Gaussian: -------
        for i in range(len(zmat)):
            for j in range(len(zmat[i])):
                fid.writelines("%s " % zmat[i][j])
            fid.write("\n")

        #-----------
        fid.write("*\n\n")
        # JOB7 B3P86 6-311+G(d,p) 0 1
        if True:
            fid.write('$new_job\n')
            fid.write('! dft\n') #orcaline2
            fid.write('%base \"'+ filename +'\"\n')
            #polarizability
            #fid.write('%elprop Polar 1\n end\n')
            #basisset
            fid.write('%basis\n  Basis ' + basis2 +'\n')
            if pol2: fid.write('  Pol ' + pol2 + '\n')
            if diff2: fid.write('  Diff ' + diff2 + '\n')
            fid.write('end\n')
            #dft functional
            fid.write('%method\n  functional '+paras['functional2'] + '\n')
            if paras['functional'] in ['B3LYP','b3lyp','B3P86','B3P','b3p','b3p86']:
                fid.write('  LDAOpt C_VWN3\n')
            fid.write('end\n\n')
            # charge multiplicity
            fid.write('* xyzfile 0 1\n\n')
    return
