#!/bin/env python
''' this module contains all functions related to the Gaussian09 program '''
import re

def writegeom(mol, fid):
    if hasattr(mol, 'zmat'):
        # here the zmat
        for i in range(len(zmat)):
            for item in zmat[i]:
                fid.writelines("%s " % item)
            fid.write("\n")
        fid.write("\n")
    elif hasattr(mol, 'xyz'):
        fid.write(mol.xyz)
        fid.write("\n")
    else:
        raise AttributeError
    return



#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(mol, **paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    #------------
    # this function uses globals: identify, path
    #------------
    index = mol.index
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')

    # JOB 1
    job1 = paras['jobs'][0]
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    if not paras['nprocs']==1:
        fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
    fid.write(job1['hotline']) # first gaussianline
    fid.write("\n\n")
    fid.write(str(mol) + "\n\n")
    fid.write("{} {}\n".format(job1['charge'], job1['mult']))

    # write geom:
    writegeom(mol, fid)

    # THE OTHER JOBS
    #for i, (charge, mult, line) in enumerate(paras['gaussianlines'][1:]):
    for i, job in enumerate(paras['jobs'][1:]):
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        fid.write(job['hotline'])
        fid.write("\n\n")
        fid.write(str(index) + " {}th calc\n\n".format(i+2) )
        if not 'allcheck' in job['hotline']:
            fid.write("{} {}\n\n".format(job['charge'], job['mult']))
    fid.close()
    return

def filewriterAH(zmat, pos, index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    #------------
    # this function uses globals: identify, path
    #------------
    filename = paras['identify'] + str(index) + "_{}.com".format(str(pos))
    fid=open(paras['path'] + '/' + index + '/' + filename,'w')

    # JOB 1
    job1 = paras['stabjobs'][0]
    fid.write("%chk=" + paras['identify'] + str(index) + "_" + str(pos) + ".chk\n")
    fid.write("%mem=1500MB\n")
    if not paras['nprocs']==1:
        fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
    fid.write(job1['hotline']) # first gaussianline
    fid.write("\n\n")
    fid.write(paras['identify'] + str(index) + "\n\n")
    fid.write("{} {}\n".format(job1['charge'], job1['mult']))
    # here the zmat
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")

    # THE OTHER JOBS
    #for i, (charge, mult, line) in enumerate(paras['gaussianlines'][1:]):
    for i, job in enumerate(paras['stabjobs'][1:]):
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + "_" + str(pos) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        fid.write(job['hotline'])
        fid.write("\n\n")
        fid.write(str(index) + " {}th calc\n\n".format(i+2) )
        if not 'allcheck' in job['hotline']:
            fid.write("{} {}\n\n".format(job['charge'], job['mult']))
    fid.close()
    return

def get_paths( mols, fileparameters):
    ''' get all paths that need to be examined later '''
    files=[]
    for molecule in mols:
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
    import re
    import time
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
        for attempt in range(3):
            try:
                if not termination(path)==1:
                    print "Error termination:",path
                    bnewfile = errortermination(path,debug)
            except IOError as e:
                time.sleep(10)
            else:
                break
        else:
            raise e

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

    print "mols_toread:", mols_toread
    return mols_toread

def errortermination(path,debug=False):
    import re
    import time

    from CINDES4.cclib.parser.gaussianparser import Gaussian
    mymol = Gaussian(path).parse()
    from CINDES4.utils import utils
    t=utils.PeriodicTable()
    if hasattr(mymol,'atomcoords'):
        coords = map(list, mymol.atomcoords[-1])
        print coords
        for sym,xyz in zip(mymol.atomnos,coords):
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
                            #xyz = mymol.atomcoords[-1]
                            xyz_f = [ item[0] + ' '.join( map( "{:12.6f}".format, item[1:])) + '\n' for item in coords ]
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






