#!/bin/env python
''' this module contains all functions related to the NWChem program '''
import re
multiplicity = {1:'singlet', 2:'doublet', 3:'triplet', 4:'quartet'}

def write_subjob(fid, job):
    hotline=job['hotline']
    func, basis = [ w for w in hotline.split() if '/' in w ][0].split('/')
    if func in ['rhf','uhf']:
        theory='scf'
    elif func in ['mp2','scf','ccsd','ccsd(t)']:
        theory=func
    else:
        theory='dft'
        if func[0]=='u': #so functional is defined as ub3lyp or ub3pw91
            odft=True
            func=func[1:] # i.e. remove the u from the functional
        else: odft=False
    fid.write("charge {charge}\n".format(charge=job['charge']))
    #basisset/functional
    fid.write("basis\n * library {}\nend\n".format(basis))
    if theory=='scf':
        fid.write('scf\n {}\n maxiter 50\n {}\nend\n'.format(multiplicity[int(job['mult'])],func))
    if theory=='dft':
        fid.write('dft\n')
        fid.write(' iterations 100\n')
        if odft: fid.write(' odft\n') #makes it an open-shell unrestricted calculation
        # for quadratic convergence:
        if 'scf=xqc' in hotline: fid.write(' cgmin\n')
        if 'fukui' in hotline: fid.write(' fukui\n')
        fid.write(' mult {}\n'.format(job['mult']))
        #if func=='b3lyp':xc='vwn_1_rpa 0.19 lyp 0.81 HFexch 0.20 slater 0.80 becke88 nonlocal 0.72'
        if func=='b3lyp':xc='b3lyp'
        elif func=='b3pw91':xc='acm'
        elif func=='bp86':xc='becke88 perdew86'
        elif func=='b3p86':xc='vwn_1_rpa 1.00 perdew86 0.81 HFexch 0.20 slater 0.80 becke88 nonlocal 0.72'
        elif func=='cam-b3lyp':xc='xcamb88 1.00 lyp 0.81 vwn_5 0.19 hfexch 1.00 cam 0.33 cam_alpha 0.19 cam_beta 0.46'
        elif func=='blyp':xc='becke88 lyp'
        else: raise NameError('no valid functional')
        fid.write(' xc {}\n'.format(xc))
        fid.write('end\n')
    if 'cosmo' in hotline:
        dielec = filter(lambda x:'cosmo' in x, hotline.split())[0].split('=')[1]
        fid.write("cosmo\n {}\nend\n".format(dielec))
    if 'smd' in hotline:
        solvent = filter(lambda x:'cosmo' in x, hotline.split())[0].split('=')[1]
        fid.write("cosmo\n do_cosmo_smd true\n solvent {}\nend\n".format(solvent))
    if 'opt' in hotline:
        fid.write('driver\n maxiter 100\nend\n')
    if 'pop=npa' in hotline:
        fid.write('property\n nbofile\nend\n')

    fid.write("task {}".format(theory))
    # functional mult
    if 'opt' in hotline: fid.write(" optimize")
    if 'freq' in hotline: fid.write(" frequencies")
    fid.write("\n\n")
    return


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(zmat,index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    #------------
    # this function uses globals: identify, path
    #------------
    filename = paras['identify'] + str(index)
    fid=open(paras['path'] + '/' + filename,'w')

    # JOB 1
    # extract info
    job1 = paras['jobs'][0]

    # write info
    fid.write("echo\nstart {filename}\n".format(filename=filename))
    fid.write("memory 1500 mb\n")
    fid.write("title \"{filename}\"\n".format(filename=filename))
    # here the zmat
    fid.write("geometry\n zmatrix\n")
    for i in range(len(zmat)):
        fid.write("  ")
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write(" end\nend\n")
    write_subjob(fid, job1)

    # THE OTHER JOBS
    for i, job in enumerate(paras['jobs'][1:]):
        # write title
        fid.write("title \"{filename}\"\n".format(filename=filename))
        write_subjob(fid, job)

    fid.close()
    return

def filewriterAH(zmat, pos, index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    #------------
    # this function uses globals: identify, path
    #------------
    filename = paras['identify'] + str(index) + "_{}".format(str(pos))
    fid=open(paras['path'] + '/' + index + '/' + filename,'w')
    job1 = paras['stabjobs'][0]

    # write info
    fid.write("echo\nstart {filename}\n".format(filename=filename))
    fid.write("memory 1500 mb\n")
    fid.write("title \"{filename}\"\n".format(filename=filename))
    # here the zmat
    fid.write("geometry\n zmatrix\n")
    for i in range(len(zmat)):
        fid.write("  ")
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write(" end\nend\n")
    write_subjob(fid, job1)

    # JOB 1
    for i, job in enumerate(paras['stabjobs'][1:]):
        # write title
        fid.write("title \"{filename}\"\n".format(filename=filename))
        write_subjob(fid, job)

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
           text = fid.readlines()[-4:]
           if re.search('Zhang',''.join(text)):
               fid.close()
               return 1
           elif re.search('IGNORE',''.join(text)):
               fid.close()
               return 2
    #-----
    filepaths = get_paths( mols_tocal, fileparameters)
    copyfilepaths = filepaths[:] #copy to be able to append to it while looping over it
    for path in copyfilepaths: #test all for information which jobs crashed
        if not termination(path)==1:
            print "Error termination:",path
            bnewfile = errortermination(path, fileparameters)
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

def errortermination(path, fileparameters):
    import re
    import time

    from CINDES4.cclib.parser.nwchemparser import NWChem
    mymol = NWChem(path).parse()
    from CINDES4.utils import utils
    t=utils.PeriodicTable()
    if hasattr(mymol,'atomcoords'):
        coords = map(list, mymol.atomcoords[-1])
        #print coords
        for sym,xyz in zip(mymol.atomnos,coords):
            xyz.insert(0,t.element[sym])
        print "atomcoords and added elements:"
        for item in mymol.atomcoords[-1]:
            print ' '.join(map(str,item))
        if fileparameters['debug']==True:
            import submitter
            if hasattr(mymol,'optdone'):
                if mymol.optdone==False:
                    print "Optimizations not converged!"
                elif mymol.optdone==True:
                    print "Optimization is converged!"
            fid = open(path[:-4],'r') #change .log in .com extension and read input file
            multcharge = re.compile(' zmatrix') #a regex for the mult charge line
            newfile=[]
            once=0 #only find that line once
            for line in fid: #copy file exept for the zmat found in the inputfile
                if multcharge.match(line) and once==0: #when found 
                    print "match!"
                    once+=1
                    #newfile.append(line) #the line with the match itself has to be included in the newfile
                    while True:
                        line= next(fid) #take al new lines
                        if line==' end\n': #end of zmat
                            #now instead of this zmat that is now completely skipped place in newfile
                            #the last coordinates of the crashed run
                            #newfile.extend([' '.join( map("{12.6f}".format, item))+'\n' for item in mymol.atomcoords[-1]])
                            #xyz = mymol.atomcoords[-1]
                            print "coords:", coords
                            # filter 'X':
                            coords = filter(lambda x:not x[0] is None, coords)
                            xyz_f = [ " " + item[0] + ' '.join( map( "{:12.6f}".format, item[1:])) + '\n' for item in coords ]
                            #xyz_f = [ item[0] + ' '.join(
                            #                              map(
                            #                                   str, item[1:]
                            #                                 )
                            #                            ) + '\n' for item in xyz ]
                            print xyz_f
                            newfile.extend(xyz_f)
                            break
                else:
                    newfile.append(line) #copy that line because it is not the zmat found in the inputfile
            print "="*20
            open(path[:-4]+'zzz','w').writelines(newfile)
            print "newfile written in: ", path[:-4] + 'zzz'

            #---- preparation for submit command ---
            splitpath = path.split('/')
            filename = splitpath[-1]
            folder = '/'.join(splitpath[:-1])
            filenamesplit = filename[:-4].split('_')
            identify= filenamesplit[0]+'_'
            index = '_'.join(filenamesplit[1:])+'zzz'
            print "folder", folder
            print "index:", index
            print "identi", identify
            #----- keywords constructed so:
            submitter.submit(folder,index,identify, script='ID_NWChem')
            return True
    else:
        print "has no coords in file"
    return False






