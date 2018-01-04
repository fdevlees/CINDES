#!/bin/env python
''' this module contains all functions related to the Gaussian09 program '''
import re


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter(zmat,index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat
    The name of the file contains the index in the name
    '''
    #------------
    # this function uses globals: identify, path
    #------------
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
    logging.debug(pprint.pformat(zmat))
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
            if str(item[1]) == spos:
                if item[0] == 'H': h=1
                hline=item[:]
                item[6] = '126.0'
                hline[6] = '234.0'
                hline[0] = 'H'
                #print "hline:",hline
    zmatnew.append(hline)
    filewriterAH(zmatnew,spos,index,**fileparameters) #now it is important where this will be written.
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
    filewriterAH(zmatnew,spos + '_2',index,**fileparameters)
    return
