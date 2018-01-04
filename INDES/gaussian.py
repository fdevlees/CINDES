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

