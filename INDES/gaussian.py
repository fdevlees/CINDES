#!/bin/env python
''' this module contains all functions related to the Gaussian09 program '''
import re


#---- START FILEWRITER2 THIS ONLY FOR MAKERS TRY TO MAKE THIS ONE UNIVERSAL ----#
def filewriter2(zmat,index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat 
    The name of the file contains the index in the name
    Still a lot of hardcoded information:
        - charge
        - multiplicity
        - memory
        - number of processors
        - calculation procedure:
            > geometry optimization
            > unrestricted DFT - B3LYP functional
            > basisset: 6-31G(d)    '''
    #------------
    # this function uses globals: identify, path, gaussianline, gaussianline2, twojob
    # maybe something like:
    # filewriter(zmat,*args,**kwargs):
    # and then calling it with
    # filewriter(zmat, gaussianline, gaussianline2, path=path, identify=identify, charge=0, mult=1)
    # or with
    # filewriter(zmat, **filedic)
    # with filedic is:
    # filedic = {"charge":0,"mult":1,"identify":identify}
    #------------
    if paras['jobs']:
        filewriter4(zmat, index, **paras)
        return
    else:
        print "WARNING: you use a deprecated functionality. use jobs keyword for up-to-date program"

    if paras['gaussianlines']:
        filewriter3(zmat, index, **paras)
        return


    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    if not paras['nprocs']==1:
        fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
    fid.write(paras['gaussianline'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    # here the zmat
    # core
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    if paras['solv']:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        fid.write(paras['gaussianline_solv0'])
        fid.write("\n")
        fid.write(str(index) + " no solvent calc\n")
        fid.write("\n")
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        fid.write(paras['gaussianline_solv1'])
        fid.write("\n")
        fid.write(str(index) + " with solvent calc\n")
        fid.write("\n")
    elif paras['polar'] == 1 and paras['multiplejobs']>=1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
    elif paras['ip'] == 1 or paras['aip']==1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
        fid.write(str(paras['charge']+1) + " " + str(paras['mult']+1) + "\n")
        fid.write("\n")
    if paras['ea'] == 1 or paras['aea']==1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
        fid.write(str(paras['charge']-1) + " " + str(paras['mult']+1) + "\n")
        fid.write("\n")
    if paras['twojob']>=2:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianline2'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
        fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
        fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
    return
#----- END FILEWRITER2 ----#

def filewriter3(zmat,index,**paras): #paras is short for fileparameters
    '''    This function creates a file with the geometry contained in zmat 
    The name of the file contains the index in the name
    Still a lot of hardcoded information:
        - charge
        - multiplicity
        - memory
        - number of processors
        - calculation procedure:
            > geometry optimization
            > unrestricted DFT - B3LYP functional
            > basisset: 6-31G(d)    '''
    #------------
    # this function uses globals: identify, path, gaussianline, gaussianline2, twojob
    # maybe something like:
    # filewriter(zmat,*args,**kwargs):
    # and then calling it with
    # filewriter(zmat, gaussianline, gaussianline2, path=path, identify=identify, charge=0, mult=1)
    # or with
    # filewriter(zmat, **filedic)
    # with filedic is:
    # filedic = {"charge":0,"mult":1,"identify":identify}
    #------------
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')

    # JOB 1
    charge, mult, line = paras['gaussianlines'][0]
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    if not paras['nprocs']==1:
        fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
    fid.write(line) # first gaussianline
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    #fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    fid.write(str(charge) + " " + str(mult) + "\n")
    # here the zmat
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")

    # THE OTHER JOBS
    for i, (charge, mult, line) in enumerate(paras['gaussianlines'][1:]):
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        fid.write(line)
        fid.write("\n")
        fid.write(str(index) + " {}th calc\n".format(i+2) )
        fid.write("\n")
        if not 'allcheck' in line:
            fid.write(str(charge) + " " + str(mult) + "\n")
            fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
    return

def filewriter4(zmat,index,**paras): #paras is short for fileparameters
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
    return

#----- BEGIN FILEWRITER JOB TYPE A BDE-MODEL ----#
def filewriterA(zmat,index,**paras): #paras is short for fileparameters
    paras['charge']=0
    paras['mult']=1
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    # gaussianline now is:
    # "# opt=calcfc ub3lyp/6-31g(d) pop=npa freq'
    fid.write(paras['gaussianline1'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    # here the zmat
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    # JOB2 larger basis for E(A) to calculate energy for E(A)BDE
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    #fid.write("%nproc={:d}\n".format(paras['nprocs'])
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    #gausline3 is:
    # '# geom=check guess=read b3p86/6-311+G(d,p)'
    fid.write(paras['gaussianline3'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    fid.write("\n")
    # JOB3 E0 E(A)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    #gausline2 is:
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']+1) + "\n")
    fid.write("\n")
    # JOB4 IP E(A+)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    #gausline4 is same as gausline2:
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']+1) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    # JOB5 EA E(A-)
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    #gausline5 is same as gausline2
    # '# geom=check guess=read b3lyp/6-311+G(d,p)'
    fid.write(paras['gaussianline2'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']-1) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
    fid.close()
    return

def filewriterAH(zmat,pos,index,**paras): #paras is short for fileparameters
    paras['charge']=0
    paras['mult']=1
    filename = paras['identify'] + str(index) + "_" + pos + ".com"
    print "in filewriter AH; filename:",filename
    fid=open(paras['path'] + '/' + index + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + "_" + pos + ".chk\n")
    fid.write("%mem=1500MB\n")
    #fid.write("%nprocshared=2\n")
    fid.write("%nprocshared={:d}\n".format(paras['nprocs']) )
    fid.write(paras['gaussianline1'])
    fid.write("\n")
    fid.write(paras['identify'] + str(index) + "\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    # zmat:
    for i in range(len(zmat)):
        for item in zmat[i]:
            fid.writelines("%s " % item)
        fid.write("\n")
    fid.write("\n")
    fid.write("--link1--\n")
    fid.write("%chk=" + paras['identify'] + str(index) + "_" + pos + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
    fid.write(paras['gaussianline3'])
    fid.write("\n")
    fid.write(str(index) + " 2nd calc\n")
    fid.write("\n")
    fid.write(str(paras['charge']) + " " + str(paras['mult']) + "\n")
    fid.write("\n")
    #print "---- FILE PRINTED SUCCESFULLY -----"
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
