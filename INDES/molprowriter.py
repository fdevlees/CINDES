
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
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    #fid.write("%chk=" + identify + str(index) + ".chk\n")
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
    if paras['polar'] == 1:
        fid.write("--link1--\n")
        fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
        fid.write("%mem=1500MB\n")
        if not paras['nprocs']==1:
            fid.write("%nprocshared="+str(paras['nprocs'])+"\n")
        #fid.write("%nprocshared=2\n")
        fid.write(paras['gaussianlinepolar'])
        fid.write("\n")
        fid.write(str(index) + " 2nd calc\n")
        fid.write("\n")
    if paras['ip'] == 1 or paras['nosub']==1:
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
    if paras['ea'] == 1:
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
    print "---- FILE PRINTED SUCCESFULLY -----"
    return
#----- END FILEWRITER ----#
#----- BEGIN FILEWRITER JOB TYPE A BDE-MODEL ----#
def filewriterA(zmat,index,**paras): #paras is short for fileparameters
    paras['charge']=0
    paras['mult']=1
    filename = paras['identify'] + str(index) + ".com"
    fid=open(paras['path'] + '/' + filename,'w')
    fid.write("%chk=" + paras['identify'] + str(index) + ".chk\n")
    fid.write("%mem=1500MB\n")
    fid.write("%nprocshared=2\n")
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
    fid.write("%nproc=2\n")
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
    fid.write("%nprocshared=2\n")
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
    fid.write("%nprocshared=2\n")
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
    fid.write("%nprocshared=2\n")
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
    fid.write("%nprocshared=2\n")
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
