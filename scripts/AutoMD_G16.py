#!/bin/env python3
#!!!!WIP!!!NOT FULLY TESTED OR FINISHED
#NEED TO GET RID OF ID_guess in MD folder w/ chmod777 dependency!! -> see if can change to ID_Gauss in CALC
#TEMP REMOVED RUNNINGYET COMPLETED
"""

Created on Thu Oct 31 20:03:09 2021
Last Modified on Thu Mar 24 11:39:09 2022


@author: davidbartsmets

"""

#-------------------------------FILL-IN-THESE-VARIABLES--------------------------------------------
InputFileSuffix="_AutoMD"#Put what you want the inputfile name to end on in here excluding '.gjf'
#Calcline='# opt=(calcfc,tight) freq blyp def2svp empiricaldispersion=gd3bj scf=(xqc,novaracc,noincfock)'
#Specified_Memory='30GB'#Put specified memory including unit ('5GB','2000MB',...) in this string or leave empty ('') for 4GB default
#Specified_Cores='8'#Put specified cores ('2','3',...) in this string or leave empty ('') for 4 cores default 
#--------------------------------------------------------------------------------------------------

import datetime
import os
import cclib
import numpy as np
import subprocess
import time

def MakeG16MDGjf(JobName,CalculationDetailsLine,Memory_given,Cores_given,DebugPrint=False):
    #here=os. getcwd()
    now = datetime.datetime.now()
    SubmissionScript='./ID_gauss'
    NewFileSuffix=InputFileSuffix+".gjf"

    Imput = JobName+'.log' #
    print('pwd: ',os.getcwd())
    os.chdir("./CALC/")
    calcpath=os.getcwd()
    try:
        os.mkdir('../ImagLogs')
    except:
        None
    #f=open(Imput)#Opens Logfile and extracts coordinates, charge, and multiplicity
    parser = cclib.io.ccopen(Imput)
    data = parser.parse()
    HasFreqs=True #assume log has freqs
    ImagFreqs=0
    try:
        for i in data.vibfreqs:
            if i < 0:
                ImagFreqs=ImagFreqs+1
        print('#ImagFreqs='+str(ImagFreqs))
        for i in range(ImagFreqs):
            print('Adding MD for #Imag freq: '+str(i+1))
            summed=np.add(data.vibdisps[i], data.atomcoords[-1])#add manual displacement to imag freq
    except:
        HasFreqs=False
        #outp = open("LogList.txt",'a')#Starts Writing Outputfile
        #outp.write(Imput+"\n")
        #outp.close
        if DebugPrint:
            print('Skipping '+Imput+' because it does not contain frequencies')

    try:
        os.mkdir('./MDs')
    except:
        None
    if HasFreqs and data.vibfreqs[0]<0 and (not Runningyet(Imput[:-4]+NewFileSuffix)):
        os.system('cp '+Imput+' ../ImagLogs/'+Imput)#ImagLogs/ 
        #print('starting file write')
        outp = open("./MDs/"+Imput[:-4]+NewFileSuffix,'w')#Starts Writing Outputfile
        #outp.write('%rwf='+Imput[:-4]+'.rwf \n')
        #outp.write("%NoSave \n")
        outp.write('%chk='+Imput[:-4]+InputFileSuffix+'.chk \n')
        outp.write('%nprocshared='+Cores_given+' \n')
        outp.write('%mem='+Memory_given+' \n')
        outp.write(CalculationDetailsLine+'\n')
        outp.write('\n'+ Imput[:-4]+NewFileSuffix+' Inputfile made on: '+now.strftime("%Y-%m-%d %H:%M:%S")+' \n')
        outp.write('\n '+str(data.charge)+' '+str(data.mult)+' \n')                
        for i in range(len(data.atomnos)):
            outp.write (str(data.atomnos[i])+' '+"{:f}".format(summed[i][0])+' '+"{:f}".format(summed[i][1])+' '+"{:f}".format(summed[i][2])+'\n')
        outp.write("\n\n\n")
        outp.close()
        print('ending file write')
        os.chdir("./MDs/")
        #os.system('DFTBA_here '+Imput[:-4]+NewFileSuffix+' '+'11:00')
        try:
            jobsub = subprocess.check_output([SubmissionScript, "./MDs/"+Imput[:-4]+NewFileSuffix,'11:00'],cwd=calcpath)
        except:
            print("landed in expect of job submission: script might not be executable")
            print("trying to make script executable")
            scriptpath = "{}/{}".format(calcpath, SubmissionScript)
            st = os.stat(scriptpath)
            print("permissions of {} is:".format(scriptpath), st)
            import stat
            os.chmod(scriptpath, st.st_mode | stat.S_IEXEC)
            print("attempting submission again")
            jobsub = subprocess.check_output([SubmissionScript, "./MDs/"+Imput[:-4]+NewFileSuffix,'11:00'],cwd=calcpath)#!!!WIP: fixing ID_gauss dependancy!!!       
        print('Output from job submission: ',jobsub)
        WaitTillDoneThenReplace(Imput[:-4]+NewFileSuffix)
        #os.system('python3 WaitForMD.py -i '+Imput[:-4]+NewFileSuffix+' &')
        os.chdir("..")
        
    elif HasFreqs and data.vibfreqs[0]>0:
        #outp = open("LogList.txt",'a')#Starts Writing Outputfile
        print(Imput+"\n")
        #outp.close
        if DebugPrint:
            print('No imag freq for '+Imput)
    elif Runningyet(Imput[:-4]+NewFileSuffix):
        print('Still/Already running/completed MD calc: '+Imput[:-4]+NewFileSuffix)
    
    os.chdir("..")

def Runningyet(calculation):#Retruns True if requested calculation is runnung or has already completed recently
    qsta = str(subprocess.check_output(['qsta'])).split('\\n')
    for i in range(len(qsta)-1):
        if qsta[i].split()[4]==calculation:
            if qsta[i].split()[2]=="R" or qsta[i].split()[2]=="Q" or qsta[i].split()[2]=="C":
                return True
            else:
                print("unxepected jobstate in RunningYet -> returns False: ",qsta[i])
    return False
 
def WaitTillDoneThenReplace(calcname):
    print('in WaitTillDoneThenReplace begin')
    Done = False
    while not(Done):
        now = datetime.datetime.now()
        qsta = str(subprocess.check_output(['qsta'])).split('\\n')
        for i in range(len(qsta)-1):
            if qsta[i].split()[4] == './MDs/'+calcname:
                if qsta[i].split()[2] == 'C':
                    Done = True
                elif qsta[i].split()[2] == 'R' or qsta[i].split()[2] == 'Q':
                    Done = False
                else:
                    print("unxepected jobstate in WaitTillDoneThenReplace -> continues as if not encountered: ",qsta[i])
        #with open('MD_waitlog.txt', 'a') as f:
        print(now.strftime("%Y-%m-%d %H:%M:%S")+': ' +
            calcname+': still waiting on calculation')
        if not(Done):
            time.sleep(300)
    if Done:
        print('in WaitTillDoneThenReplace: calc defined as Done')
        try:
            Imput = calcname[:-4]+'.log'
            parser = cclib.io.ccopen(Imput)
            data = parser.parse()
            a=data.vibfreqs[0]
            HasFreqs=True
        except:
            HasFreqs=False
        Imput = calcname[:-4]+'.log'
        parser = cclib.io.ccopen(Imput)
        data = parser.parse()
        if HasFreqs and data.vibfreqs[0] > 0:
            #with open('MD_waitlog.txt', 'a') as f:
            print(now.strftime("%Y-%m-%d %H:%M:%S")+': ' +
                calcname + ': calculation Done and no more imag freq. Copying log file to ../'+calcname[:-11]+'.log \n')
            os.system('cp '+calcname[:-4]+'.log'+' ../'+calcname[:-11]+'.log')
        elif HasFreqs:
            #with open('MD_waitlog.txt', 'a') as f:
            print(now.strftime("%Y-%m-%d %H:%M:%S")+': ' +
                    calcname+': calculation Done, but still imag freq for '+Imput+'\n Copying Anyway, but NEEDS MANUAL ACTION! (or will perform 2nd Manual Dosplacement after +- 1day, but manual action is strongly recommended)')
            os.system('cp '+calcname[:-4]+'.log'+' ../'+calcname[:-11]+'.log')
        elif not HasFreqs:
            #with open('MD_waitlog.txt', 'a') as f:
            print(now.strftime("%Y-%m-%d %H:%M:%S")+': ' + calcname+': calculation Error, done but does not contain ANY freqs for '+Imput+'. (DISABLED) Trying to restart '+calcname+' \n')

 

"""
    log_files = [f for f in os.listdir('.') if f.endswith('.log')]
    
    try:
        f=open("LogList.txt")
    except:
        f=open("LogList.txt",'a+')
        f.close
        f=open("LogList.txt")
    
    for line in f:
        try:
            log_files.remove(line.strip())
            #outp = open("TestLog.txt",'a')#Starts Writing Outputfile
            #outp.write(line.strip()+"removed, should be same as \n")
            #outp.close

        except:
            None
    #outp = open("TestLog.txt",'a')#Starts Writing Outputfile
    print('Log Files That Will Be Checked '+str(log_files))
    #outp.close
    for i in range(len(log_files)):
"""   