# -*- coding: utf-8 -*-
"""
Created on Mon Nov 22 01:42:26 2021
Last Modified on Mon Nov 22 12:40:26 2021 

@author: davidbartsmets
"""

def get_wbo_data(logfile):#logfile=logfile starting from Wiberg Bond Order keyphrase
    Index=0
    Subindex=0
    data=[]

    Tati=True
    #with open(logfile) as f:
    for line in logfile:
        if line.startswith(" Wiberg bond index, Totals by atom:"):
            #print('Woller')
            break
        elif Tati:
            splitted=line.split()
            if len(splitted)==0 or splitted[0]=='----':
                continue
            elif splitted[0]=='Atom':
                Index=Index+1
                Subindex=0
                continue
            elif Index==1:
                Subindex=Subindex+1
                x=[]
                for i in range(len(splitted)):
                    if i>1:
                        x=x+[float(splitted[i])]
                data.append(x)
            elif Index>1:
                Subindex=Subindex+1
                x=[]
                for i in range(len(splitted)):
                    if i>1:
                        data[Subindex-1].append(float(splitted[i]))
         
    #print('data')
    #print(data)

    b = dict()
    #print(b)
    for i in range(len(data)):
        for j in range(len(data)):
            b[str(i+1)+"-"+str(j+1)]=data[i][j]
            #print(b)
    #print(b)

    return b
