# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 16:13:26 2017

@author: Tatiana
"""
import numpy as np

def getfloat(x):
    return float(x.replace('D','E'))

def skippen(c,f):
    # Skip c lines in file f
    for i in xrange(c):
        next(f)

def splitsen(f):
    # returns the float value in the second column if the next line in file f
    a=next(f).split()
    d=getfloat(a[1])
    return d

def getPosVal(f,column=1):
    # move to the next line in the file and split it at the spaces
    lineParts = next(f).split()
    # return the first part (i.e. title of the line) and the number at the specified column
    return lineParts[0], getfloat(lineParts[column])

def populateMatrix(m, pos, val):
    # define a dictionary to translate the string characters of 'pos'
    transl = {'x':0, 'y':1, 'z':2}
    # use the dictionary to translate each character in an index and add to list 'i'
    i = [transl[char] for char in pos]
    if len(i) == 2:
        # 2D matrix
        m[i[0], i[1]] = val
    elif len(i) == 3:
        # 3D matrix
        m[i[0], i[1], i[2]] = val
    return m

def copyLowerToUpperTriangle(m):
    lowerTriangle = np.tril(m, -1)
    m += np.transpose(lowerTriangle)
    return m

def constructSymmetricalMatrix(f):
    # Create a new 3x3 matrix filled with zeros
    m = np.zeros([3,3])
    # Read the next 6 lines of the file and store in matrix m
    for j in range(6):
        # position in an xy or xyz string representation and the float 'value' to store at the position of m
        position, value = getPosVal(f)
        m = populateMatrix(m, position, value)
    m = copyLowerToUpperTriangle(m)
    return m

def vperm(b):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if i==j and i==k:
                    continue
                elif i!=j and i==k:
                    if j==0 and i+j+k==2:
                        b[j,i,i]=b[1,0,1]
                        b[i,i,j]=b[1,0,1]
                    elif j==0 and i+j+k==4:
                        b[j,i,i]=b[2,0,2]
                        b[i,i,j]=b[2,0,2]
                    elif j==1 and i+j+k==5:
                        b[j,i,i]=b[2,1,2]
                        b[i,i,j]=b[2,1,2]
                elif i!=k and i==j:
                    if k==1 and i+j+k==1:
                        b[k,i,i]=b[0,0,1]
                        b[i,k,i]=b[0,0,1]
                    elif k==2 and i+j+k==2:
                        b[k,i,i]=b[0,0,2]
                        b[i,k,i]=b[0,0,2]
                    elif k==2 and i+j+k==4:
                        b[k,i,i]=b[1,1,2]
                        b[i,k,i]=b[1,1,2]
                elif i!=k and i!=j and k!=j:
                    b[i,j,k] =b[1,0,2]
    return b

def mperm(b):
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if i!=j and j!=k:
                    if i+j+k==1 and i==0 and j==1:
                        b[i,k,j]=b[0,1,0]
                    elif i+j+k==2 and i==0 and j==2:
                        b[i,k,j]=b[0,2,0]
                    elif i+j+k==4 and i==1 and k==1:
                        b[i,k,j]=b[1,2,1]
                    else:
                        b[i,j,k]=b[i,k,j]
    return b
#--------------------------------------------------
def alpha_iso(a):
   aiso=(a[0,0]+a[1,1]+a[2,2])/3
   return aiso

def alpha_aniso(a):
   aniso=np.sqrt(((a[0,0]-a[1,1])**2+(a[0,0]-a[2,2])**2+(a[2,2]-a[1,1])**2+ (6*(a[1,2]**2+ a[0,2]**2+a[0,1]**2)))/2)
   return aniso

def beta_per(b):
    bper=0
    for i in range(3):
        bper += 2*b[2,i,i] - 3*b[i,2,i] + 2*b[i,i,2]
    bper /= 5
    return bper

def beta_par(b):
   bpar=0
   for i in range(3):
        bpar += b[2,i,i] + b[i,2,i] + b[i,i,2]
   bpar /= 5
   return bpar

def t1(tw):
        ab=0
        for i in range(3):
            ab += tw[i,i,i]**2
        return ab

def t2(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[i,i,j]**2
        return ab

def t3(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[i,i,i] * tw[i,j,j]
        return ab

def t4(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[j,i,i] * tw[i,i,j]
        return ab

def t5(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[i,i,i] * tw[j,j,i]
        return ab

def t6(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[j,i,i]**2
        return ab

def t7(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and k!=j:
                        ab += tw[i,i,j] * tw[j,k,k]
        return ab

def t8(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and k!=j:
                        ab += tw[j,i,i] * tw[j,k,k]
        return ab

def t9(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and k!=j:
                        ab += tw[i,i,j] * tw[k,k,j]
        return ab

def t10(tw):
       ab=0
       for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and j!=k:
                        ab += tw[i,j,k]**2
       return ab

def t11(tw):
       ab=0
       for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and j!=k:
                        ab += tw[i,j,k] * tw[j,i,k]
       return ab

def t12(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                if j!=i:
                    ab += tw[i,j,j]**2
        return ab

def t13(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and k!=j:
                        ab += tw[i,j,j] * tw[i,k,k]
        return ab

def t14(tw):
        ab=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if j!=i and k!=i and k!=j:
                        ab += tw[i,i,k] * tw[j,j,k]
        return ab

def bzzz(tw):
    bz1 = t1(tw)/7 + (4*t2(tw) + 2*t3(tw) + 4*t4(tw) + 4*t5(tw) + t6(tw))/35 + (4*t7(tw) + t8(tw) + 4*t9(tw) + 2*t10(tw) + 4*t11(tw))/105
    return bz1

def bzxx(tw):
    bx1 = t1(tw)/35 + 4*t3(tw)/105 - 2*t5(tw)/35 + 8*t2(tw)/105 + 3*t12(tw)/35 - 2*t4(tw)/35 + t13(tw)/35 - 2*t14(tw)/105 - 2*t7(tw)/105 + 2*t10(tw)/35 - 2*t11(tw)/105
    return bx1

def get_alpha_data(f):
    a = dict()
    skippen(5,f)
    # Construct a 3x3 symmetrical matrix from the next 6 lines in the file
    a1 = constructSymmetricalMatrix(f)
    skippen(17,f) # si Gaussian est embetant et imprime tt deux fois
    # Polarizability Dynamic
    # Construct a 3x3 symmetrical matrix from the next 6 lines in the file
    a2 = constructSymmetricalMatrix(f)
    a["Static isotropic polarizability"]= alpha_iso(a1)
    a["Static anisotropy of the polarizability"]= alpha_aniso(a1)
    a["Dynamic isotropic polarizability"]=alpha_iso(a2)
    a["Dynamic anisotropy of the polarizability"]= alpha_aniso(a2)
    return a

def get_beta_data(f):
    b = dict()
    # Skip 11 lines in the file
    skippen(11,f)
    # Create a new 3x3x3 matrix filled with zeros
    b1 = np.zeros([3,3,3])
    # Read the next 10 lines of the file and store in the matrix
    for j in range(10):
        # position in an xyz string representation and the float 'value' to store at the position of the matrix
        position, value = getPosVal(f)
        b1 = populateMatrix(b1, position, value)
        b1 =vperm(b1)
    skippen(57,f)
    b2 = np.zeros([3,3,3])
    for j in range(18):
    # position in an xyz string representation and the float 'value' to store at the position of the matrix
        position, value = getPosVal(f)
        b2 = populateMatrix(b2, position, value)
        b2=mperm(b2)
    b["Static longitudinal first hyperpolarizability"]=beta_par(b1)
    b["Static transversal first hyperpolarizability"]=beta_per(b1)
    b["Static HRS first hyperpolarizability"]=np.sqrt(bzzz(b1)+bzxx(b1))
    b["Static depolarization ratio"]=bzzz(b1)/bzxx(b1)
    b["Dynamic longitudinal first hyperpolarizability"]=beta_par(b2)
    b["Dynamic transversal first hyperpolarizability"]=beta_per(b2)
    b["Dynamic HRS first hyperpolarizability"]=np.sqrt(bzzz(b2)+bzxx(b2))
    b["Dynamic depolarization ratio"]=bzzz(b2)/bzxx(b2)
    return b
