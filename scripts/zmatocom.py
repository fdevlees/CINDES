#!/bin/env python
''' this one converts a zmatrix file to a .com file'''
import sys
filename = sys.argv[1]
import re
print "filename:",filename

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

print "ZMAT:"
print zmat


with open(filename + '.com','w') as fid:
    fid.write('# opt \n\ntitle\n\n0 1\n')
    #row 1
    fid.write(zmat[0][0] + '\n')
    #row 2
    fid.write( ' '.join( [ zmat[1][i] for i in [0,1,4] ] ) + '\n' )
    #row 3
    fid.write( ' '.join( [ zmat[2][i] for i in [0,1,4,2,5] ] ) + '\n' )
    #other rows
    for j in range(3,len(zmat)):
        fid.write( ' '.join( [ zmat[j][i] for i in [0,1,4,2,5,3,6] ] ) + '\n' )
    fid.write('\n')



