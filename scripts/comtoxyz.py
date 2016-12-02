#!/bin/env python
''' this one converts a .com file to a .com.xyz file to be opened by MOLDEN par example'''
import sys
filename = sys.argv[1]
import re
print "filename:",filename

multcharge = re.compile('^[0-9]\s[0-9]') #only set the compiler
try:
  with open(filename) as fid: #again open as fid
    for line in fid:
        if multcharge.match(line): #from where there is a match it reads the subsequant lines as the zmat
            print "match!"
            print multcharge.match(line).group()
            xyzs=[]
            line=next(fid)
            while not line == '\n': #until empty line
                xyzs.append(line.split())
                line=next(fid)
            break #so that only the first match is used. after the other matches there is no zmat
except StopIteration:
  pass
print "xyzs", xyzs

with open(filename + '.xyz','w') as fid:
    fid.write("{:4d}".format(len(xyzs)))
    fid.write("\n\n")
    #other rows
    for j in range(len(xyzs)):
        fid.write( ' '.join( [ item for item in xyzs[j] ]) + '\n' )
    fid.write('\n')
    
    



