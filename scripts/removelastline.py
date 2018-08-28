#!/bin/env python
""" this script checks if file ends with IGNORED and if so removes it from the file"""
import os

# get the logfiles
logfiles = [ f for f in os.listdir('.') if '.log' in f and not 'zzz' in f ]

for logfile in logfiles:
    print "logfile:", logfile

    # open in read mode
    with open(logfile, 'r') as f:

        # go to EOF
        f.seek(0, os.SEEK_END)
        nlines = f.tell()

        # check if not almost empty file
        if nlines<5: continue
        print "nlines:", nlines

        # read backwards until I find a newline character
        chars=[]
        i=0
        while True:
            f.seek(nlines-i)
            char = f.read(1)
            if char=='\n':
                break
            else:
                chars.append(char)
            i+=1

        # end make a string of it
        lastpartoffile = "".join(chars[::-1])
        print "last part of file:", lastpartoffile

        # test if 'IGNORED' in string
        if 'IGNORED' in lastpartoffile:
            print "IGNORED in ", logfile
            # now we need the whole file
            f.seek(0)
            totalfile = f.readlines()
        else:
            continue

    # and write the whole file except the last line
    with open(logfile, 'w') as f:
        f.writelines(totalfile[:-1])



