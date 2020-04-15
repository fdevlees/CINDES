#!/bin/env python
import argparse
parser = argparse.ArgumentParser(description="converts changed out.log to cyclesinfo format")
parser.add_argument('outputfile', nargs='?', default='cyclestest',help='filename default is "cyclestest"')
parser.add_argument('inputfile', nargs='?', default='zzz3',help='filename default is "zzz3"')
args=parser.parse_args()

out = open(args.outputfile,'w')

newfile = []
maxcount = 100

with open(args.inputfile,'r') as fid:
   line = next(fid)
   print("firstline:", line)
   while True:  #maxcount
        assert 'COUNT' in line, line # first line
        counter = line.split()[1]
        if int(counter) >= maxcount: break
        print("counter:", counter)
        line = next(fid) #go to next line
        while True:
            if 'COUNT' in line: break
            assert 'k=' in line, line
            splitline = line.split()
            (k,l) = (splitline[1],splitline[3])
            line = next(fid)
            while True:
                if 'k=' in line or 'COUNT' in line: break
                assert not 'COUNT' in line, line #no C,O,N are also in the word 'COUNT'
                splitline = line.split()
                confline = '_'.join(splitline[:-1]) + ' ' + splitline[-1]
                printline =  confline + ' {} {} {}\n'.format(counter,k,l) 
                out.write(printline)
                try:
                    line = next(fid)
                except StopIteration:
                    raise SystemExit('End of File')
