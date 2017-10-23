#!/bin/env python
import sys
import pickle
from pprint import pprint

class args:
    column=3


class Tablebin(object):
    ''' class for tablebin binary pickled data files '''

    def __init__(self, filename='tablebin', column=1, data=None):
        if data:
            self.data = data
        else:
            self.data = read_file(filename)
        self.set_confs()
        self.set_values(column)
        return

    def __str__(self):
        return "Tablebin object with:\n{} elements".format(len(self.data))

    def set_confs(self, factor=1.0, max_value=None):
        ''' get configurations '''
        self.indices = [ item[0] for item in self.data ]
        self.confs = [ item.split('_') for item in self.indices ]
        return

    def set_values(self, column, ftype=float):
        self.Y = [ ftype(item[column]) for item in self.data ]
        return

    def write(self, filename='tablebin_new'):
        data_out = [ line for line in map(list,zip(self.indices, self.Y)) ]
        #print "data_out:", data_out
        with open(filename,'wb') as fout:
            pickle.dump(data_out,fout)
        print "new file written with name {} with {} lines".format(filename,str(len(data_out)))
        return

    def print_table(self):
        for item in self.confs: print item
        for i, (confje, y) in enumerate(zip(self.confs, self.Y)):
            conf = ' '.join( [ '{:9s}'.format(group) for group in confje ] )
            data = '{:15.8f}'.format(y)
            print '{:4} {} {}'.format(i,conf,data)
        return

    def print_table_2(self, column=[1]):
        if column:
            maxlen = max( [ len(item[0]) for item in sortev[::-1] ] )
            for item in sortev[::-1]:
                print  '{ind:{width}}{data}'.format(width= maxlen,ind = item[0], data = ' '.join( [ '{:15.8}'.format(item[int(i)]) for i in args.column ] ) )
        else:
            for item in sortev[::-1]:
                print '{} {:.6}'.format(item[0],item[1])

    def get_seq(self):
        '''seq is sequence of all types of functional groups and dopants present'''
        seq = set()
        for conf in self.confs:
            for group in conf:
                if not group in seq:
                    seq.add(group)
        # now make sure the order is first tert than second:
        seq = list(seq)
        for secondary_group in ['CO','O','S']:
            if secondary_group in seq:
                seq.remove(secondary_group)
                seq.append(secondary_group) # take it out and append it back at the end
        if True: #place CH as first element
            if 'CH' in seq:
                seq.remove('CH')
                seq.insert(0,'CH')
        return seq

    def get_ngps(self):
        '''ngps = tuple of number of groups per site'''
        #assuming that every conf has the same number of sites:
        nsites = len(self.confs[0])
        gps= [ set() for _ in range(nsites) ]
        for conf in self.confs:
            for i,group in enumerate(conf):
                if not group in gps[i]:
                    gps[i].add(group)
        #print "gps:", gps
        ngps = map(len,gps)
        return ngps

    def diversity_filter(self, n=10, divindex=1):
        from CINDES4.utils.diversity1 import Diversifier
        self.diversifier = Diversifier(index=divindex)
        self.diversifier.set_filter(False)
        divalues, occupancy = self.diversifier.table_run(self, index=divindex)

        # if max number is given:
        if n:
            # to order simultaneously
            l = sorted(zip(divalues, self.confs, self.indices, self.Y))
            # get the firt n elements of all the lists
            divalues, self.confs, self.indices, self.Y = zip(*l[:n])
            if args.verbose>=1:
                print "5 most diverse samples"
                for i in xrange(5):print divalues[i], self.indices[i], self.Y[i]
        return






def read_file(filename):
    data = []
    with open(filename,'rb') as f:
        try:
            data.append(pickle.load(f))
        except EOFError:
            print "empty tablebin?"
            raise
    return data[-1]

def main(args):
    table = Tablebin(args.filename)

    if args.formatted:
        table.print_table()
    print "jos"
    if args.diversify:
        table.diversity_filter(divindex=args.diversify, n=args.nmax)

    if args.new_table:
        table.write(args.new_table)

    return


#if args.plot:
#    SCFs = [ row[1] for row in confs ]
#    import matplotlib.pyplot as plt
#    plt.plot(SCFs,'ro')
#    plt.grid(True)
#    plt.show()


def get_homo(filename,identify='ada_'):
    filetje = path + '/databc/' + identify + filename + '.log'
    #filetje = path + '/s15/dia_' + filename + '.log'
    f = ccopen(filetje)
    f.logger.setLevel(logging.ERROR)
    datatje = f.parse()
    HOMO = datatje.myhomos[1]
    #datatje.mymos[0]['alpha'][0][HOMO] #0 is guess orbital energies
    Ehomo=datatje.mymos[1]['alpha'][0][HOMO] #1 is after first optimization
    Elumo=datatje.mymos[1]['alpha'][0][HOMO+1]
    if 0:
        print "filename, Ehomo:", filename, Ehomo
    elif 1:
        print "filename, Ehomo, Elumo, Egap:", filename, Ehomo, Elumo, Elumo-Ehomo
    else:
        sys.stdout.write('#')
    return Ehomo

class Unbuffered(object):
    def __init__(self,stream):
        self.stream = stream
    def write(self,data):
        self.stream.write(data)
        self.stream.flush()
    def __getattr__(self,attr):
        return getattr(self.stream, attr)
sys.stdout = Unbuffered(sys.stdout)

if False:
    if args.regxy: # plots two properties vs each other. and fits a straight line through it
        from cclib.parser import ccopen
        import logging
        import os
        import matplotlib.pyplot as plt
        import fnmatch
        path = os.getcwd()
        from scipy import stats
        propx = [ item[1]*27.21138 for item in sortev ] 
        propy = [ -1*get_homo(item[0],identify=args.identify) for item in sortev ]
        plt.plot(propx,propy,'.r')
        slope, intersept, r_value, p_value, std_err = stats.linregress(propx,propy)
        print "slope:", slope
        print "intersept:", intersept
        print "p_value:", p_value
        print "std_err:", std_err
        print "R=",r_value**2
        x = sorted(propx)
        y = [ slope*xje+intersept for xje in x ]
        plt.plot(x, y, '-')
        plt.xlabel('property1 IP(eV)')
        plt.ylabel('property2 -Ehomoe(eV)')
        plt.title('two properties regression')
        plt.show()


#if args.regrs:
#    #plot two properties vs each other. indicated via command line
#    propx = [ item[int(args.regrs[0])] for item in sortev ]
#    propy = [ item[int(args.regrs[1])] for item in sortev ]
#    import matplotlib.pyplot as plt
#    from scipy import stats
#    plt.plot(propx,propy,'.r')
#    slope, intersept, r_value, p_value, std_err = stats.linregress(propx,propy)
#    print "slope:", slope
#    print "intersept:", intersept
#    print "p_value:", p_value
#    print "std_err:", std_err
#    print "R=",r_value**2
#    x = sorted(propx)
#    y = [ slope*xje+intersept for xje in x ]
#    plt.plot(x, y, '-')
#    plt.xlabel('property1')
#    plt.ylabel('property2')
#    plt.title('two properties regression')
#    plt.show()





if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="unpickles data stored with pickle module")
    parser.add_argument("-p","--plot",action="store_true",help="make also a plot of the data")
    parser.add_argument("-m","--max", type=float,default = 1.0e99, help="sets a maximum for the absolute values taken into account")
    parser.add_argument("-x","--regxy",action="store_true",help="make a plot of two propeties against each other and test linear correlation")
    parser.add_argument("-N","--new_table",type=str,default='',nargs='?',const='tablebin_new',help="write new table")
    parser.add_argument("-n","--nmax",type=int,default=100,help="nconfs in new table")
    parser.add_argument("-f","--formatted",action="store_true",help="print all in formatted order")
    parser.add_argument("-F","--diversify",type=int,default=0,help="print all in formatted order")
    parser.add_argument("-i","--identify",action="store",type=str,help="identify for gethomo")
    parser.add_argument("-y","--regrs",nargs=2,help="make a plot of value1 and value2 against each other and test linear correlation")
    parser.add_argument("-c","--column",nargs='+',default = [1], help="specify which columns are printed")
    parser.add_argument("-v","--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument('filename', default='tablebin',nargs='?')
    args=parser.parse_args()
    #args_dict = vars(args)
    #print args_dict

    main(args)
