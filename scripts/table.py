#!/bin/env python
import sys
import pickle
from pprint import pprint

class args:
    column=3


class Tablebin(object):
    ''' class for tablebin binary pickled data files '''

    def __init__(filename):
        self.data = read_file(filename)
        self.set_confs()
        return

    def set_confs(self, factor=1.0, max_value=None):
        '''sets the molecules in self.confs sorted based on their first property value

            - one can apply a cutoff value
        '''
        if not max_value is None:
            confs = [ item for item in self.data if abs(item[1]) < max_value ]  # last entry in tablebin
        else:
            confs = self.data
        self.confs = sorted(confs,key=lambda entry: entry[1]) #sort based on second column
        return

    def print_table(self):
        for i,item in enumerate(self.confs):
            conf = ' '.join( [ '{:9s}'.format(group) for group in item[0].split('_') ] )
            data = ' '.join( [ '{:15.8f}'.format(i) for i in item[1:] ] )
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

        return seq

def read_file(filename):
    data = []
    while True:
        try:
            data.append(pickle.load(args.file))
        except EOFError:
            break
    return data[-1]

def main():
    table = Tablebin(args.filename)

    if args.formatted:
        table.print_table()


    return


print "number: ", len(data[-1]) # no. of entries

if args.plot:
    SCFs = [ row[1] for row in confs ]
    import matplotlib.pyplot as plt
    plt.plot(SCFs,'ro')
    plt.grid(True)
    plt.show()


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

if True:
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


if args.regrs:
    #plot two properties vs each other. indicated via command line
    propx = [ item[int(args.regrs[0])] for item in sortev ]
    propy = [ item[int(args.regrs[1])] for item in sortev ]
    import matplotlib.pyplot as plt
    from scipy import stats
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
    plt.xlabel('property1')
    plt.ylabel('property2')
    plt.title('two properties regression')
    plt.show()





if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description="unpickles data stored with pickle module")
    parser.add_argument("-p","--plot",action="store_true",help="make also a plot of the data")
    parser.add_argument("-m","--max", type=float,default = 1.0e99, help="sets a maximum for the absolute values taken into account")
    parser.add_argument("-x","--regxy",action="store_true",help="make a plot of two propeties against each other and test linear correlation")
    parser.add_argument("-f","--formatted",action="store_true",help="print all in formatted order")
    parser.add_argument("-i","--identify",action="store",type=str,help="identify for gethomo")
    parser.add_argument("-y","--regrs",nargs=2,help="make a plot of value1 and value2 against each other and test linear correlation")
    parser.add_argument("-c","--column",nargs='+',default = [1], help="specify which columns are printed")
    parser.add_argument('file', type=argparse.FileType('rb'), default='tablebin',nargs='?')
    args=parser.parse_args()
    args_dict = vars(args)
    print args_dict

    main(**args_dict)
