#!/bin/env python

import matplotlib.pyplot as plt
import numpy as np

def read_line(line):
    row = []
    for item in line.split():
        if item in ['nan', 'None']:
            row.append(np.nan)
        else:
            try:
                item = float(item)
            except ValueError:
                pass
            row.append(item)
    return np.asarray(row)

def read_input(inputfile):
    with open(inputfile,'r') as f:
        datar = [ read_line(line) for line in f.read().splitlines() ]
    data = np.asarray(datar)
    return data

def plot1(data, ycolumns): # just plot some data columns
    for column in ycolumns:
        plt.plot(data[:,column], label=str(column))
    plt.legend()
    plt.show()

def plot2(data, ycolumns, xcolumn): # just plot some data columns vs another column
    for column in ycolumns:
        plt.plot(data[:,xcolumn], data[:,column], label=str(column))
    plt.legend()
    plt.show()

def set_nxy(n):
    xys = { '1':(1,1),
            '2':(2,1), '3':(1,3), '4':(2,2), '5':(3,2), '6':(3,2), '7':(4,2), '8':(4,2), '9':(3,3), '10':(4,3)
          }
    return xys[str(n)]

def plot3(data, ycolumns, xcolumn=None): # plot some data columns vs another column in different subplots
    if xcolumn is None:
        x = np.arange(data.shape[0])
    else:
        x = data[:, xcolumn]

    nx,ny = set_nxy(len(ycolumns))
    print nx, ny
    f, axs = plt.subplots(ny,nx)
    axs2d = [ item for sublist in axs for item in sublist ]
    print "axs2d:",axs2d
    for i,column in enumerate(ycolumns):
        a=axs2d[i]
        a.plot(x, data[:,column], label=str(column))
        a.legend()
    plt.show()

def main(inputfile, ycolumns, plottypes, xcolumn, **kwargs):

    # read input
    data = read_input(inputfile)

    if 1 in plottypes:
        plot1(data, ycolumns)
    if 2 in plottypes:
        plot2(data, ycolumns, xcolumn)
    if 3 in plottypes:
        plot3(data, ycolumns, xcolumn)

    return


if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser(description="Plot Utility")
    parser.add_argument("-i","--inputfile",type = str,default='PRs',help="name of the input file. default name: PRs")
    #parser.add_argument("-z","--zmatrixfile",type = str,default='ZMAT',help="name of the zmatrix file. default name: ZMAT")
    #parser.add_argument("-v","--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-Y", "--ycolumns", nargs='*', default=[1], type=int)
    parser.add_argument("-P", "--plottypes", nargs='*', default=[], type=int)
    parser.add_argument("-x", "--xcolumn", type=int, help='the column that represents the x axis if -1 no x is used')
    args=parser.parse_args()

    print args.ycolumns
    print "args:", args

    main(**vars(args))
