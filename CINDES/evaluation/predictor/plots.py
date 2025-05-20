#!/bin/env python
ttert = False # keep a distinction between tert apical and termedial. this will give double amount of coefficients
debug=True


########################
#####   IMPORTS    #####
########################
if False:
    import seaborn as sns
    sns.set(style="white")
    #pass two degree values of hues. 
    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    #cmap = sns.diverging_palette(20, 220, as_cmap=True)

import pickle
from pprint import pprint
import matplotlib.pyplot as plt
import matplotlib
from re import findall, split
import numpy as np
from sklearn import linear_model, metrics
from sklearn.model_selection import train_test_split
from abc import ABCMeta, abstractmethod
import pandas as pd
from CINDES.utils.writings import log_io, sprint

from itertools import islice

def take(n, iterable):
    "Return first n items of the iterable as a list"
    return list(islice(iterable, n))

#################################
#####    BASIC FUNCTIONS    #####
#################################

def cycle(iterable=('orangered','tomato','red','indianred','darkred','deeppink')):
    # cycle('ABCD') --> A B C D A B C D A B C D ...
    saved = []
    for element in iterable:
        yield element
        saved.append(element)
    while saved:
        for element in saved:
              yield element

def residuals(p,y,x):
    err = y - sum([ item * x for item in p ] )
    return err

def peval(x,p):
    return sum( [ item*x for item in p ] )

def slice_it(li, splits):
    '''splits the large coefficients vector in small vectors per funct. group'''
    start = 0
    lis = []
    if args.equalsites:
        nkinds = 2
        splits = (12,15)
    else:
        nkinds = len(splits)
    for i in range(nkinds):
        stop = start+splits[i]
        lis.append(li[start:stop])
        start = stop
    return lis

def sprint_old(n,*args,**kwargs):
    '''tries to prints the first n items of iterable objects'''
    #if kwargs:
    #    dictlist = kwargs.items()
    #    args = tuple(dictlist) + args
    for i in range(n):
        for item in args:
            try:
                print(item[i], end=' ')
            except IndexError:
                break
        print()
    return

def print_coef(coef,seq):
    '''prints all the fitted coefficients in a nice row column way'''
    if True:
        means = []
        stds = []
        for i in range(len(seq)): #max fiveteen functional groups
            print('{:7}'.format(seq[i]), end=' ')
            all = []
            if args.equalsites: n=2
            else: n=10
            for j in range(n): #max 10 sites
                try:
                    print('{:10.5f}'.format(coef[j][i]), end=' ')
                    all.append(coef[j][i])
                except IndexError:
                    print(" "*10, end=' ')
            print('{:10.3f}'.format(np.mean(all)), end=' ')
            print('{:10.5f}'.format(np.std(all)), end=' ')
            means.append(np.mean(all))
            stds.append(np.std(all))
            print()
        print("all:",all)
    return means,stds

def get_color(index,colors=['r', 'g', 'b', 'y']):
    i = index%len(colors)
    return colors[i]


def multibar_plot(X,seq,std=0,fig=0,ax=0):
    '''makes a bar plot for each site each coefficient'''
    N = len(seq)
    offset = 0.15
    ind = offset + np.arange(N)
    width=0.7
    if fig==0:
        fig1 = plt.figure()
        ax = fig1.add_subplot(111)
    from cycler import cycler
    colors = sns.color_palette("deep",n_colors=10)
    #colors = ('orangered','darkcyan','red','indianred','darkred','deeppink','g','b','y','c')
    #labels = ('site1','site2','site3','site4','site5','site6')
    labels = tuple( 'site '+str(i+1) for i in range(len(X[0])) )
    N = len(X)
    #ax.set_prop_cycle('color',cycle(['b','r','g','c','k','y','m']))
    print("I'M HERE!")
    #ax.set_prop_cycle(cycler('color',['b','r','g','c','k','y','m']))
    rects = len(seq) * [None]
    for i in range(len(seq)):
        try:
            nsit = len(X[i])
            inds = i + offset + (width/nsit) * np.arange(nsit)
            print("inds:", inds)
            rects[i] = ax.bar(inds, X[i], width/nsit,color=colors,label=labels)
            #df = DataFrame(inds,X[i],width/nsit,columns=labels)
            #df.plot(type='bar')
            #if i==0:
            #    ax.legend()
        except IndexError:
            pass
    legend = ax.legend(rects,labels,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    for i in range(len(X[0])):
       legend.legendHandles[i].set_color(colors[i])
    ax.yaxis.grid(True)
    #handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles[::-1], labels[::-1])
    property = 'Ionization Potential'
    label = 'average contribution to the ' + property + ' (eV)'

    ax.set_ylabel(label)
    ax.set_xlabel('functional group')
    plt.xticks(ind+width/2,seq)
    #plt.title('average property contribution vs. functional group')
    plt.xlim(0,len(seq))
    plt.show()
    return

def bar_plot(X,seq,std=0):
    '''makes a simple bar plot for the mean and sdv of each func. group'''
    N = len(seq)
    offset = 0.15
    ind = offset + np.arange(N)
    width=0.7
    fig,ax = plt.subplots(1)
    if std==0:
        rects1 = ax.bar(ind, X, width, color='lightpink')
    else:
        rects1 = ax.bar(ind, X, width, color='lightpink',yerr=std)
    
    ax.set_ylabel('average property contribution')
    ax.set_xlabel('functional group')
    plt.xticks(ind+width/2,seq)
    plt.title('average property contribution vs. functional group')
    plt.xlim(0,len(seq))

    if False:
        def autolabel(rects):
            for rect in rects:
                height = rect.get_height()
                ax.text(rect.get_x() + rect.get_width()/2., 
                        1.05*height,
                        '%d' % int(height),
                        ha='center', va = 'bottom')
        autolabel(rects1)
    return fig,ax

def bar_plot2(X,std=0):
    '''simple plot of mean and SD for each site'''
    N = len(X)
    offset = 0.15
    Ns = np.arange(1,N+1)
    ind = offset + Ns
    width=0.7
    fig,ax = plt.subplots(1)
    if std==0:
        rects1 = ax.bar(ind, X, width, color='red')
    else:
        rects1 = ax.bar(ind, X, width, color='red',yerr=std)
    
    ax.set_ylabel('average property contribution')
    ax.set_xlabel('site')
    plt.xticks(ind+width/2,Ns)
    plt.title('average property contribution per site')
    #plt.xlim(0,15)
    return fig,ax

def TDPlot(data,hits,labels,args):
    '''two dimensional data plot. data should be a numpy.array.

    the plotted matrix: data
    the annotations matrix: hits #this can be the same as data
    the labels for the heatmap axis: labels
    the colorbar label: args.label
    the colormap used: cmap #global variable
    '''
    import matplotlib.ticker as ticker
    data[data == 0.00000] = np.nan
    print("matrix:", data)
    #pprint(map(list,list(data)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if False:
        cax = sns.heatmap(data, cmap=cmap,
            square=True, xticklabels=labels, yticklabels=labels,
            linewidths=.5, fmt="d", cbar_kws={"shrink": .5}, ax=ax )
        cax2= sns.heatmap(hits,annot=True,alpha=0.0,fmt="d",cbar=False)
        ax.set_xticklabels(labels[0:16],rotation='vertical')
        ax.set_yticklabels(labels[:13][::-1],rotation='horizontal')
    elif False:
        cax = ax.matshow(data, interpolation='nearest')
        fig.colorbar(cax)
    else:
        #convert the data to a nice pandas.DataFrame. Seaborn likes dataframes.
        import pandas as pd
        from matplotlib.colors import ListedColormap
        if ttert:
            labels_columns_nottert = labels[:data.shape[1]]
            labels_columns = np.concatenate( (labels_columns_nottert, labels_columns_nottert ) )
            print("labels_columns:", labels_columns)
            data_pd = pd.DataFrame(data=data,index=labels[:data.shape[0]],columns=labels_columns)
        else: data_pd = pd.DataFrame(data=data,index=labels[:data.shape[0]],columns=labels[:data.shape[1]])
        cax2=sns.heatmap(hits,annot=True,alpha=0.0,fmt="d",cbar=False,annot_kws={"color":'k'})
        cax= sns.heatmap(data_pd,
                         cmap=cmap,
                         #cmap = ListedColormap( cmap.colors[::-1] ),
                         square=True,
                         linewidths=1.0,
                         center=args.center,
                         cbar_kws={ "label" : "contribution to " + args.label }
                         )
        plt.xticks(rotation=45)
        plt.yticks(rotation=45)
        plt.xlabel('secondary positions')
        plt.ylabel('tertiary positions')
    print("seq:", labels)
    if False and not args.fraction:
        for (i, j), z in np.ndenumerate(hits):
            if not z==0:
                ax.text(j, i, '{:4d}'.format(z), ha='center', va='center')
    plt.show()
    return

def plots_lin(slist):
    def plot1(data,i,vmin=-1.2,vmax=0.8):
        if False:
            center = (vmin+vmax) / 2.0
        else: center = 0.0
        bla = (311,312,313)
        titles= ('HOMO-LUMO gap','HOMO','LUMO')
        ax = fig.add_subplot(bla[i])
        ax.set_title(titles[i])
        labels =  ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S']
        flabels=  [ funcs[label] for label in labels]
        if ttert:columns=['tert. apical','tert. medial','secondary']
        else:   columns=['tertiary'    ,'secondary']
        data_pd = pd.DataFrame(data=data,index=flabels[:data.shape[0]],columns=columns).transpose()
        print("data_pd:", data_pd)
        im = sns.heatmap(data_pd,annot=True,square=True,cbar=i==0,
                 ax=ax,
                 cmap=cmap,
                 center=center,
                 linewidths=.5,
                 fmt="5.2f",
                 cbar_ax= None if i else cbar_ax,
                 cbar_kws=None if i else { "label" : "contribution to property (eV)" },
                 vmin =  vmin,
                 vmax =  vmax)
        plt.yticks(rotation=0,fontsize=12)
        plt.xticks(rotation=45,fontsize=12)
        return im

    fig = plt.figure()
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #determine min and max
    A = np.asarray(slist)
    vmin = min(A.flatten())
    vmax = max(A.flatten())

    for i,item in enumerate(slist):
        npdata = np.asarray(item)
        npdata[npdata == 0.00000] = np.nan
        if ttert:#for diamantane.
            npdata = npdata[:,[0,2,4]]
        print("npdata:", npdata)
        im = plot1(npdata,i,vmin=vmin,vmax=vmax)
    fig.tight_layout( rect=[ 0, 0, .8,1])
    plt.show()

def indtocon(index):
    #return [list(item) for item in index.split('_')]
    #return  [ findall('[A-Z][^A-Z]*',item) for item in index.split('_') ]
    return index.split('_')

def contoind(conf):
    return '_'.join([''.join(item) for item in conf])

#++++++++++++++++++++++++++
#++++     GLOBALS    ++++++
#++++++++++++++++++++++++++

eV=27.2113838

matplotlib.rcParams['mathtext.default']='regular'
funcs = {'CCFFF': '$C-CF_3$',
         'CCHHH': '$C-CH_3$',
         'CCN': '$C-C\\equiv N$',
         'CCl': '$C-Cl$',
         'CF': '$C-F$',
         'CH': '$C-H$',
         'CNHH': '$C-NH_2$',
         'CNOO': '$C-NO_2$',
         'CCOOH': '$C-COOH$',
         'CO': '$C=O$',
         'COH': '$C-OH$',
         'CSH': '$C-SH$',
         'N': '$N$',
         'O': '$O$',
         'S': '$S$'}
seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S']
fseq= [ funcs[func] for func in seq ] #whahaha :)



