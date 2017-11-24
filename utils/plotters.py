#!/bin/env python
'''This file contains all the plotting functions I have build

   This file starts with all the utility functions the other plot functions use

   Plotters should all have:
       - a clear input format
       - no dependance on global variables
       - **kwargs support
       - a clear docstring

'''

######## IMPORT STATEMENTS

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

######## GLOBAL VARIABLES

ttert = False
#cmap = sns.diverging_palette(220, 20, as_cmap=True)
funcs_r = {'CCFFF': '$C-CF_3$',
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
         'S': '$S$',
         'CNHCHHH' : '$CNHCH_3$',
         'COCHHH'  : '$COCH_3$' }
class Funcs(object):
    def __init__(self):
        self.funcs = funcs_r
    def __getitem__(self, key):
        try:
            ret = self.funcs[key]
        except KeyError:
            ret = '${}$'.format(key)
        return ret
funcs = Funcs()

######## UTILITY FUNCTIONS


####### PLOTTERS

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

def bar_plot2(X,std=0, xlabels=None):
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
    ax.set_xticklabels(xlabels)
    plt.xticks(ind+width/2,Ns)
    plt.title('average property contribution per site')
    #plt.xlim(0,15)
    return fig,ax

def heatmap_2d(data,hits,labels,args,cmap='viridis'):
    '''two dimensional data plot. data should be a numpy.array.

    the plotted matrix: data
    the annotations matrix: hits #this can be the same as data
    the labels for the heatmap axis: labels
    the colorbar label: args.label
    the colormap used: cmap #global variable
    '''
    import matplotlib.ticker as ticker
    data[data == 0.00000] = np.nan
    print "matrix:", data
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
        print "BLA"
        import pandas as pd
        from matplotlib.colors import ListedColormap
        if ttert:
            labels_columns_nottert = labels[:data.shape[1]]
            labels_columns = np.concatenate( (labels_columns_nottert, labels_columns_nottert ) )
            print "labels_columns:", labels_columns
            data_pd = pd.DataFrame(data=data,index=labels[:data.shape[0]],columns=labels_columns)
        else: data_pd = pd.DataFrame(data=data,index=labels[:data.shape[0]],columns=labels[:data.shape[1]])
        #cax2=sns.heatmap(hits,annot=True,alpha=0.0,fmt="d",cbar=False,annot_kws={"color":'k'})
        cax= sns.heatmap(data_pd,
                         cmap=cmap,
                         #cmap = ListedColormap( cmap.colors[::-1] ),
                         square=True,
                         linewidths=1.0,
                         center=None,
                         cbar_kws={ "label" : "contribution to " + args.label }
                         )
        plt.xticks(rotation=45)
        plt.yticks(rotation=45)
        plt.xlabel('secondary positions')
        plt.ylabel('tertiary positions')
    print "seq:", labels
    if False and not args.fraction:
        for (i, j), z in np.ndenumerate(hits):
            if not z==0:
                ax.text(j, i, '{:4d}'.format(z), ha='center', va='center')
    plt.show()
    return

def heatmap_1d(raw_data, ylabels=None, xlabels=None):
    '''1d regression plot'''
    debug=True

    print raw_data

    # make sure data is a numpy array
    data = np.asarray(raw_data)
    if ylabels is None:
        nsites = data.shape[1]
        ylabels = [ "site {}".format(i+1) for i in range(nsites) ]
    print "xlabels", xlabels
    print "ylabels", ylabels
    data_pd = pd.DataFrame(data=data.T, index=ylabels, columns=xlabels).astype(float)
    print "data_pd:"
    print data_pd

    # set zeros to NAN.
    data[data==0.0]=np.nan

    # DOES NOT WORK: data_pd = data_pd.sort_values('secondary',axis='columns',ascending=False)

    # start figure:

    if False:
        fig=plt.figure()
        ax=fig.add_subplot(111)
        cax=sns.heatmap( data_pd,
                         cmap='viridis',
                         square=True,
                         linewidths=1.0,
                         #center=args.center
                         )
    else:
        fig, ax= plt.subplots(1,1)
        cbar_ax = fig.add_axes([.905, .3, .05, .3])
        sns.heatmap(data_pd,
                square=True,
                ax=ax,
                cmap='viridis',
                center=None,
                fmt='.2g',
                cbar_ax=cbar_ax,
                annot=True)
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.ylabel('sites')
    plt.xlabel('functional groups')
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
        print "data_pd:", data_pd
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
        print "npdata:", npdata
        im = plot1(npdata,i,vmin=vmin,vmax=vmax)
    fig.tight_layout( rect=[ 0, 0, .8,1])
    plt.show()

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
    print "I'M HERE!"
    #ax.set_prop_cycle(cycler('color',['b','r','g','c','k','y','m']))
    rects = len(seq) * [None]
    for i in range(len(seq)):
        try:
            nsit = len(X[i])
            inds = i + offset + (width/nsit) * np.arange(nsit)
            print "inds:", inds
            rects[i] = ax.bar(inds, X[i], width/nsit,color=colors,label=labels)
            #df = DataFrame(inds,X[i],width/nsit,columns=labels)
            #df.plot(type='bar')
            #if i==0:
            #    ax.legend()
        except IndexError:
            pass
    legend = ax.legend(rects,labels,bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
    for i in range(len(X[0])):
        try:
            legend.legendHandles[i].set_color(colors[i])
        except IndexError:
            print "i:", i
            print "X[0]:", X[0]
            print "colors:", colors
            pass
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

def prop_substituent_last_cycle( totalsites, maxnsites=10, datacolumn=0):
    ''' This plots the dependence of the property vs changing the substituent
    for each site.

    Input:
        - totalsites
        - maxnsites (the maximum number of sites)
        - datacolumn
    '''
    import seaborn as sb
    import re
    tags_r=['ro','bs','g^','c*','mp','y|','k+','rd','bv','gh']
    tags_r=['-ro','-bs','-g^','-c*','-mp','-y|','-k+','-rd','-bv','-gh']
    tags_r=['-o','-s','-^','-*','-p','-<','->','-d','-v','-h']
    tags = Cycle(tags_r)

    colors = sb.hls_palette(nsites+1,l=.4) #l=lightness the smaller the darker. 

    i = len(totalsites)-1
    #do for each site:
    for j in range(len(totalsites[i])):
        run=totalsites[i][j]
        x = np.array(range(len(run)))
        if len(run) == maxnsites:
            my_xticks = [ funcs[ re.split('[0-9]',item[0])[0] ] for item in run ]
            plt.xticks(x,my_xticks)
            plt.xticks(rotation=45)
            print "my_xticks", my_xticks
        y = [ item[1][datacolumn] for item in run ]
        itje = j
        plt.plot(x,y,tags[itje],label=' site:' + str(j+1), color=colors[itje])

    if True:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_dpi(100)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        legend=ax.legend(loc='center left', fancybox=True, framealpha=0.5, bbox_to_anchor=(1,.5),fontsize=12)

    #plt.ylabel('ionization potential (a.u.)')
    #plt.ylabel('ionization potential (eV)')
    plt.ylabel(args.label)
    plt.xlabel('substituent')
    plt.title('property vs substituent')
    return
