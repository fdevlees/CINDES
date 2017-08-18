#!/bin/env python
ttert = False # keep a distinction between tert apical and termedial. this will give double amount of coefficients
debug=True


########################
#####   IMPORTS    #####
########################
if True:
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
from CINDES4.utils.writings import log_io, sprint

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
    for i in xrange(nkinds):
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
                print item[i],
            except IndexError:
                break
        print
    return

def print_coef(coef,seq):
    '''prints all the fitted coefficients in a nice row column way'''
    if True:
        means = []
        stds = []
        for i in range(len(seq)): #max fiveteen functional groups
            print '{:7}'.format(seq[i]),
            all = []
            if args.equalsites: n=2
            else: n=10
            for j in range(n): #max 10 sites
                try:
                    print '{:10.5f}'.format(coef[j][i]),
                    all.append(coef[j][i])
                except IndexError:
                    print " "*10,
            print '{:10.3f}'.format(np.mean(all)),
            print '{:10.5f}'.format(np.std(all)),
            means.append(np.mean(all))
            stds.append(np.std(all))
            print
        print "all:",all
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
        import pandas as pd
        from matplotlib.colors import ListedColormap
        if ttert:
            labels_columns_nottert = labels[:data.shape[1]]
            labels_columns = np.concatenate( (labels_columns_nottert, labels_columns_nottert ) )
            print "labels_columns:", labels_columns
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
    print "seq:", labels
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
seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S','CCHO','CPh', 'CSOOOH']
seq = ['CH','CNHH','CNOO','CSH','COH','CCHO','CPh', 'CSOOOH']  # for substituent study of eline
seq = ['CH','N','B','O','S','P']  # for substituent study of eline
fseq= [ funcs.get(func, func) for func in seq ] #whahaha :)

###########################
####      CLASSES     #####
###########################

class Dataset(): #abstract data class
    #__metaclass__ = ABCMeta
    def __init__(self,name,*rgs,**kwargs):
        self.name = name
        return

    def extract(self,table=[],**kwargs):
        #confs, data = self.readfile(cutoff=20.0)
        if table==[]:
            print "read table from tablebin"
            column=args.column
            print("column: ", column)
            self.confs, data = self.readfile(column=column,**kwargs)
            #print "self.confs[20:22]:", self.confs[19:22]
            self.Y = self.extractY(data,column=1)
        else:
            print "read table from call"
            inds = [ item[0] for item in table ]
            self.confs = [ indtocon(item) for item in inds]
            data = [ item[1:] for item in table ]
            column = args.column-1
            self.Y = self.extractY(data,column=column)
            #self.Y = 27.2113838*self.Y

        # here symmetry operations: input: args.symmetry | self.confs | self.Y | self.syms
        if args.symmetry==True:
            datadict = dict() #dictionary to avoid duplicates
            for conf, ytje in zip(self.confs, self.Y): #loop over configurations + values
                conf_np = np.array(conf) #convert to np.array for boolean indexing
                if debug: print conf_np
                for isomer in self.syms: #for each symmetrical isomer
                    newconf = conf_np[isomer] #create a new conf. 
                    newind = contoind(newconf) #make the conf an index to be used in the dict
                    datadict[newind]=ytje #update the dict. if already present nothing will happen. 
            self.inds, self.Y = zip( *datadict.items() ) #convert the dictionary back to lists. 
            self.confs = [ indtocon(index) for index in self.inds ] #convert the indices back to configurations. 
            if debug:
                print "len confs:", len(self.confs)
                print "len Y:", len(self.Y)
                print "len dict", len(datadict)
        # output: self.confs(updated) | self.Y(updated) | self.inds | 

        self.X = self.extractX(self.confs)
        #print "extraction succesfull"
        if args.verbose<0:
            print "X:",self.X.shape
            print "Y:",self.Y.shape

        test= 0
        if test==1:
            n = min(2,len(self.Y))
            sprint(n,self.X,self.Y)
        return

    def extract2(self):
        self.X2 = self.extract2DX(self.confs)
        #PHENALENE:
        #hits = np.sum(self.X2, axis=0).reshape([11,11]).astype(int)
        # DIAMONDOIDS:
        if ttert: hits = np.sum(self.X2, axis=0).reshape([12,15*2]).astype(int)
        else: hits = np.sum(self.X2, axis=0).reshape([12,15]).astype(int)
        return hits

    def analyze(self,data):
        '''analyzes the structure of the data file'''
        print "data analyzation"
        print "verbosity:", args.verbose
        ldat = len(data)
        print "data consists of {} elements".format(ldat)
        for i in range(ldat):
            datel = data[i]
            dcolumn = [ item[args.column] for item in datel ]
            print "element {} has {} datapoints".format(i,len(datel)),
            avg = np.mean(dcolumn)
            std = np.std(dcolumn)
            print "with average of {} and std of {}".format(avg,std)
        if args.verbose>2:
            for i,item in enumerate(data[0]):
                print i, item
        return

    def openfile(self,file='tablebin'):
        ruwdata = []
        with open(file,'rb') as fid:
          if True:
            while True:
                try:
                    ruwdata.append(pickle.load(fid))
                except EOFError:
                    break
        if args.analyseinput:
            self.analyze(ruwdata)
        if args.verbose>1:
            print "first 10 of ruwdata:"
            for i in range(10):
                print ruwdata[0][i]
        return ruwdata

    def readfile(self,cutoff=False,column=1,**kwargs):
        '''extracts the data from filename: name '''
        ruwdata = self.openfile(**kwargs)
        datar=ruwdata[-1][:] # last entry in tablebin
        #########
        if False:
            data = datar[:488] + datar[1616:]
        else:
            data = datar
        if args.verbose>1:
            print "first 10 of data:"
            sprint(10,data)
        ########
        if cutoff:
            print "cutoff applied 15 eV"
            cutoff = 15.0
            data = [ [ item[0], float(item[column])] for item in datar if float(item[1])<cutoff ] 
        if False:
            data = [ [ item[0], float(item[column])] for item in datar if not 'B' in item[0] ]
        self.data = data
        confs = [ indtocon(item[0]) for item in data ]
        indices = [ item[0] for item in data ]
        self.indices = indices
        self.confs = confs
        return confs,data   

    def extractY(self,data,column=0):
        Y = np.array([ item[column] for item in data ])
        return Y

    def write_X(self):
        '''print the coefficient matrix. can only be called when
        the extractX of the child is already used''' 
        with open('Amatrix','w') as fid:
           for item in self.X:
               fid.write(str(item))
               fid.write('\n')
        return

    def difmodel(self,indices,refconf=[]):
        ''' this guesses all data only from the differences of a selected configuration '''
        #this is the procedure as used in the MC model
        # uses:
        # - self.data
        # - self.indices
        # - refconf = reference conformation. 
        #define reference configuration. normally optimum. max or min value
        param = {'reference':'arg'}
        if param['reference']== 'minimum':
            reference = min(self.data[-100:],key = lambda x:x[1])
        elif param['reference']=='maximum':
            reference = max(self.data,key = lambda x:x[1])
        elif param['reference']=='arg':
            reference = refconf
            assert not refconf==[]
        print "reference:", reference
        creference = indtocon(reference[0])
        print "creference:", creference
        #print "self.data:"
        sprint(5, self.data)
        pprint(self.data[:5])
        #Dtable = dict(self.data)
        Dtable = dict([ item[:2] for item in self.data] )
        #print "Dtable:"
        #print(Dtable)
        predictions = []
        for index in indices:
            deltaetje = 0
            rconf = indtocon(index)
            print "rconf:" , rconf, "real value:", Dtable[index]
            for i in range(len(creference)): # now we want to have a value erandom for this configuration and test it with a certain probability
                if not rconf[i] == creference[i]:
                    confje = creference[0:i] + [rconf[i]] + creference[i+1:]
                    indje= contoind(confje)
                    #print "Dtable[indje]:", Dtable[indje]
                    try:
                        deltaetje += Dtable[indje] - reference[1]
                    except KeyError:
                        print "KeyError"
                        break
                    #print "deltaetje:", deltaetje
                    print "indje:", indje, "deltaetje:", deltaetje
            erandom = float (reference[1] + deltaetje)
            #erandom = float ( deltaetje)
            print "erandom:", erandom, "reference[1]:", reference[1]
            predictions.append(erandom)
        sprint(100,predictions, self.data[-100:])
        return predictions

    def linreg_analyse2(self,clf,hits=[],model='OLS',**kwargs):
        Rscore = clf.score(self.X2,self.Y)
        print "total score {}:".format(model), Rscore 
        if args.verbose>2: print clf.coef_ #also very large coefficients

        if ttert: C = clf.coef_.reshape([12,15*2])[:]
        else: C = clf.coef_.reshape([12,15])[:]
        #C = clf.coef_.reshape([11,11])[:]
        if args.verbose>1 and args.plot:
            labels = [ funcs[func] for func in self.seq ]
            TDPlot(C,hits,labels,args)

        print "small test:"
        sprint(5,clf.predict(self.X2), self.Y)
        if args.verbose>0:
            if args.fraction:
                preds_train = clf.predict(self.X_train)
                preds_test  = clf.predict(self.X_test)
                rmse_train = metrics.mean_squared_error(preds_train, self.Y_train)
                mae_train = metrics.mean_absolute_error(preds_train, self.Y_train)
                rmse_test = metrics.mean_squared_error(preds_test, self.Y_test)
                mae_test = metrics.mean_absolute_error(preds_test, self.Y_test)
                print "small comparison prediction vs real target value (training):"
                sprint(5,preds_train, self.Y_train)
                print "small comparison prediction vs real target value (testing):"
                sprint(5,preds_test, self.Y_test)
                print "RMSE training:", rmse_train
                print "MAE training:",  mae_train
                print "RMSE test:", rmse_test
            plt.show()

        #plot the dataset vs the predictions. has to be straight line for good fit
        if args.xyplot and args.verbose>0:
            if args.fraction:
                print "training red / test bleu"
                plt.plot(preds_train,self.Y_train,'ro',alpha=0.5)
                plt.plot(preds_test,self.Y_test,'bo',alpha=0.25)
            else:
                plt.plot(clf.predict(self.X2),self.Y,'ro')
            plt.show()
        if args.fraction:
            return rmse_train, mae_train, rmse_test, rmse_train, Rscore
        else:
            return

    def linreg_analyse(self,clf,model='OLS',**kwargs):

        # SCORES
        print "total score {}:".format(model), clf.score(self.X,self.Y)
        if args.verbose>0:
            if args.fraction:
                print "training score {}:".format(model), clf.score(self.X_train,self.Y_train)
                print "testing  score {}:".format(model), clf.score(self.X_test,self.Y_test)
            try:
                print "alpha:", clf.alpha
            except AttributeError:
                pass
            try:
                print "alpha_:", clf.alpha_
                print "Cross Validation Values_:", clf.cv_values_
            except AttributeError:
                pass

        if args.verbose>2: print clf.coef_ #also very large coefficients

        if args.verbose>0:
            #prints the first 10 datapoints and the prediction.
            if args.fraction:
                preds_train = clf.predict(self.X_train)
                preds_test  = clf.predict(self.X_test)
                rmse_train = metrics.mean_squared_error(preds_train, self.Y_train)
                mae_train = metrics.mean_absolute_error(preds_train, self.Y_train)
                rmse_test = metrics.mean_squared_error(preds_test, self.Y_test)
                mae_test = metrics.mean_absolute_error(preds_test, self.Y_test)
                print "small comparison prediction vs real target value (training):"
                sprint(5,preds_train, self.Y_train)
                print "small comparison prediction vs real target value (testing):"
                sprint(5,preds_test, self.Y_test)
                print "RMSE training:", rmse_train
                print "MAE training:",  mae_train
                print "RMSE test:", rmse_test
                print "MAE test:",  mae_test
            else:
                print "small test:"
                sprint(5,clf.predict(self.X), self.Y)

        #print "outlier test"
        #print [ item for item in enumerate(self.Y) if item[1]>14.0 ]        
        #print [ (self.confs[i], self.Y[i]) for i in range(620,622) ]  

        #plot the dataset vs the predictions. has to be straight line for good fit
        if args.xyplot and args.verbose>0:
            if args.fraction:
                print "training red / test bleu"
                plt.plot(preds_train,self.Y_train,'ro',alpha=0.5)
                plt.plot(preds_test,self.Y_test,'bo',alpha=0.25)
            else:
                plt.plot(clf.predict(self.X),self.Y,'ro')
            plt.show()

        #cut the coefficient vector per functional group
        if args.verbose>1 or args.plot>0:
            C = clf.coef_
            coef = slice_it(C,self.ngps)
            #print "coef", coef
            # print them. gives back the mean and stds per functional groups
            means,stds = print_coef(coef,seq)

            # for a transposed format per site
            coeft= []
            coeft = map(lambda *row: list(row), *coef)
            coeft = [[ item if item else 0.0 for item in sublist ] for sublist in coeft ]
            self.coeft = coeft

        if args.plot>2:
            smeans=[]
            sstds=[]
            for i in range(len(self.ngps)): # for al sites do:
                total = [ item[i] for item in coeft ] #list of all for same site
                smeans.append(np.mean(total))
                sstds.append(np.std(total))
                bar_plot2(smeans,std=sstds)

            plt.show()

        if args.plot>1:
            fig,ax = bar_plot(means,seq,std=stds)
            #plt.xticks(range(len(seq)),seq)
        if args.plot>0:
            multibar_plot(coeft,seq,fig=0,ax=0)
        if args.fraction:
            return rmse_train, mae_train, rmse_test, rmse_train
        else:
            return

    def linreg(self,model,alpha=0,intercept=False,printing=0,twosite=False,**kwargs):
        '''does the linear regression and finds the useful parameters'''

        if model == 'LinearRegression':
            clf = linear_model.LinearRegression(fit_intercept=intercept)
        elif model in ['Ridge']:
            clf = linear_model.Ridge(alpha=alpha,fit_intercept=intercept,tol=0.001,solver='auto')
        elif model in ['RidgeCV','ridgecv']:
            clf = linear_model.RidgeCV(alphas=alpha, fit_intercept=intercept, store_cv_values=True)
        elif model in ['Lasso']:
            clf = linear_model.Lasso(alpha=alpha,fit_intercept=intercept,tol=0.001)
        elif model in ['ElasticNet']:
            l1_ratio = 0.1 #default 0.5
            alpha = 1e-2
            clf = linear_model.ElasticNet(alpha=alpha,l1_ratio=l1_ratio, fit_intercept=intercept,tol=0.001)

        if args.verbose>2:
            if twosite:
                sprint(5,self.X2)
            else:
                sprint(5,self.X)
            sprint(5,self.Y)
        if args.fraction:
          if twosite:
            self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.X2,self.Y, train_size = args.fraction)
            print "size training set:", np.shape(self.Y_train)
            print "size test set:", np.shape(self.Y_test)
            clf.fit(self.X_train, self.Y_train)
          else:
            self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.X,self.Y, train_size = args.fraction)
            print "size training set:", np.shape(self.Y_train)
            print "size test set:", np.shape(self.Y_test)
            clf.fit(self.X_train, self.Y_train)
        else:
            if twosite:
                clf.fit(self.X2,self.Y)
            else:
                clf.fit(self.X,self.Y)

        return clf

    def predict(self,indices,clf):
        pre_confs = [ indtocon(item) for item in indices ]
        pre_Xs = self.extractX(pre_confs)
        #pre_2Xs= self.extract2D(pre_confs)
        print "prediction. 2 Xs:"
        preds = clf.predict(pre_Xs)
        sprint(2,pre_Xs,preds)
        return preds


    def predict2(self,indices,clf):
        pre_confs = [ indtocon(item) for item in indices ]
        if debug: print "indices:", indices
        #pre_Xs = self.extractX(pre_confs)
        pre_2Xs= self.extract2DX(pre_confs)
        print "two dimensional prediction. 2 Xs:"
        preds = clf.predict(pre_2Xs)
        sprint(2,pre_2Xs,preds)
        return preds


class Adamantane(Dataset):
    '''  TYPE OF MOLECULE SUBCLASS OF DATASET   '''

    #specific adamantane parameters
    ngps = (12,12,12,12,15,15,15,15,15,15) #for every instance this is same

    def __init__(self,*rgs,**kwargs):
        self.seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S']
        self.bonds = ( (0,9),(0,5),(0,7),
                     (3,9),(3,4),(3,8),
                     (2,6),(2,8),(2,7),
                     (1,4),(1,5),(1,6) )
        if args.symmetry:
            Adasym = [
             [ 1, 2, 3, 4, 5, 6, 7, 8, 9,10],
             [ 1, 3, 4, 2, 7, 8, 9,10, 5, 6],
             [ 1, 4, 2, 3, 9,10, 5, 6, 7, 8],
             [ 3, 2, 4, 1, 6, 7, 5, 9,10, 8],
             [ 4, 2, 1, 3, 7, 5, 6,10, 8, 9],
             [ 4, 1, 3, 2, 6,10, 8, 9, 7, 5],
             [ 2, 4, 3, 1,10, 5, 9, 7, 8, 6],
             [ 3, 1, 2, 4,10, 8, 6, 7, 5, 9],
             [ 2, 3, 1, 4, 9, 7, 8, 6,10, 5],
             [ 4, 3, 2, 1, 8, 9, 7, 5, 6,10],
             [ 3, 4, 1, 2, 5, 9,10, 8, 6, 7],
             [ 2, 1, 4, 3, 8, 6,10, 5, 9, 7] ]
            self.syms = np.array(Adasym) -1
            #ADA example
            #conf = ['a','b','c','d','e','f','g','h','i','j']
            #conf_new = np.array(conf)[Asyms[3]] #for example
        return

    def extractX(self,confs):
        nC= len(confs)
        if args.verbose>2:
            for i in range(20):
                print confs[i]
            print "seq:", seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        site0 = np.zeros([nC,12])
        site1 = np.zeros([nC,12])
        site2 = np.zeros([nC,12])
        site3 = np.zeros([nC,12])
        site4 = np.zeros([nC,15])
        site5 = np.zeros([nC,15])
        site6 = np.zeros([nC,15])
        site7 = np.zeros([nC,15])
        site8 = np.zeros([nC,15])
        site9 = np.zeros([nC,15])
        if args.verbose>1:
            print "dim site0", site0.shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                if group=='CNOO60':
                    group='CNOO'
                j = seq.index(group) #find the index of the group of that sequence
                if i == 0:
                    site0[k,j]=1 #index is number of configuration , number of group
                if i == 1:
                    site1[k,j]=1
                if i == 2:
                    site2[k,j]=1
                if i == 3:
                    site3[k,j]=1
                if i == 4:
                    site4[k,j]=1
                if i == 5:
                    site5[k,j]=1
                if i == 6:
                    site6[k,j]=1
                if i == 7:
                    site7[k,j]=1
                if i == 8:
                    site8[k,j]=1
                if i == 9:
                    site9[k,j]=1
                if i == 10:
                    print 'too long'
                if i in [0,1,2,3] and j > 11:
                    print 'too far index'
                if i in [4,5] and j > 14:
                    print 'too far index'
        if args.equalsites:
            tertiair = site0 + site1 + site2 + site3
            secondair = site4 + site5 + site6 + site7 + site8 + site9
            X = np.concatenate((tertiair,secondair),axis=1)
        else:
            X = np.concatenate((site0,site1,site2,site3,site4,site5,site6,site7,site8,site9),axis=1)
        return X

    def extract2DX(self,confs):
        nC= len(confs)
        nter = 12
        nsec = 15
        X = np.zeros( [ nC , nter*nsec ] )
        #bonds = ( (0,4),(2,4),(1,5),(3,5) )
        for k in range(len(confs)): # for all the configurations:
            B = np.zeros([12,15]) #so tertiary * secondairy sites
            #B = ( B04, B24, B15, B35 )
            #loop over all combinations 
            for combi in self.bonds:
                i1,i2 = combi
                group1 = confs[k][i1]
                group2 = confs[k][i2]
                if group1=='CNOO60':
                    group1='CNOO'
                if group2=='CNOO60':
                    group2='CNOO'
                gr1 = seq.index(group1)
                gr2 = seq.index(group2)
                B[ gr1, gr2 ] += 1
            if True: # so make one total matrix were all combis are combined
                #Btotal = np.sum( B , axis = 0)
                Bflatten = B.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        return X

    def test_conf(confs):
        results = []
        for conf in confs:
            for i1,i2 in zip(*self.bonds):
                if i1=='CSH' and i2=='CF':
                    results.append(conf)
        return results

############################ END CLASS ADAMANTANE

############################ START CLASS DIAMANTANE
class Diamantane(Dataset):
    '''Diamantane class'''
    ngps = (12,12,12,12,15,15) #for every instance this is same

    def __init__(self,*rgs,**kwargs):
        self.seq = ['CH','CCHHH','CCFFF','N','CF','CCl','CNHH','CNOO','CCN','CSH','COH','CCOOH','CO','O','S','CCHO','CPh','CSOOOH']
        if args.symmetry:
            if False:
                self.syms = [ [ 0,1,2,3,4,5 ], [1,0,3,2,5,4] , 
                              [ 2,1,0,3,4,5 ], [1,2,3,0,5,4] ,
                              [ 0,3,2,1,4,5 ], [3,0,1,2,5,4] ,
                              [ 2,3,0,1,4,5 ], [3,2,1,0,5,4] ]
            else:
                self.syms = [ [ 0,1,2,3,4,5 ], [1,0,3,2,5,4] ]
        return

    def extract2DX(self,confs):
        nC= len(confs)
        nter = 12
        nsec = 15
        if ttert: X = np.zeros( [ nC, nter*nsec*2 ] )
        else: X = np.zeros( [ nC , nter*nsec ] )
        bonds = ( (0,4),(2,4),(1,5),(3,5) )  # bonds 0,4 and 1,5 are symmetrically similar and 2,4 and 3,5 idem dito
        
        for k in range(len(confs)): # for all the configurations:
            B04 = np.zeros([12,15]) #so tertiary * secondairy sites
            B24 = np.zeros([12,15]) #so tertiary * secondairy sites
            B15 = np.zeros([12,15]) #so tertiary * secondairy sites
            B35 = np.zeros([12,15]) #so tertiary * secondairy sites
            B = ( B04, B24, B15, B35 )
            #loop over all combinations 
            for combi,bmatrix in zip(bonds,B):
                i1,i2 = combi
                group1 = confs[k][i1]
                group2 = confs[k][i2]
                gr1 = seq.index(group1)
                gr2 = seq.index(group2)
                bmatrix[ gr1, gr2 ] = 1
            if not ttert: # so make one total matrix were all combis are combined
                Btotal = np.sum( B , axis = 0)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
            elif ttert: # keep a distinction between tert apical and tert medial sites. 
                Btertapical = np.sum( (B04,B15), axis=0)
                Btertmedial = np.sum( (B24,B35), axis=0)
                Btotal = np.concatenate( (Btertapical,Btertmedial),axis=1)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        print "Btotal constructed; shape X2:", np.shape(Btotal)
        return X

    def extractX(self,confs):
        nC= len(confs)
        #global seq
        #global newseq
        #newseq=[]
        #for item in seq:
        #    for conf in confs:
        #        if item in conf:
        #            newseq.append(item)
        #            break
        #print "newseq:", newseq
        #seq = newseq
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "seq:", seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        site0 = np.zeros([nC,12])
        site1 = np.zeros([nC,12])
        site2 = np.zeros([nC,12])
        site3 = np.zeros([nC,12])
        site4 = np.zeros([nC,15])
        site5 = np.zeros([nC,15])
        if args.verbose>1: print "dim site0", site0.shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                j = seq.index(confs[k][i]) #find the index of the group of that sequence
                if i == 0:
                    site0[k,j]=1 #index is number of configuration , number of group
                if i == 1:
                    site1[k,j]=1
                if i == 2:
                    site2[k,j]=1
                if i == 3:
                    site3[k,j]=1
                if i == 4:
                    site4[k,j]=1
                if i == 5:
                    site5[k,j]=1
                if i == 6:
                    print 'too long'
                if j > 14:
                    print 'too far index'
        X = np.concatenate((site0,site1,site2,site3,site4,site5),axis=1)
        if False: # try to make sites equal
            a1 = site0 + site1
            a2 = site2 + site3
            a3 = site4 + site5
            X = np.concatenate((a1,a2,a3),axis=1)
        return X

class Phenalene(Dataset):
    def __init__(self,nsites,*args,**kwargs):
        self.ngps = nsites * (12,)
        global seq
        seq = self.seq
        return

    def extractX(self,confs):
        nC= len(confs)
        nsites = len(self.ngps)
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "self.seq:", self.seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        LoS = []
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        if args.verbose>1: print "dim site0", LoS[0].shape
        for k in range(len(confs)): # for all the configurations:
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                cleangroup = split('(\d+)',group)[0]
                j = self.seq.index(cleangroup) #find the index of the group of that sequence
                LoS[i][k,j] = 1
        X = np.concatenate(LoS,axis=1)
        return X

    def extract2DX(self,confs):
        nC= len(confs)
        X = np.zeros( [ nC , 121 ] )
        bonds = ( (0,6),(5,6),(4,8),(3,8),(1,7),(2,7),
                  (0,1),(2,3),(4,5) )

        for k in range(len(confs)): # for all the configurations:
            B = []
            for combi in bonds:
                B.append( np.zeros( [11,11] ) )
            #loop over all combinations 
            for combi,bmatrix in zip(bonds,B):
                i1,i2 = combi
                group1 = confs[k][i1]
                cleangroup1 = split('(\d+)',group1)[0]
                group2 = confs[k][i2]
                cleangroup2 = split('(\d+)',group2)[0]
                gr1 = self.seq.index(cleangroup1)
                gr2 = self.seq.index(cleangroup2)
                bmatrix[ gr1, gr2 ] = 1
                bmatrix[ gr2, gr1 ] = 1
            if True: # so make one total matrix were all combis are combined
                Btotal = np.sum( B , axis = 0)
                Bflatten = Btotal.flatten()
                X[k] = Bflatten
        print "X2 constructed; shape X2:", np.shape(X)
        return X

class Thiadiazinyl(Dataset):
    def __init__(self,*args,**kwargs):
        #self.seq = [ 'N', 'CH', 'CF', 'CCFFF', 'CCHHH', 'COH', 'CNOO', 'CNHH', 'CCOOH', 'COCHHH', 'CNHCHHH' ]
        self.seq = [ 'CNHCH3', 'CSOCH3', 'COCH3', 'CSCH3', 'CSO3H', 'N', 'CCOOH', 'CCF3', 'CCH3', 'CHCO', 'CFCO', 'COOH','CSOH',
                     'CNH2'  , 'COH', 'CSH', 'CCN', 'CH', 'CF', 'CBr', 'CCl']
        global seq
        seq = self.seq
        self.ngps = 5 * (len(self.seq),)
        return

    def extractX(self,confs):
        nC= len(confs)
        nsites = len(self.ngps)
        if args.verbose>2:
            for i in range(10):
                print confs[i]
            print "self.seq:", self.seq
            print "len confs:", len(confs)
            print "len confs[0]:", len(confs[0])
        LoS = []
        for i in range(nsites):
            LoS.append( np.zeros( [nC, self.ngps[i] ] ) )
        if args.verbose>1: print "dim site0", LoS[0].shape
        for k in range(len(confs)): # for all the configurations:
            #if k<30: print confs[k]
            for i in range(len(confs[k])): #for all the groups in the configuration
                group = confs[k][i]
                #cleangroup = split('(\d+)',group)[0]
                if group=='CSO2OH':
                    cleangroup= 'CSO3H'
                elif group=='CHNCH3':
                    cleangroup= 'CNHCH3'
                else:
                    cleangroup = group
                j = self.seq.index(cleangroup) #find the index of the group of that sequence
                LoS[i][k,j] = 1
        X = np.concatenate(LoS,axis=1)
        return X



#####################################
#####       MAIN PROGRAM       ######
#####################################

def supermain(args):
    if False:
        resultlist = []
        for column in [1,2,3]:
            args.column=column
            run=main(args)
            resultlist.append(run.coeft)
        print "resultlist:"
        print resultlist
        plots_lin(resultlist)
    else:
        main(args)
    return 

def main(args):
    identify = args.identify
    if 'ada' in identify:
        myrun = Adamantane(args.file)
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane(args.file)
    elif 'phe' in identify:
        myrun = Phenalene(name=args.file,nsites=9)
    elif any(item in identify for item in ['fre','thia']):
        myrun = Thiadiazinyl(args.file)
    elif any(item in identify for item in ['fre','naph']):
        pass
        #myrun = Naphtol(args.file)
    else:
        raise SystemExit('no identify was identified')
    myrun.extract(file=args.file)
    linmodels = [ 'LinearRegression', 'Ridge']
    allerrors = []
    for i in range(args.times[0]):
        if args.ols:
            clf_LS = myrun.linreg(model=linmodels[0], intercept=args.intercept)
            if args.intercept:
                print "intercept:", clf_LS.intercept_
            if args.analyse:
                myrun.linreg_analyse(clf_LS)
        if args.ridge:
            print "args.ridge:", args.ridge
            clf_Ridge = myrun.linreg(model='Ridge', alpha=args.ridge, intercept=args.intercept)
            if args.intercept:
                print "intercept:", clf_Ridge.intercept_
            if args.analyse:
                errors = myrun.linreg_analyse(clf_Ridge,model='Ridge')
                allerrors.append(errors)
        if args.ridgecv:
            print "args.ridgeCV"
            alpha = [ 10**i for i in np.arange(-10,10,0.5) ]
            clf_RidgeCV = myrun.linreg(model='RidgeCV',alpha=alpha)
            if args.intercept:
                print "intercept:", clf_RidgeCV.intercept_
            if args.analyse:
                myrun.linreg_analyse(clf_RidgeCV,model='RidgeCV')
        if args.lasso:
            print "args.lasso"
            alpha = args.lasso
            clf_Lasso = myrun.linreg(model='Lasso',alpha=alpha)
            if args.analyse:
                myrun.linreg_analyse(clf_Lasso,model='Lasso')
        if args.twosite:
            hits = myrun.extract2()
            print hits
            if True:
                alpha=50
                print "args.ridge:", alpha
                clf_Ridge2D = myrun.linreg(model='Ridge', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_Ridge2D,hits=hits,model='Ridge')
                    allerrors.append(errors)
            if False:
                alpha=1e-2
                print "args.Lasso:", alpha
                clf_Lasso2D = myrun.linreg(model='Lasso', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_Lasso2D,hits=hits,model='Lasso')
                    allerrors.append(errors)
            if False:
                alpha=1e-2
                print "args.ElasticNet:", alpha
                clf_EN2D = myrun.linreg(model='ElasticNet', alpha=alpha,twosite=True)
                if args.analyse:
                    errors = myrun.linreg_analyse2(clf_EN2D,hits=hits,model='ElasticNet')
                    allerrors.append(errors)
    if args.ridge or args.ols or args.ridgecv or args.twosite:
        try:
            for item in allerrors:
                print " ".join(map(str,item))
            print "means RMSE_train/MAE_train/RMSE_test/MAE_test:", np.mean(allerrors,axis=0)
        except (ValueError,TypeError):
            print "error error"
            for item in allerrors:
                print item
            pass

    if False:
        new_indices = [
        'CSH_CH_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CCHHH_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CCFFF_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CNHH_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CNOO_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_COH_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CSH_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_N_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CF_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ,
        'CSH_CCl_COH_CH_CO_CNOO_CCHHH_COH_CCl_N' ]
        prediction = myrun.predict2(new_indices, clf_Ridge2D)
        print "prediction:", prediction


    #myrun.difmodel()
    return myrun


#####################################
#####     END MAIN PROGRAM     ######
#####################################


#####################################
#####   START CALL FROM CINDES   ####
#####################################
#for default args:
class defaults(object):
   def __init__(self):
       self.verbose=1
       self.plot=0
       self.twosite=True
       self.xyplot=0
       self.symmetry=False
       self.column=1
       self.analyse=False
       self.equalsites=False
       self.fraction=None
       self.intersect=False
args = defaults()


def regression(table, indices,identify,column=2, **kwargs):
    print "In call in fitter.py"
    global args
    args.column=column
    if any(item in identify for item in ['ada', 'adhoma']):
        myrun = Adamantane('ada')
    #elif 'dia' in identify:
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane('dia')
    elif 'pro' in identify:
        myrun = Propane(name='pro',nsites=3)
    else:
        raise SystemExit('No identify_ identified')
    #print "table:", table
    myrun.extract(table=table)
    #print "done extraction"

    # DETERMINE ALPHA:
    if True:
        alphas = [ 1*10**i for i in [ -4, -2, -1, 0, 1, 2, 4 ] ]
        clf_RidgeCV = myrun.linreg(model='RidgeCV', alpha=alphas)
        alpha = clf_RidgeCV.alpha_
        print "alpha used:", alpha
    else:
        alpha=1e-4
    print "args.ridge alpha parameter:", alpha

    clf = myrun.linreg(model='Ridge',alpha=alpha)
    myrun.linreg_analyse(clf,model='Ridge',**kwargs)

    predictions = myrun.predict(indices,clf)
    print "one_dimensional predictions:", predictions
    return predictions

@log_io()
def twodim_regression(table, indices,identify,column=2, **kwargs):
    print "In call in fitter.py"
    global args
    args.column=column
    if any(item in identify for item in ['ada', 'adhoma']):
        myrun = Adamantane('ada')
    elif any(item in identify for item in ['dia','dilu','diho','dimi','dima','dilumi']):
        myrun = Diamantane('dia')
    else:
        raise SystemExit('No identify_ identified')
    myrun.extract(table=table)
    hits = myrun.extract2()   #different. 
    print hits
    # DETERMINE ALPHA:
    if True:
        alphas = [ 1*10**i for i in [ -4, -2, -1, 0, 1, 2, 4 ] ]
        clf_RidgeCV = myrun.linreg(model='RidgeCV', twosite=True, alpha=alphas)
        alpha = clf_RidgeCV.alpha_
        print "alpha used:", alpha
    else:
        alpha=1e-4
    print "args.ridge alpha parameter:", alpha
    clf_Ridge2D = myrun.linreg(model='Ridge', alpha=alpha,twosite=True)
    errors = myrun.linreg_analyse2(clf_Ridge2D, model='Ridge')
    predictions = myrun.predict2(indices,clf_Ridge2D)
    print "two_dimensional predictions:", predictions
    return predictions

#####################################
#####    END CALL FROM CINDES   #####
#####################################


####################################
#####   DIFFERENCE MODEL (MC)   ####
####################################
#two steps. first get the instance and return it
def get_instance():
    myrun = Dataset('name')
    myrun.readfile()
    return myrun
def dif_predict(indices,refconf, instance):
    return instance.difmodel(indices,refconf)
##############################################
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="unpickles data stored with pickle module")
    parser.add_argument("-a","--analyse",action="store_true",help="analyze the results")
    parser.add_argument("-R","--ridgecv",action="store_true",help="do a rigde CV")
    parser.add_argument("-s","--equalsites",action="store_true",help="restrict that symmetrical similar sites have similar coefficients")
    parser.add_argument("-S","--symmetry",action="store_true",help="symmetrically similar substitutions are added")
    parser.add_argument("-c","--column",action="store",type=int,default=1,help="which column to be used. default first")
    parser.add_argument("-C","--center",action="store",type=float,default=None,help="which column to be used. default first")
    parser.add_argument("-e","--electronvolt",action="store_true",help="use electron volts instead of hartree")
    parser.add_argument("-i","--identify",type = str,default='ada',help="identify. ada or dia")
    parser.add_argument("-l","--label",type = str,default='property',help="sets the property label for plots")
    parser.add_argument("-n","--analyseinput",action="store_true",help="do a short analysis of the input")
    parser.add_argument("-N","--times",action="store",nargs=1,default=[1], type=int, help="do N times")
    parser.add_argument("-o","--ols",action="store_true",help="do an ordinary least square regression")
    parser.add_argument("-r","--ridge",action="store",nargs='?',type=float,const=1e-4,help="do a ridge regression")
    parser.add_argument("-L","--lasso",action="store",nargs='?',type=float,const=1e-4,help="do a lasso regression")
    parser.add_argument("-t","--twosite",action="store_true",help="do an ordinary least square regression")
    parser.add_argument("-p","--plot",action="count",help="make also a plot of the data")
    parser.add_argument("-f","--fraction",action="store",nargs='?',type=float,const=0.5,help="take only a fraction of data set")
    parser.add_argument("-v","--verbose", action="count", default=0, help="increase output verbosity")
    parser.add_argument("-x","--xyplot",action="store_true",help="make also an y vs predict(x) plot of the data")
    parser.add_argument("-I","--intercept",action="store_true",help="use an intercept. else no intercept is used data expected centered")
    #parser.add_argument('file', type=argparse.FileType('rb'),help="a pickled tablebin file")
    parser.add_argument('file',help="a pickled tablebin file")
    args=parser.parse_args()

    main(args)


