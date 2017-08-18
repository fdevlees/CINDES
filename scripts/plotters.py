'''This file contains all the plotting functions I have build

   This file starts with all the utility functions the other plot functions use

   Plotters should all have:
       - a clear input format
       - no dependance on global variables
       - **kwargs support
       - a clear docstring

'''

######## IMPORT STATEMENTS



######## GLOBAL VARIABLES

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


def prop_substituent_last_cycle( totalsites, maxnsites=10, datacolumn=0)
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
