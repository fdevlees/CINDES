import numpy as np
from copy import deepcopy

from CINDES.algorithms.pyevolve import G1DList , GSimpleGA, GAllele, Mutators, Initializators, Selectors, Consts, DBAdapters, Crossovers
import CINDES.algorithms.pyevolve as pyevolve

DE=False

# main function:
def reduce_conflicts(molecule, core=[], active=[], passive=[]):
    print "\n{0} {1} {0}".format("$"*20, molecule.index)
    conf = molecule.conf
    if DE:
        print "initial conf:", conf
        print "initial index", molecule.index

    new_conf = []
    for group in conf:
        if is_float(group[-1]):
            new_conf.append(group)
        elif len(group)>2 or group in [['C','Ph'], ['C','Th']]:
            group.append(1.0)
            new_conf.append(group)
        else:
            pass
            # no relevant dihedral for:", group
    print "new_conf:", new_conf
    molecule.dihedrals = [item[-1] for item in new_conf]
    print "dihedrals:", molecule.dihedrals

    import pickle
    TZMat={'core':core, 'passive':passive, 'active':active}
    #with open('TZMat','w') as f: pickle.dump(TZMat, f)
    try:
        run_pyevolve(molecule, **TZMat)
    except ValueError:
        print "\t Error in run_pyevolve to reduce conflicts"
    if DE:
        print "final conf:", molecule.conf
        print "final index", molecule.index
    return

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def slice_it(li, splits, ngps=None):
    '''splits the large coefficients vector in small vectors per funct. group'''
    start = 0
    lis = []
    nkinds = len(splits)
    for i in xrange(nkinds):
        if splits[i]==0: continue
        stop = start+splits[i]
        lis.append(li[start:stop])
        start = stop
    return lis

#def get_xyz(molecule, core,active,passive):
def get_xyz(molecule, core, active, passive):
    #print "core,active,passive:", core, active, passive
    from CINDES.evaluation.construction import constructor2
    import pickle
    #with open('TZMat','r') as f: TZMAT=pickle.load(f)
    ncore = len(core)

    #TZMATGLOBAL = deepcopy(TZMat)

    #from pprint import pprint
    #print "in get_xyz; TZMATGLOBAL:"
    #pprint(TZMATGLOBAL)
    c = deepcopy(core)
    a = deepcopy(active)
    p = deepcopy(passive)
    #c = deepcopy(TZMAT['core'])
    #a = deepcopy(TZMAT['active'])
    #p = deepcopy(TZMAT['passive'])
    #molecule.set_zmat( zcon.constructor2(molecule.conf,c,a,p, links=myrun.symlinks) )
    try:
        molecule.set_zmat( constructor2(molecule.conf,c,a,p) )
    except IndexError:
        print "indexerror:", molecule.conf
        raise
    molecule.zmatoxyz()
    #print molecule.xyz
    xyz = np.array([ item[1] for item in molecule.xyz ])
    #print "xyz:", xyz
    return xyz[ncore:]

def get_distance_matrix(xyz):
    from scipy.spatial.distance import pdist
    dist = pdist(xyz, metric='euclidean')
    #print "dist:", dist
    return dist

class Fitness_Function():
    '''A class for the fitness functions. The class contains the attr's needed for evaluation. Here this will be
    the database that does not change. maybe even the kernel. see which part stays here and what part has to be.
    done in learning.py.
    Aim is that i can use:
        evaluator = FitnessFunction()
    and subsequently in each iteration
        evaluator.evaluate(population)
    '''
    def __init__(self, molecule, **kwargs):
        '''for evaluation i need at least to have the database and the core / active / passive (all in zmatrix)
        i probably should also already get a self.kernel here such that the evaluatefunction only should call predict
        '''
        self.molecule = molecule
        self.kwargs = kwargs
        self.dihedralindices = []
        for i, group in enumerate(self.molecule.conf):
            if len(group)>2 or group in [['C','Ph'], ['C','Th']]:
                self.dihedralindices.append(i)
            else:
                #print "no relevant dihedral for:", group
                pass
        return

    def eval_function(self, dihedrals):
        dihedrals=map(int,dihedrals)
        #return sum(self.dihedrals_to_dist2(dihedrals))
        return self.dihedrals_to_minvalue(dihedrals)

    def dihedrals_to_dist(self, dihedrals):
        self.dihedrals_to_conf(dihedrals)
        # set xyz for that conf:
        #TZMat = deepcopy(self.TZMat.copy())
        #xyz = get_xyz(self.molecule, **TZMat)
        xyz = get_xyz(self.molecule, **self.kwargs)
        #del TZMat
        # now focus only on xyz of not the core.
        #xyz = xyz[self.ncore:]
        dist= get_distance_matrix(xyz)
        #print "len dist:", len(dist)
        dist=dist[dist<2.0]
        #print "len reduced dist:", len(dist)
        # and focus only on distances < 2.0 angstrom
        #print "dist:", dist
        #raise SystemExit('stop')
        return dist

    def dihedrals_to_minvalue(self, dihedrals):
        ''' a more clever evaluation function hopefully '''
        self.dihedrals_to_conf(dihedrals)
        # set xyz for that conf:
        #TZMat = deepcopy(self.TZMat.copy())
        #global TZMat
        #import pprint
        #print "in dihedrals_to_minvalue: TZMAT:", pprint.pformat(TZMat)
        xyz = get_xyz(self.molecule, **self.kwargs)
        #xyz = get_xyz(self.molecule)
        # now focus only on xyz of not the core.
        #print "len xyz:", len(xyz)
        #xyz = xyz[self.ncore:]

        # maybe juse myrun.adj? NO
        #print "self.conf:", self.molecule.conf
        #print "len xyz:", len(xyz)

        # get list of natoms site
        natoms_site = []
        for group in self.molecule.conf:
            if is_float(group[-1]): group=group[:-1]
            if group == [ 'C', 'Ph']: natoms_site.append(11)
            else:
                natoms_site.append(len(group)-1)
        #print "natoms_site:", natoms_site
        XYZs = slice_it(xyz, natoms_site)
        #print XYZs
        minvalues=[]
        for i, xyz1 in enumerate(XYZs):
            for j, xyz2 in enumerate(XYZs):
                if i>=j: continue
                # get distance matrix between xyz1 and xyz2
                from scipy.spatial.distance import cdist
                dist12 = cdist(xyz1, xyz2)
                # take the minimum value of that distance matrix
                mind = np.min(dist12)
                # append that minimum value to a list
                minvalues.append(mind)
        #print "minvalues:", minvalues
        #print len(minvalues)
        #raise SystemExit('raise')
        # return sum of that list
        return np.min(minvalues)

    def dihedrals_to_conf(self, dihedrals):
        # set conf related to the given dihedrals
        for index, dihedral in zip(self.dihedralindices, dihedrals):
            if is_float(self.molecule.conf[index][-1]):
                self.molecule.conf[index][-1]=dihedral
            else:
                self.molecule.conf[index].append(dihedral)
        #print "self.molecule.conf:", self.molecule.conf
        return self.molecule.conf

    def get_final_molecule(self, dihedrals):
        dihedrals=map(int,dihedrals)
        self.dihedrals_to_conf(dihedrals)
        print "molecule.conf:", self.molecule.conf
        dihedralindex = '_'.join([''.join(map(str,item)) for item in self.molecule.conf])
        print "dihedralindex", self.molecule.index
        # index has to be without dihedrals? index has to stay the same for the whole event
        from string import digits
        #self.molecule.index= self.molecule.index.translate(None, digits)
        print "self.molecule.index:", self.molecule.index
        get_xyz(self.molecule, **self.kwargs)
        #get_xyz(self.molecule)
        return

def GoodEnoughDistance(ga_engine):
    '''This function determines the convergence. After 50 generations we release the criterium
    for a faster convergence'''
    if ga_engine.currentGeneration > 250:
        distance_goal = 0.5
    elif ga_engine.currentGeneration > 150:
        distance_goal = 0.6
    elif ga_engine.currentGeneration >100:
        distance_goal = 0.9
    elif ga_engine.currentGeneration > 50:
        distance_goal = 1.2
    else:
        distance_goal = 1.5
    return ga_engine.bestIndividual().score > distance_goal

#def run_pyevolve(molecule, TZMat):
def run_pyevolve(molecule, **kwargs):
    '''options should be a Run instance having at least:
        options.nsites
        options.

    '''
    # 0.
    #print "TZMat:", TZMat
    function = Fitness_Function(molecule, **kwargs)

    # 1. Enable the logging system:
    pyevolve.logEnable()

    # 2. Set Genome instance using as allelles the sites with the different functionalisations.
    genome = G1DList.G1DList(len(molecule.dihedrals))
    genome.setParams(rangemin=1, rangemax=360, gauss_mu=0.0, gauss_sigma=30.0)

    # 3. Set evaluator function (objective function) or set precalculation is True! this circumvents serial evaluation
    genome.evaluator.set(function.eval_function)
    # 4. Set mutator function
    genome.mutator.set(Mutators.G1DListMutatorRealGaussian)
    #genome.mutator.set(Mutators.G1DListMutatorIntegerRange)

    # 5. Set initalizator function
    genome.initializator.set(Initializators.G1DListInitializatorReal)
    # 6. set Crossover function: types: G1DListCrossoverUniform, G1DListCrossoverSinglePoint, G1DListCrossoverTwoPoint
    genome.crossover.set( Crossovers.G1DListCrossoverUniform)
    #print "genome:\n", genome

    # 7. set Genetic Algorithm Instance using a defined random.seed()
    ga = GSimpleGA.GSimpleGA(genome)
    # 8. set Selector
    #ga.selector.set(Selectors.GRouletteWheel)
    ga.selector.set(Selectors.GTournamentSelector)
    pop = ga.getPopulation()
    pop.setParams(tournamentPool=3)
    # 9. set NGEN (number of generations)
    ga.setGenerations(500)
    # 10. set min / max (optimize to a maximum or to a minimum)
    ga.setMinimax(Consts.minimaxType["maximize"])
    # 11. set MUP (mutation probability)
    ga.setMutationRate(0.2) #i added this from another example
    # 12. set CXP (crossover probability)
    ga.setCrossoverRate(0.5)
    # 13. set termination at convergence?:
    #ga.terminationCriteria.set(GSimpleGA.ConvergenceCriteria)
    ga.terminationCriteria.set(GoodEnoughDistance)
    # 14. set population size
    ga.setPopulationSize(50)
    # 15. set elitism
    ga.setElitism(True)
    ga.nElitismReplacement = 1

    # 17. for plotting / logging
    sqlite_adapter = DBAdapters.DBSQLite(dbname='dihedralstats.db', identify=molecule.index, resetDB=False, resetIdentify=True, commit_freq=10)
    ga.setDBAdapter(sqlite_adapter)

    #print "GenAlg:", ga

    # Do the evolution, with stats dump
    ga.evolve(freq_stats=5)
    best = ga.bestIndividual()
    print "best individual:", best.score
    #from pprint import pprint
    #pprint(vars(best))
    print "genomeList:", best.genomeList
    function.get_final_molecule(best.genomeList)
    return


if __name__=="__main__":
    TZMat = {'active': [[['C', '6', '1.54903973', '5', '108.91333933', '2', '-179.98390819', '0'],
             ['H', '55', '1.09176330', '6', '112.06998651', '5', '-171.13704403', '0'],
             ['H', '55', '1.09167289', '6', '112.03579178', '5', '-51.15297205', '0'],
             ['H', '55', '1.09175360', '6', '112.01194966', '5', '68.83644586', '0']],
            [['C', '10', '1.54921201', '9', '108.94362986', '3', '179.95687730', '0'],
             ['H', '11', '1.09176866', '10', '112.05073775', '9', '-51.47725968', '0'],
             ['H', '11', '1.09174395', '10', '112.09504358', '9', '68.50461826', '0'],
             ['H', '11', '1.09178666', '10', '112.05054335', '9', '-171.47444129', '0']],
            [['C', '2', '1.54923438', '1', '108.93464720', '3', '-179.96945401', '0'],
             ['H', '23', '1.09176233', '2', '112.06064717', '1', '-51.61060335', '0'],
             ['H', '23', '1.09182317', '2', '112.04098011', '1', '68.39930015', '0'],
             ['H', '23', '1.09173390', '2', '112.04062392', '1', '-171.60135636', '0']],
            [['C', '3', '1.54925065', '1', '108.92146008', '10', '150.00531303', '0'],
             ['H', '67', '1.09182245', '3', '112.07689138', '1', '-171.28945217', '0'],
             ['H', '67', '1.09180604', '3', '112.06271316', '1', '-51.30893450', '0'],
             ['H', '67', '1.09177572', '3', '112.02838756', '1', '68.70013607', '0']],
            [['C', '9', '1.57632685', '3', '117.95148184', '1', '175.31105103', '0'],
             ['H', '71', '1.09432277', '9', '107.90565685', '3', '-158.69820353', '0'],
             ['H', '71', '1.09277204', '9', '110.96721855', '3', '-43.33576984', '0'],
             ['H', '71', '1.07772118', '9', '117.58478077', '3', '82.09742691', '0']],
            [['C', '7', '1.57630940', '6', '118.01550169', '5', '-63.32673757', '0'],
             ['H', '47', '1.09429991', '7', '107.89775061', '6', '-158.95097787', '0'],
             ['H', '47', '1.09272470', '7', '111.01193082', '6', '-43.61710874', '0'],
             ['H', '47', '1.07758475', '7', '117.55750858', '6', '81.86218998', '0']],
            [['C', '8', '1.57650954', '2', '108.76813744', '1', '-169.95870308', '0'],
             ['H', '19', '1.09432845', '8', '107.93700896', '2', '77.44018991', '0'],
             ['H', '19', '1.09267759', '8', '110.96045425', '2', '-167.21638404', '0'],
             ['H', '19', '1.07761425', '8', '117.53455849', '2', '-41.74304368', '0']],
            [['C', '5', '1.57628544', '2', '108.74135116', '1', '68.70601974', '0'],
             ['H', '27', '1.09432481', '5', '107.93485193', '2', '77.35483350', '0'],
             ['H', '27', '1.09269374', '5', '111.00297384', '2', '-167.28620767', '0'],
             ['H', '27', '1.07751605', '5', '117.49554865', '2', '-41.79063633', '0']],
            [['C', '1', '1.57637137', '3', '108.73256302', '9', '-170.00509577', '0'],
             ['H', '43', '1.09438599', '1', '107.94587694', '3', '77.31004739', '0'],
             ['H', '43', '1.09263001', '1', '110.94708547', '3', '-167.33774953', '0'],
             ['H', '43', '1.07772202', '1', '117.53248236', '3', '-41.90845000', '0']],
            [['C', '4', '1.57619379', '3', '108.71980063', '1', '-169.97617626', '0'],
             ['H', '59', '1.09439984', '4', '107.87748469', '3', '77.22376614', '0'],
             ['H', '59', '1.09277593', '4', '110.99835097', '3', '-167.39663710', '0'],
             ['H', '59', '1.07746434', '4', '117.54423667', '3', '-41.93550658', '0']]],
     'core': [['C'],
          ['C', '1', '1.62538422'],
          ['C', '1', '1.62518435', '2', '108.41725190'],
          ['C', '3', '1.62533764', '1', '109.97947022', '2', '-60.67402437', '0'],
          ['C', '2', '1.62505501', '1', '109.99495138', '3', '60.68391879', '0'],
          ['C', '5', '1.62543914', '2', '108.39090520', '1', '-60.66167680', '0'],
          ['C', '6', '1.62488302', '5', '109.98989640', '2', '-60.69409440', '0'],
          ['C', '2', '1.62506885', '1', '109.99892884', '3', '-60.62571678', '0'],
          ['C', '3', '1.62511114', '1', '110.00476594', '2', '60.60759658', '0'],
          ['C', '9', '1.62468104', '3', '108.44083009', '1', '-60.63318114', '0']],
      'passive': [[['C', '9', '1.57639961', '3', '108.70708366', '1', '68.73609527', '0'],
              ['H', '15', '1.09433292', '9', '107.92685895', '3', '77.32537907', '0'],
              ['H', '15', '1.09267526', '9', '110.98727735', '3', '-167.32828385', '0'],
              ['H', '15', '1.07772229', '9', '117.56808334', '3', '-41.87386066', '0']],
             [['C', '5', '1.57658004', '2', '117.92456583', '1', '175.27553001', '0'],
              ['H', '31', '1.09435401', '5', '107.92998868', '2', '-158.80768173', '0'],
              ['H', '31', '1.09263553', '5', '110.97909064', '2', '-43.46580175', '0'],
              ['H', '31', '1.07759149', '5', '117.57724681', '2', '81.96277040', '0']],
             [['C', '8', '1.57661690', '2', '117.92670608', '1', '-63.39566224', '0'],
              ['H', '35', '1.09438449', '8', '107.91534511', '2', '-158.77567614', '0'],
              ['H', '35', '1.09265608', '8', '110.95877434', '2', '-43.40832916', '0'],
              ['H', '35', '1.07762570', '8', '117.57442589', '2', '82.02411699', '0']],
             [['C', '1', '1.57649273', '3', '117.92818651', '9', '-63.44953939', '0'],
              ['H', '39', '1.09439609', '1', '107.97552554', '3', '-158.59694490', '0'],
              ['H', '39', '1.09265403', '1', '110.94489950', '3', '-43.25643296', '0'],
              ['H', '39', '1.07770490', '1', '117.52444170', '3', '82.15281411', '0']],
             [['C', '7', '1.57655075', '6', '108.74191446', '5', '-169.93805677', '0'],
              ['H', '51', '1.09441815', '7', '107.95528472', '6', '77.36719139', '0'],
              ['H', '51', '1.09266567', '7', '110.96550016', '6', '-167.27023022', '0'],
              ['H', '51', '1.07768601', '7', '117.55768167', '6', '-41.86231877', '0']],
             [['C', '4', '1.57627078', '3', '117.98756267', '1', '-63.37988600', '0'],
              ['H', '63', '1.09435666', '4', '107.86906878', '3', '-158.87631015', '0'],
              ['H', '63', '1.09267571', '4', '111.01425701', '3', '-43.52129584', '0'],
              ['H', '63', '1.07753033', '4', '117.56731587', '3', '81.97274476', '0']]]}

    conf = [['C', 'Ph', 1.0],
             ['C', 'Ph'],
              ['C', 'S', 'H', 60.],
               ['C', 'H'],
                ['C', 'S', 'H', 135.],
                 ['O'],
                  ['C', 'C', 'H', 'O'],
                   ['C', 'N', 'O', 'O'],
                    ['C', 'N', 'H', 'H'],
                     ['C', 'N', 'H', 'H']]
    conf = [['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph'], ['C', 'Ph']]
    conf = [['C', 'Ph', 69], ['C', 'Ph', 168], ['C', 'Ph', 352], ['C', 'Ph', 188], ['C', 'Ph', 253], ['C', 'Ph', 20], ['C', 'Ph', 343], ['C', 'Ph', 238], ['C',
        'Ph', 289], ['C', 'Ph', 250]]
    from molecule import Molecule
    mol = Molecule(conf=conf)
    if True:
        #[ reduce_conflicts(mol, **TZMat) for _ in xrange(10) ]
        reduce_conflicts(mol, **TZMat)

