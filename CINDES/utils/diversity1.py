
# initially written by F. De Vleeschouwer
# 
verbose=False

import math

class Diversifier(object):
    def __init__(self, index=1, nmax=100, discardCH=False):
        self.filter=None
        self.index=index
        self.nmax=nmax
        self.discardCH=discardCH
        return

    def set_filter(self, filter=False):
        self.filter=filter
        return

    def filter(self):
        #for k in range(population_size):
        #    if (float(homo[k]) <= -7.000) and (float(lumo[k]) <= -5.000) and (float(solv[k]) <= -30.000) and (float(dipool[k]) >= 5.000):
        #        site2[0][size2]=site[0][k]
        #        site2[1][size2]=site[1][k]
        #        site2[2][size2]=site[2][k]
        #        site2[3][size2]=site[3][k]
        #        homo2[size2]=homo[k]
        #        lumo2[size2]=lumo[k]
        #        dipool2[size2]=dipool[k]
        #        solv2[size2]=solv[k]
        #        functie2[size2]=functie[k]
        #        size2 += 1
        pass

    def table_run(self, Table, index=1, seq=None):
        if self.filter:
            self.filter(Table)
        if not seq:
            seq=Table.get_seq()
        divalues, occupancy =self.diversify(
                                        X=Table.confs,
                                        Y=Table.Y,
                                        seq=seq,
                                        index=index)
        return divalues, occupancy

    def get_occupancy12(self,seq,confs, discardCH=False):

        if discardCH:
            if 'CH' in seq:
                seq.remove('CH')
                print("    NOTA BENE: CH groups are discarded in the diversity analysis")
        #n groups
        M=len(seq)
        #n sites
        N=len(confs[0])
        #n confs
        nconfs=len(confs)

        occupancy_single = [[0.000 for j in range(M)] for i in range(N)]
        occupancy_percentage_single = [[0.000 for j in range(M)] for i in range(N)]
        sum_occupancy_single = [0.000 for j in range(M)]
        sum_occupancy_percentage_single = [0.000 for j in range(M)]

        # for every site
        for i in range(N):
            # for every group
            for j in range(M):
                # for every configuration
                for conf in confs:
                    # if in the new selection of confs, the ith site of the gth conf has the jth substituent: add 1 to occupancy counter
                    #if site2[i][g] == substituent[j]:
                    try:
                        if conf[i] == seq[j]:
                            # that is the occupancy of the ith site the jth substituent:
                            occupancy_single[i][j]+=1
                    except IndexError:
                        print("conf:", conf)
                        print("seq:", seq)
                        print("i,j", i,j)
                        raise

        # get summed occupancy 1D
        # now the occupancies for each site are added up. so we get a 1D list of occupancies of every substituent
        for j in range(M):
            for i in range(N):
                sum_occupancy_single[j] += occupancy_single[i][j]

        # 1
        # in 2D occupancy. get percentages by 
        for i in range(N):
            for j in range(M):
                occupancy_percentage_single[i][j]=occupancy_single[i][j]*100./nconfs

        # 2
        # in 1D occupancy. get percentages
        for j in range(M):
            sum_occupancy_percentage_single[j]=sum_occupancy_single[j]*100./nconfs/N

        return sum_occupancy_percentage_single, occupancy_percentage_single

    def get_occupancy3(self, seq, confs, combis=None):
        # 3
        # counter
        teller=0
        M=len(seq)
        #n sites
        N=len(confs[0])
        #n confs
        nconfs=len(confs)


        # count all combi of group occurances
        combis = ( (0,9),(0,5),(0,7),
                   (1,4),(1,5),(1,6),
                   (2,6),(2,8),(2,7),
                   (3,9),(3,4),(3,8) )
        tert_positions, sec_positions = list(zip(*combis))
        n_tert, n_sec = (len(tert_positions), len(sec_positions))

        if True:
            occupancy_double = [[0.000 for j in range(M)] for i in range(M)]
            occupancy_percentage_double = [[0.000 for j in range(M)] for i in range(M)]
            for g,conf in enumerate(confs):
                # for every tert site
                for t_p in tert_positions:
                    # for every secondary site
                    for s_p in sec_positions:
                        # for every group on site i
                        for k in range(M):
                            # for every group on site j
                            for l in range(M):
                                # if site i has substituent k and site j has substituent l or vice versa
                                if (conf[t_p] == seq[k]) and (conf[s_p] == seq[l]):
                                    # counter of combination (2D list) group k and group l: add one
                                    occupancy_double[k][l]+=1
                                    # total number of combis
                                    teller+=1

        else:
            occupancy_double = [[0.000 for j in range(M)] for i in range(M)]
            occupancy_percentage_double = [[0.000 for j in range(M)] for i in range(M)]
            for g,conf in enumerate(X):
                # for every site
                for i in range(N):
                    # for every next site
                    for j in range(i+1,N):
                        # for every group on site i
                        for k in range(M):
                            # for every group on site j
                            for l in range(k,M):
                                # if site i has substituent k and site j has substituent l or vice versa
                                if (conf[i] == seq[k]) and (conf[j] == seq[l]) or (conf[i] == seq[l]) and (conf[j] == seq[k]):
                                    # counter of combination (2D list) group k and group l: add one
                                    occupancy_double[k][l]+=1
                                    # total number of combis
                                    teller+=1

        # normalize every combi occurance to percentages per combi
        for k in range(M):
            for l in range(k,M):
                occupancy_percentage_double[k][l]=occupancy_double[k][l]*100./teller
        print("n bonds in confs:", teller)

        return occupancy_double, occupancy_percentage_double

    def diversify(self, X, Y, seq=None, index=1, **kwargs):

        # Y has to be a list of prop values for each conf. not just one propvalue. so:
        if not isinstance(Y[0],list):
            Y = [[item] for item in Y]

        # sequence
        if not seq:
            seq =["CCCH", "CCFFF", "CCHHH", "CCHO", "CCl", "CCN", "CF", "CH", "CNC", "CNHH", "CNOO80", "COH", "CPh80", "CSH", "CSOOOH", "CTh80", "N"]

        # column = substituent M  and row = site N
        # n substituents
        M= len(seq)
        # nsites
        N= len(X[0]) # assuming every X has the same length. 

        # nconfs
        nconfs = len(X)

        def diversify12(index):
            occupancy_single = [[0.000 for j in range(M)] for i in range(N)]
            occupancy_percentage_single = [[0.000 for j in range(M)] for i in range(N)]
            sum_occupancy_single = [0.000 for j in range(M)]
            sum_occupancy_percentage_single = [0.000 for j in range(M)]

            # for every site
            for i in range(N):
                # for every group
                for j in range(M):
                    # for every configuration
                    for conf in X:
                        # if in the new selection of confs, the ith site of the gth conf has the jth substituent: add 1 to occupancy counter
                        #if site2[i][g] == substituent[j]:
                        if conf[i] == seq[j]:
                            # that is the occupancy of the ith site the jth substituent:
                            occupancy_single[i][j]+=1

            # get summed occupancy 1D
            # now the occupancies for each site are added up. so we get a 1D list of occupancies of every substituent
            for j in range(M):
                for i in range(N):
                    sum_occupancy_single[j] += occupancy_single[i][j]

            # 1
            # in 2D occupancy. get percentages by 
            for i in range(N):
                for j in range(M):
                    occupancy_percentage_single[i][j]=occupancy_single[i][j]*100./nconfs

            # 2
            # in 1D occupancy. get percentages
            for j in range(M):
                sum_occupancy_percentage_single[j]=sum_occupancy_single[j]*100./nconfs/N

            # make a list to store the diversity values
            # 1
            diversity_value_single = [0 for g in range(nconfs)]
            # 2
            diversity_value_singlebis = [0 for g in range(nconfs)]

            # 1 2 fill list with div values for each conf
            for g,conf in enumerate(X):
                for i in range(N):
                    for j in range(M):
                        if self.discardCH:
                            if seq[j]=='CH':
                                # div value 1
                                #diversity_value_single[g]+= 100. * nconfs/N ?
                                # div value 2
                                #diversity_value_singlebis[g]+= ?
                                continue
                        if conf[i] == seq[j]:
                            if False:
                                # div value 1
                                diversity_value_single[g]+=occupancy_percentage_single[i][j]
                                # div value 2
                                diversity_value_singlebis[g]+=sum_occupancy_percentage_single[j]
                            else:
                                # div value 1
                                diversity_value_single[g]+=occupancy_percentage_single[i][j]**2
                                # div value 2
                                diversity_value_singlebis[g]+=sum_occupancy_percentage_single[j]**2

            if verbose:
                # print div1
                for j in range(M):
                    print(seq[j], occupancy_single[0][j], occupancy_single[1][j], occupancy_single[2][j], occupancy_single[3][j])
                # print div1 percentages
                for j in range(M):
                    print(seq[j], occupancy_percentage_single[0][j], occupancy_percentage_single[1][j], occupancy_percentage_single[2][j], occupancy_percentage_single[3][j])

            # return
            if index==1:
                return diversity_value_single, occupancy_percentage_single
            elif index==2:
                return diversity_value_singlebis, sum_occupancy_percentage_single
            else:
                return

        def diversify3():
            # 3
            # counter
            teller=0

            occupancy_double = [[0.000 for j in range(M)] for i in range(M)]
            occupancy_percentage_double = [[0.000 for j in range(M)] for i in range(M)]

            # count all combi of group occurances
            for g,conf in enumerate(X):
                # for every site
                for i in range(N):
                    # for every next site
                    for j in range(i+1,N):
                        # for every group on site i
                        for k in range(M):
                            # for every group on site j
                            for l in range(k,M):
                                # if site i has substituent k and site j has substituent l or vice versa
                                if (conf[i] == seq[k]) and (conf[j] == seq[l]) or (conf[i] == seq[l]) and (conf[j] == seq[k]):
                                    # counter of combination (2D list) group k and group l: add one
                                    occupancy_double[k][l]+=1
                                    # total number of combis
                                    teller+=1

            # normalize every combi occurance to percentages per combi
            for k in range(M):
                for l in range(k,M):
                    occupancy_percentage_double[k][l]=occupancy_double[k][l]*100./teller
            print(teller)

            # list of div3 values
            diversity_value_double = [0 for g in range(nconfs)]

            # fill the list
            for g,conf in enumerate(X):
                for i in range(N):
                    for j in range(i+1,N):
                        for k in range(M):
                            for l in range(k,M):
                                if (conf[i] == seq[k]) and (conf[j] == seq[l]) or (conf[i] == seq[l]) and (conf[j] == seq[k]):
                                    diversity_value_double[g]+=occupancy_percentage_double[k][l]
            # print div3 
            for k in range(M):
                for l in range(k,M):
                    print(seq[k], seq[l], occupancy_double[k][l])

            return diversity_value_double, occupancy_double

        # 1 2
        if index==1 or index==2: divalues, occupancy = diversify12(index)
        else: divalues, occupancy = diversify3()

        # div indices 1,2,3 for each conf
        if verbose:
            print("DIVERSITY:")
            for g,conf in enumerate(X):
                #print site2[0][g], site2[1][g], site2[2][g], site2[3][g], homo2[g], lumo2[g], dipool2[g], solv2[g], functie2[g], diversity_value_single[g], diversity_value_singlebis[g], diversity_value_double[g]
                print('_'.join(conf), Y[g], divalues[g])

        return divalues, occupancy
