import math

# sequence
substituent=["CCCH", "CCFFF", "CCHHH", "CCHO", "CCl", "CCN", "CF", "CH", "CNC", "CNHH", "CNOO80", "COH", "CPh80", "CSH", "CSOOOH", "CTh80", "N"]

# column = substituent M  and row = site N

# n substituents
M=17

# nsites
N=4

# nconfs
population_size=0

size=0
size2=0
wordlist = [ ]

# n seek
max=100

# counter
teller=0

# 1
occupancy_single = [[0.000 for j in range(M)] for i in range(N)]
occupancy_percentage_single = [[0.000 for j in range(M)] for i in range(N)]
sum_occupancy_single = [0.000 for j in range(M)]
sum_occupancy_percentage_single = [0.000 for j in range(M)]

# 2
occupancy_double = [[0.000 for j in range(M)] for i in range(M)]
occupancy_percentage_double = [[0.000 for j in range(M)] for i in range(M)]

# 3
occupancy_triple = [ ]
occupancy_percentage_triple = [ ]

# read file list of words!
file=open("BFS_VH_alles")
for line in file:
    for word in line.split():
        wordlist.append(word)
        size += 1
file.close


# meaning each line in the file had 9 items site1,site2,site3,site4, homo, lumo, dipool, solv, functie
population_size = size/9
print population_size

# each item going to one of these empty boxes:

# not here that site is site is a list of sites and than every site has the substituent on that site of all confs. So it is counterintuitive! I have it
# transposed.
site = [["" for g in range(population_size)] for i in range(N)]
homo = [0 for g in range(population_size)]
lumo = [0 for g in range(population_size)]
dipool = [0 for g in range(population_size)]
solv = [0 for g in range(population_size)]
functie = [0 for g in range(population_size)]

# we are going to refine these boxes and fill these boxes:
site2 = [["" for g in range(population_size)] for i in range(N)]
homo2 = [0 for g in range(population_size)]
lumo2 = [0 for g in range(population_size)]
dipool2 = [0 for g in range(population_size)]
solv2 = [0 for g in range(population_size)]
functie2 = [0 for g in range(population_size)]

# now fill the boxes
for g in range(population_size):
    site[0][g]=wordlist[g*9]
    site[1][g]=wordlist[g*9+1]
    site[2][g]=wordlist[g*9+2]
    site[3][g]=wordlist[g*9+3]
    homo[g]=wordlist[g*9+4]
    lumo[g]=wordlist[g*9+5]
    dipool[g]=wordlist[g*9+6]
    solv[g]=wordlist[g*9+7]
    functie[g]=wordlist[g*9+8]

# we now only take one specific part of the total population matching the criteria;
# HOMO < -7
# LUMO < -5
# solv < -30
# dipool > 5
# and put the results in new boxes
for k in range(population_size):
    if (float(homo[k]) <= -7.000) and (float(lumo[k]) <= -5.000) and (float(solv[k]) <= -30.000) and (float(dipool[k]) >= 5.000):
        site2[0][size2]=site[0][k]
        site2[1][size2]=site[1][k]
        site2[2][size2]=site[2][k]
        site2[3][size2]=site[3][k]
        homo2[size2]=homo[k]
        lumo2[size2]=lumo[k]
        dipool2[size2]=dipool[k]
        solv2[size2]=solv[k]
        functie2[size2]=functie[k]
        size2 += 1

# print the number of confs fullfilling the criteria and print it
print size2
for g in range(size2):
    print site2[0][g], site2[1][g], site2[2][g], site2[3][g], homo2[g], lumo2[g], dipool2[g], solv2[g], functie2[g]


# for every site
for i in range(N):
    # for every group
    for j in range(M):
        # for every configuration
        for g in range(size2):
            # if in the new selection of confs, the ith site of the gth conf has the jth substituent: add 1 to occupancy counter
            if site2[i][g] == substituent[j]:
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
        occupancy_percentage_single[i][j]=occupancy_single[i][j]*100./size2

# 2
# in 1D occupancy. get percentages
for j in range(M):
    sum_occupancy_percentage_single[j]=sum_occupancy_single[j]*100./size2/N

# make a list to store the diversity values
# 1
diversity_value_single = [0 for g in range(size2)]
# 2
diversity_value_singlebis = [0 for g in range(size2)]

# 1 2 fill list with div values for each conf
for g in range(size2):
    for i in range(N):
        for j in range(M):
            if site2[i][g] == substituent[j]:
                # div value 1
                diversity_value_single[g]+=occupancy_percentage_single[i][j]
                # div value 2
                diversity_value_singlebis[g]+=sum_occupancy_percentage_single[j]

# 3
# count all combi of group occurances
for g in range(size2):
    # for every site
    for i in range(N):
        # for every next site
        for j in range(i+1,N):
            # for every group on site i
            for k in range(M):
                # for every group on site j
                for l in range(k,M):
                    # if site i has substituent k and site j has substituent l or vice versa
                    if (site2[i][g] == substituent[k]) and (site2[j][g] == substituent[l]) or (site2[i][g] == substituent[l]) and (site2[j][g] == substituent[k]):
                        # counter of combination (2D list) group k and group l: add one
                        occupancy_double[k][l]+=1
                        # total number of combis
                        teller+=1

# normalize every combi occurance to percentages per combi
for k in range(M):
    for l in range(k,M):
        occupancy_percentage_double[k][l]=occupancy_double[k][l]*100./teller
print teller

# list of div3 values
diversity_value_double = [0 for g in range(size2)]

# fill the list
for g in range(size2):
    for i in range(N):
        for j in range(i+1,N):
            for k in range(M):
                for l in range(k,M):
                    if (site2[i][g] == substituent[k]) and (site2[j][g] == substituent[l]) or (site2[i][g] == substituent[l]) and (site2[j][g] == substituent[k]):
                        diversity_value_double[g]+=occupancy_percentage_double[k][l]

# list of geometric sum of div1 and div3
general_diversity = [0 for g in range(size2)]

# fill list
for g in range(size2):
    general_diversity[g]=math.sqrt(diversity_value_single[g]*diversity_value_single[g]+diversity_value_double[g]*diversity_value_double[g])

# print div1
for j in range(M):
    print substituent[j], occupancy_single[0][j], occupancy_single[1][j], occupancy_single[2][j], occupancy_single[3][j] 

# print div1 percentage
for j in range(M):
    print substituent[j], occupancy_percentage_single[0][j], occupancy_percentage_single[1][j], occupancy_percentage_single[2][j], occupancy_percentage_single[3][j]

# print div3 
for k in range(M):
    for l in range(k,M):
        print substituent[k], substituent[l], occupancy_double[k][l]

# div indices 1,2,3 for each conf
print "DIVERSITY:"
for g in range(size2):
    print site2[0][g], site2[1][g], site2[2][g], site2[3][g], homo2[g], lumo2[g], dipool2[g], solv2[g], functie2[g], diversity_value_single[g], diversity_value_singlebis[g], diversity_value_double[g]

