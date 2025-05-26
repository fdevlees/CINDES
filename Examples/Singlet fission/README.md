# Singlet fission project

In this project, we explored the combinatorial chemical compound space of **fulvene** derivatives using the **best-first search** algorithm to identify optimal substitution patterns satisfying the primary thermodynamic singlet fission condition (i.e., **E(S1)∼ 2E(T1)**).

In the given example, 6 sites were considered, symmetrically paired, with a substituent library of size 15 containing NO2, CN, CF3, Cl, F, BF2, BH2, H, CH3, SH, OCH3, OH, NH2, SiH3, N(CH3)2.

## Related publication

Irene Casademont-Reig, Roger Monreal-Corona, Eline Desmedt, Freija De Vleeschouwer, and Mercedes Alonso. Pursuit of Singlet Fission Fulvenes Candidates using Inverse Design. Submitted to Digital Discovery.

## Input file

The **INPUT** file is explained below:

program gaussian    <span style="color:purple"> _computational chemistry software_ </span>
procedure bfs
property function
 tdenergy_S[0] - 2 * tdenergy_T[0] 
optimum maximum
hasimagfreq
startcalcs
  1.1
    nprocs 16
    walltimelimit 48h
    geom ZMAT_fulvene
    xyz_fulvene g_fulvene hasimagfreq
    0 1 # opt=calcfc freq M062X/6-311+G(d,p)
  2.1
    nprocs 16
    walltimelimit 8h
    geom _fulvene
    etenergies_S tdenergy_S 
    0 1 # M062X/6-311+G(d,p) td=(singlets,nstates=3)
  2.1
    nprocs 16
    walltimelimit 8h
    geom _fulvene
    etenergies_T tdenergy_T
    0 1 # M062X/6-311+G(d,p) td=(triplets,nstates=3)
endcalcs
timestep 300
zmatrixfile ZMAT_fulvene
extrawaittime 10
debug
#startind CNMeMe_CSiHHH_CNHH
#sequences 1
#1 0 2
try_ready
readtable
#SYMLINKS every nr in left column needs to be lower than every nr in right column; may need to change sites order to achieve this
nprocs 8
sites 7 11 15 27 23 19
symlinks 3
1 4
2 5
3 6
identify ID_fulvene
natomscore 6
nch3 6
END
nsubsit1 15
C N O O
C C N
C C F F F
C Cl
C F
C B F F
C B H H
C H
C C H H H
C S H
C O C H H H
C O H
C N H H
C Si H H H
C N Me Me
nsubsit2 15
C N O O
C C N
C C F F F
C Cl
C F
C B F F
C B H H
C H
C Si H H H
C C H H H
C S H
C O C H H H
C O H
C N H H
C N Me Me
nsubsit3 15
C N O O
C C N
C C F F F
C Cl
C F
C B F F
C B H H
C H
C Si H H H
C C H H H
C S H
C O C H H H
C O H
C N H H
C N Me Me

