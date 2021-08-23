# HydrogenBondAnalysis
Pipeline for comparing Hydrogen Bonds with wildtype and other mutants
Hydrogen Bond Analysis script by Ebru Cetin
the bonds that occur simultaneously contributed to the total duration 
only once, while bonds occurring at different times were summed to the
total duration. In this script we track those hydrogen bonds on 
single mutants of DHFR whose occupancy deviate from that of the WT by ±30%.

# Consists of 4 scripts
1. index_extraction.tcl : gets the index number of the atoms
2. merging_bonds.py     : merges bonds according to their duration
                          writes down the h-bonds with specific threshold
3. ifunique.py          : compares the h-bonds with other mutants
4. plot_unique.py       : plots the comparison

## merging_bonds.py 
                 requires two input files for both of the Mutant and Wild Type:
                  1. name-index.dat : which can be retrived by index_extraction.tcl
                  2. name-hbonds.dat: which can be retrived by VMD > Extensions > Timeline > Hbonds > writetoafile

                gives 3 output files:
                 1.ligand h-bond info (if any, please change the ligand name to the one on your system)
                 2.merged h-bonds info
                 3.info of bonds differs +- 30 % from WT

## ifunique.py
                 Requires merged-bonds.txt file which can be produced by the previous script for all mutants in the
                 mutation matrix

## plot_unique.py
                 Requires the output of ifunique.py
