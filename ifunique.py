## 

from itertools import *
import numpy as np
from collections import defaultdict

#creation of empty dicts
mutant={}
base={}

target_name= 'f153s' 


full_target=str(target_name)+"-wtc-30%-bonds" #full name of the required file

#this matrix allows to find mutated residues in different mutants
mutant_matrix ={ 
                 "a26t" :['ALA26', 'THR26','None','None'], # None is reversed for a second mutation
                 "d27e" :['ASP27', 'GLU27','None','None'],
                 "f153s" :['PHE153', 'SER153','None','None'],
                 "i5f" :['ILE5', 'PHE5','None','None'],
                 "i94l" :['ILE94', 'LEU94','None','None'],
                 "l28r" :['LEU28', 'ARG28','None','None'],
                 "m20i" :['MET20', 'ILE20','None','None'],
                 "p21l" :['PRO21', 'LEU21','None','None'],
                 "r98p" :['ARG98', 'PRO98','None','None'],
                 "w30g" :['TRP30', 'GLY30','None','None'],
                 "w30r" :['TRP30', 'ARG30','None','None'],
                 "n59a" :['ASN59', 'ALA59','None','None']
}

with open ('{}.txt'.format(str(full_target)) , "r") as g:
    for line2 in g:
    	#stripping the values from the text
        line2 = line2.split(':')
        bond=line2[0].strip()
        bond_couple=bond.split('-')
        donor=bond_couple[0]
        acceptor=bond_couple[1]
        value2=line2[-1].strip()
        occ=value2.split(' ')

        m_occ=occ[0] # mutant occupancy
        wt_occ=occ[1] # wt occupancy
        diff=occ[2] # difference in occupancy

        # gets the mutant h-bond information
        # reads the wt residue name for that mutation
        if target_name in mutant_matrix:
            gy=mutant_matrix.get(target_name)
            if gy[1] in donor: # if SER153 in donor then donor name for other mutants is PHE153
               donor=gy[0]
            elif gy[1] in acceptor:
               acceptor=gy[0]
            elif gy[3] in donor: # enters if there is a second mutation
                donor=gy[2]
            elif gy[3] in acceptor:
                acceptor=gy[2]
    
        bond_key=donor+'-'+acceptor
        
        if bond_key in base:
           base[bond_key].append((float(m_occ)))
        else:
           base[bond_key]=[float(m_occ)]
        

final = defaultdict(list)
database=['wt','i5f','m20i','p21l','a26t','d27e','l28r','w30g','w30r','i94l','r98p','f153s']


for key in database:
        full_name= key+"-merged-bonds" #name in the database 
    
        with open('{}.txt'.format(str(full_name)) , "r") as f:
            for line in f:
                line = line.split(':')
                ley=line[0].strip()
                b_couple=ley.split('-')
                don=b_couple[0]
                acc=b_couple[1]
                value=line[-1].strip()

                if key == 'wt':
                    mutant[ley]=float(value)
                    continue
                else:
                    if don or acc in mutant_matrix.values(): # same string operation as in compared mutant
                        zy=mutant_matrix.get(key)
                        if zy[1] in don:
                           don=zy[0]
                        elif zy[1] in acc:
                           acc=zy[0]
                        elif zy[3] in don:
                           don=zy[2]
                        elif zy[3] in acc:
                           acc=zy[2]

                b_key=don+'-'+acc
                if b_key in mutant:
                    mutant[b_key] += float(value)
                else:
                    mutant[b_key]=float(value[0:-1])
       
        for z in base:
            if z not in mutant:
                final[z].append(0) #if the bond is not found in the %30 different bonds append "0" for that mutant
            else:
                value3=mutant.get(z)
                final[z].append(value3) # otherwise add the value to the list
        mutant={}


unique_file=str(target_name)+"-unique-bonds.txt"

with open(unique_file,'w') as f:
    for key, value in final.items():
        #f.write('%s:    %s %s %s' % (key,value[0],value[1],value[2]))
        f.write('%s:    %s %s %s %s %s %s %s %s %s %s %s %s\n' % (key, value[0],value[1],value[2],value[3],value[4],value[5],value[6],value[7],value[8],value[9],value[10],value[11]))

        