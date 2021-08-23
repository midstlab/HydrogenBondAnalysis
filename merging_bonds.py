from collections import defaultdict
from collections import OrderedDict
import matplotlib.pyplot as plt 
import itertools
import numpy as np

file='f153s'

# merges individual bond under residue pair

def merging_bonds(file):

	# getting index info from the file
    index={}
    with open(str(file)+'-index.dat') as f:
        for line in f:
            line=line.split(' ')
            index[line[0]]=[line[1],line[2][0:-1]]


    trj=defaultdict(list)

    # opens the file name-hbonds.dat reads
    # hbonds from index information
    counter=0
    b=()
    with open(str(file)+'-hbonds.dat') as g:
        for line in g:
            a=line.split(' ')
            #print(a)
            if a[0]=='freeSelLabel':
                continue
            elif a[0][0:-1]=="#":
                continue
            elif a[0]=="#":
                continue
            elif a[0] =='freeSelString':
                b=(a[3],a[4])
                trj[b]=list()
            elif line==None:
                break
            elif line:
                trj[b].append([a[0],a[1][0:-1]])
                counter+=1

    for key,val in trj.items():
        length =len(trj[key]) #gets the length of the trajectory

    occ=[]*length #creates empty array to store the bonds
    pro={} # creates protein dictionary
    
    #gets the index numbers as key and writes h-bond info throughout the trajectory
    for key,val in trj.items():
    	for line in val:
    		occ.append(float(line[1]))
    	pro[key]=occ
    	occ=[]
    
    #creates and ordered dictionary
    final=OrderedDict(sorted(pro.items(),key=lambda x:x[0]))
    
    res=defaultdict(list)
    
    #gets the index names as TRP139, SER123 and atom names 'O','HD22' and matches them to index numbers 
    for key,val in final.items():
        if key[0] in index.keys():
            new_key_don = index.get(key[0])
            new_key_don = tuple(new_key_don)
            new_key_acc = index.get(key[1][0:-1])
            new_key_acc=tuple(new_key_acc)
            new_key=new_key_don+new_key_acc
            val1=final.get((key[0],key[1]))
            res[new_key].append(val1)
    
    
    # merges all h-bond trajectories into an array of arrays
    # e.g. SER3-VAL93 : [[][][]]
    res_pair = defaultdict(list)
    for key,val in res.items():
        #print(key)
        res_pair[key[0],key[2]].append(val)
    

    occ_dict={}

    for x,y in res_pair.items():
    	list1=[]
    	for a in y:
    		list1.append(a[0])
    	occ_dict[x]=list1
    
    empty={}
    for x in occ_dict:
    	empty[x]=[0]*length #for each dictionary item in occ_dict it creates an array of '0's



    for x in empty:
        for a in occ_dict:
            if x==a: # if the dictionary keys are same
                for m1 in range(len(occ_dict[a])):
                    list1=occ_dict[a][m1] # list1 reads the h-bond info at a time instant
                    for m2 in range(len(list1)):
                        if empty[x][m2]==0 and list1[m2]==1: # merging of h-bonds durations 
                            empty[x][m2]=1
                        
                                 
    final_dict={}
    
    # counts the discrete event and normalize it to 100.
    for i in empty:
        count=0
        for j in range(len(empty[i])):
            if empty[i][j] ==1:
               count+=1
        pre_sum=count/length*100
        final_dict[i]=round(pre_sum,2)
        
    print('Hydrogen bonds found in system', str(file),'is',len(final_dict)) # gives the total of merged bonds
    return(final_dict)

#merges the h-bond data of WT trajectories
def wt_occ(wt1,wt2):
    wtc={}
    for key in wt1:
        for ley in wt2:
            val1=wt1.get(key)
            val2=wt2.get(ley)
            if key == ley:
                #print(val1)
                average=(val1+val2)/2
                wtc[key]=round(average,2)
            elif key in wt1 and not wt2:
                average=(abs(val1)/2)
                wtc[key]=round(average,2)
            elif key in wt2 and not wt1:
                average=(abs(val2)/2)
                wtc[key]=round(average,2)
    return wtc

# calculates the difference between mutant and WT trajectories
def occ_diff(mutant,wtc): 
    difference={}

    for hbond1 in mutant:
        for hbond2 in wtc:
            a=[]
            if hbond1==hbond2:
                val1=mutant.get(hbond1)
                val2=wtc.get(hbond2)
                difference[hbond1]=[float(val1),float(val2),round(float(val1)-float(val2),2)]
    add_list={}

    #add wtc hbonds not present in mutant
    for d in wtc:
        if d not in mutant:
            val1=wtc.get(d)
            z=round(val1,2)
            add_list[d]=[0,z,-z]

    #add mutant hbonds not present in wtc
    for s in mutant:
        if s not in wtc:
            val3=mutant.get(s)
            y=round(val3,2)
            add_list[s]=[y,0,y]

    difference.update(add_list)

    top_bonds={}

    #gets the h-bonds with +-30 % deviation
    for key in difference:
        for item in [difference[key][len(difference[key])-1]]:
            if abs(item) >= 30 :
                val=difference.get(key)
                if key in top_bonds:
                    top_bonds.append(key)
                else:
                    top_bonds[key]=val
    
    #collects the h-bond data of the ligand
    dhf={}
    for key,value in difference.items():
        if 'DHF'in key[0] or 'DHF' in key[1]:
            #print('bugn Salimin doum gunu')
            dhf[key]=value
    
    m_file=str(file)+"-dhf-bonds.txt"

    with open(m_file,'w') as f:
        for key, value in dhf.items():
            f.write('%s-%s :  %s %s %s\n' % (key[0],key[1], value[0],value[1],value[2]))

    return top_bonds

# writes merged bonds into a file

def merged_bonds_write(file):
    m_file=str(file)+"-merged-bonds.txt"

    with open(m_file,'w') as f:
        for key, value in mutant.items():
            f.write('%s-%s :  %s\n' % (key[0],key[1], value))

    w_file='wt'+'-merged-bonds'

    with open('{}.txt'.format(str(w_file)),'w') as k:
        for bond, percent in wt.items():
            k.write('%s-%s :  %s\n' % (bond[0],bond[1], percent))

    diff_file=str(file)+'-wtc-30%-bonds.txt'

    with open(diff_file,'w') as g:
        for key,val in wt_compared.items():
            g.write('%s-%s:    %s %s %s\n' % (key[0],key[1], val[0],val[1],val[-1]))

#merging bonds for all systems
mutant=merging_bonds(file)
wt1=merging_bonds('wt1')
wt2=merging_bonds('wt2')

#calculation of WT occupancy
wt=wt_occ(wt1,wt2)
# calculation of mutant difference
wt_compared=occ_diff(mutant,wt)
#writes merged bonds into a file.
merged_bonds_write(file)