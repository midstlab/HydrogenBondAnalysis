import matplotlib.pyplot as plt 

name = 'f153s' # input('Enter the filename:')

fullname=str(name)+"-unique-bonds"


with open('{}.txt'.format(str(fullname)) , "r") as f:
    
    
    for line in f:
        

        line = line.split(':')
        key=line[0].strip()
        occ=line[1].split(' ')
        value=occ[4:]
        value[-1]=value[-1][0:-1]
        last_list=[0]*len(value)	
        
        for k in range(len(value)):
            last_list[k]=round(float(value[k]),2)	
       
        #print(len(last_list))
        #print(last_list)
        
        names=["WT","I5F","M20I","P21L","A26T","D27E","L28R","W30G","W30R","I94L","R98P","F153S"]
        #print(len(names))
        plt.figure()
       	plt.bar(names,last_list)
        locs, labels = plt.xticks()
        plt.xticks(locs, names,  ha='right', rotation=90, fontsize=8)
        plt.ylabel('Occupancies (%)')
        plt.xlabel('Mutants')
        plt.title(str(key))
        plt.tight_layout(pad=2)
        fullname=name.upper()+'-'+str(key)+'.png'
        plt.savefig(fullname,dpi=100)
        last_list=[]