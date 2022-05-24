#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random

# use: python tree_processing.py input.trees metadata.txt prefix
# need to run in slim_3.4 environemnt (source activate slim_3.4)

# uncomment these lines when testing from command line or jupyter notebook
#sys.argv = ['tree_processing.py', './dc_slim.trees', 'metadata.txt', 'big_test', '75', '1e-8']
#args = sys.argv
treefile = sys.argv[1]
prefix = sys.argv[3]
outdir = "/mnt/home/clarkm89/DC_slim/output/"

seed=random.randint(1,1e6)
print(f"random seed is {seed}")

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")


print(f"Arguments of the script : {sys.argv[1:]=}")
print(f"Loading tree file", sys.argv[1])

# read in treefile 
orig_ts = pyslim.load(treefile)

# recapitate tree
rts = orig_ts.recapitate(recombination_rate = 1e-8, Ne=5000) #may what to change this Ne


orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Before recapitation, the max number of roots was {orig_max_roots}, "
      f"and after recapitation, it was {recap_max_roots}.")

# subset out individuals over the 1000 year time period of interest
ind_subsets = {}
nodes_subsets = {}
ind_label = {}
for x in range(121): #for each of 120 time points (100 after beg. of decline, 20 before decline)
    ind_subsets[x] = rts.individuals_alive_at(x*10)
    nSamples = len(ind_subsets[x])
    nodes_list = []
    for ind in range(nSamples):
        nodes_list.append(rts.individual(ind_subsets[x][ind]).nodes.tolist())
    nodes_subsets[x] = nodes_list
    y = []
    for i in range(nSamples):
        ind=rts.individual(i)
        label = f"{ind.metadata['pedigree_id']}" #these are tsk ids, I want slim ids
        y.append(label)
    ind_label[x] = y


# overlay mutations
ts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=1.0e-8, random_seed = seed, keep=True)) 
# should increase genome size to get more mutations

print(f"The tree sequence now has {ts.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {ts.diversity()}.")


# calculate stats for subsets
    # pwp and coal code from Gideon

pwp_dfs = {}
#coal_dfs = {}
for set in nodes_subsets:
    print(set)
    nSamples = len(ind_subsets[set]) #101
    pairs = [(i, j) for i in range(nSamples) for j in range(nSamples)]
    pwp = ts.divergence(nodes_subsets[set],indexes=pairs)
    pwp = np.reshape(pwp, (nSamples,nSamples))
    pwp_df = pd.DataFrame(data = pwp, columns = ind_label[set])
    pwp_df.to_csv(outdir+prefix+"/"+prefix+"_pwp_"+str(set*10)+".txt", sep=' ', index=True)
    print(f"done pwp")
    
    #coal = rts.divergence(nodes_subsets[set],indexes=pairs,mode="branch")
    #coal = np.reshape(coal, (nSamples,nSamples))
    #coal_df = pd.DataFrame(data = coal, columns = ind_label[set])   
    #coal_dfs[set] = coal_df
    #print(f"done coal")

print(f"finished calculating divergence statistics")



# export data for visualization
#for set in pwp_dfs: 
#    pwp_dfs[set].to_csv("mnt/home/clarkm89/DC_slim/output/"+prefix+"_pwp_"+str(set*10)+".txt", sep=' ', index=True)
#    coal_dfs[set].to_csv("mnt/home/clarkm89/DC_slim/output/"+prefix+"_coal_"+str(set*10)+".txt", sep=' ', index=True)

#print("finished exporting data")


