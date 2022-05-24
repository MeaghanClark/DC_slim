#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random

# uncomment these lines when running from command line
#sys.argv = ['tree_processing.py', '../full_output/tree_nWF_full_run_42922_75.trees', '/Users/meaghan/Desktop/DC_slim/het', 'big_test']

seed=random.randint(1,1e6)
print(f"random seed is {seed}")

treefile = sys.argv[1]
outdir = sys.argv[2]
prefix = sys.argv[3]

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")

# read in treefile 
orig_ts = pyslim.load(treefile)

print(f"Loaded tree file")

# recapitate tree
rts = pyslim.recapitate(orig_ts, recombination_rate = 1e-8, ancestral_Ne=4232)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=2.59e-5, random_seed = seed, keep=True)) 
# should increase genome size to get more mutations or set mutation rate to 2.59e-5

print(f"The tree sequence now has {mts.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mts.diversity()}.")

# define sampling periods
before = range(0, 10000, 50)
during = range(9955, 10005, 5)
after = range(10050, 10250, 50)

sampling = [*before, *during, *after]    

pedigree_id = []
het = []
gen = [] 

for n in sampling: 
    ind_nodes = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        ind_nodes.append(ind.nodes)
    # the vector of per-individual heterozygosities:
    ind_het = mts.diversity(ind_nodes, mode="site")
   # mean_het = np.mean(ind_het)

    # save output
    x = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        label = f"slim_{ind.metadata['pedigree_id']}"
        x.append(label)
    pedigree_id = np.append(pedigree_id, x) 
    het = np.append(het, ind_het)
    gen = np.append(gen, np.repeat(n, len(ind_het))) 


het_data = { 'pedigree_id':pedigree_id, 
            'het':het, 
            'gen':gen,
           }

panda_df = pd.DataFrame(data=het_data)
panda_df.to_csv(outdir+"/"+prefix+"_pi.txt", sep=' ', index=True)
