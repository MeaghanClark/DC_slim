#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random

# uncomment these lines when running from command line
#sys.argv = ['tree_processing.py', '../full_output/tree_pWF_full_run_42822_1.trees', './full_output', 'big_test']

seed=random.randint(1,1e6)
print(f"random seed is {seed}")

treefile = sys.argv[1]
outdir = sys.argv[2]
prefix = sys.argv[3]

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")

# read in treefile 
orig_ts = pyslim.load(treefile)

# recapitate tree
rts = pyslim.recapitate(orig_ts, recombination_rate = 1e-8, ancestral_Ne=4232)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=1.0e-8, random_seed = seed, keep=True)) 
# should increase genome size to get more mutations

print(f"The tree sequence now has {mts.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mts.diversity()}.")

# define sampling periods
before = range(0, 10000, 50)
during = range(9955, 10005, 5)
after = range(10050, 10250, 50)

sampling = [*before, *during, *after]
sample = {}
x = 0
for n in sampling: 
    sample[x] = mts.individuals_alive_at(n)
    x= x + 1

# calculate pi
for slice in range(len(sampling)): 
    groups = {
    'sample' : sample[slice],} # could add more groups to this
    group_order = ['sample'] # add in groups
    sampled_nodes = [[] for _ in groups]
    for j, k in enumerate(group_order):
        for ind in groups[k]:
            sampled_nodes[j].extend(rts.individual(ind).nodes)

    #ind divergence
    ind_nodes = []
    ind_group = []
    ind_ids = []
    for j, group in enumerate(group_order):
        for ind in groups[group]:
          ind_ids.append(ind)
          ind_nodes.append(mts.individual(ind).nodes)
          ind_group.append(group_order[j])

    nind = len(ind_ids)
    pairs = [(i, j) for i in range(nind) for j in range(nind)]
    ind_div = mts.divergence(ind_nodes, indexes=pairs) # this is what gives pi, diagonal is mean pi within population,
    # individual level divergence (heterozygosity) 

     # save output
    x = []
    for i in ind_ids:
       ind = mts.individual(i)
       label = f"slim_{ind.metadata['pedigree_id']}"
       x.append(label)

    b = np.reshape(ind_div, (nind,nind))
    panda_df = pd.DataFrame(data = b, columns = x)
    panda_df.to_csv(outdir+"/"+prefix+"_"+str(slice)+"_pi.txt", sep=' ', index=True)