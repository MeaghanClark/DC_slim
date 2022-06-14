#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random
np.set_printoptions(threshold=sys.maxsize)


# uncomment these lines when running from command line
#sys.argv = ['tree_processing.py', '../slim_output_05312022/tree_nWF_100_1.trees', '/Users/meaghan/Desktop/DC_slim/het', 'relatedness_test', 1e-8, 5]
# [0] -- python script name
# [1] -- tree file
# [2] -- outdir
# [3] -- prefix
# [4] -- mu
# [5] -- gen time

seed=random.randint(1,1e6)
print(f"random seed is {seed}")

treefile = sys.argv[1]
outdir = sys.argv[2]
prefix = sys.argv[3]
mu = sys.argv[4]
gen_time = sys.argv[5]

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")
print(f"mu is {mu}")
print(f"gen time is {gen_time}")

# read in treefile 
orig_ts = pyslim.load(treefile)

print(f"Loaded tree file")

# recapitate tree
rts = pyslim.recapitate(orig_ts, recombination_rate = 1e-8, ancestral_Ne=5487)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
rate = float(mu)/float(gen_time)
print(f"using rate {rate}")

mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=rate, random_seed = seed, keep=True)) 
# should increase genome size to get more mutations or set mutation rate to 2.59e-5

print(f"The tree sequence now has {mts.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mts.diversity()}.")

# define sampling periods
before = range(0, 450, 50)
during = range(450, 500, 5)
after = range(500, 750, 50)

sampling = [*before, *during, *after]
print(sampling)

# caculate heterozygosity and relatedness 
pedigree_id = []
het = []
gen = [] 
rel = []

for n in sampling: 
    ind_nodes = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        ind_nodes.append(ind.nodes)
        
    # make vector of per-individual heterozygosities:
    ind_het = mts.diversity(ind_nodes, mode="site")
   
    # make vector of pedigree_ids 
    x = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        label = f"slim_{ind.metadata['pedigree_id']}"
        x.append(label)
        
    # save het output for nth sampling point
    pedigree_id = np.append(pedigree_id, x) 
    het = np.append(het, ind_het)
    gen = np.append(gen, np.repeat(n, len(ind_het))) # generation

    # make matrix of genetic relatedness for this sampling point
    nind = len(mts.individuals_alive_at(n))
    pairs = [(i, j) for i in range(nind) for j in range(nind)]
    ind_rel = mts.genetic_relatedness(ind_nodes, indexes=pairs)
    #ind_rel = np.append(n, ind_rel)
    rel = np.append(rel, str(ind_rel)) # convert relatedness into a string so there aren't issue with different numbers of individuals
        

# Output HET data for all sampling points
# assemble into dictionary for het data
het_data = { 'pedigree_id':pedigree_id, 
            'het':het, 
            'gen':gen,
           }

# convert het to dataframe and export
het_df = pd.DataFrame(data=het_data)
het_df.to_csv(outdir+"/"+prefix+"_het.txt", sep=' ', index=True)


# export relatedness 
rel_df = pd.DataFrame(data=rel)
rel_df.to_csv(outdir+"/"+prefix+"_relatedness.txt", sep=',', index=True)
#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random

# uncomment these lines when running from command line
#sys.argv = ['tree_processing.py', '../slim_output_05312022/tree_nWF_100_1.trees', '/Users/meaghan/Desktop/DC_slim/het', 'relatedness_test', 1e-8, 5]
# [0] -- python script name
# [1] -- tree file
# [2] -- outdir
# [3] -- prefix
# [4] -- mu
# [5] -- gen time

seed=random.randint(1,1e6)
print(f"random seed is {seed}")

treefile = sys.argv[1]
outdir = sys.argv[2]
prefix = sys.argv[3]
mu = sys.argv[4]
gen_time = sys.argv[5]

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")
print(f"mu is {mu}")
print(f"gen time is {gen_time}")

# read in treefile 
orig_ts = pyslim.load(treefile)

print(f"Loaded tree file")

# recapitate tree
rts = pyslim.recapitate(orig_ts, recombination_rate = 1e-8, ancestral_Ne=5487)

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
rate = float(mu)/float(gen_time)
print(f"using rate {rate}")

mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=rate, random_seed = seed, keep=True)) 
# should increase genome size to get more mutations or set mutation rate to 2.59e-5

print(f"The tree sequence now has {mts.num_mutations} mutations, "
      f"and mean pairwise nucleotide diversity is {mts.diversity()}.")

# define sampling periods
before = range(0, 450, 50)
during = range(450, 500, 5)
after = range(500, 750, 50)

sampling = [*before, *during, *after]
print(sampling)

# caculate heterozygosity and relatedness 
pedigree_id = []
het = []
gen = [] 
rel = []

sampling = [0]
for n in sampling: 
    ind_nodes = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        ind_nodes.append(ind.nodes)
        
    # make vector of per-individual heterozygosities:
    ind_het = mts.diversity(ind_nodes, mode="site")
   
    # make vector of pedigree_ids 
    x = []
    for i in mts.individuals_alive_at(n):
        ind = mts.individual(i)
        label = f"slim_{ind.metadata['pedigree_id']}"
        x.append(label)
        
    # save het output for nth sampling point
    pedigree_id = np.append(pedigree_id, x) 
    het = np.append(het, ind_het)
    gen = np.append(gen, np.repeat(n, len(ind_het))) # generation

    # make matrix of genetic relatedness for this sampling point
    nind = len(mts.individuals_alive_at(n))
    pairs = [(i, j) for i in range(nind) for j in range(nind)]
    ind_rel = mts.genetic_relatedness(ind_nodes, indexes=pairs)
    #ind_rel = np.append(n, ind_rel)
    rel = np.append(rel, str(ind_rel)) # convert relatedness into a string so there aren't issue with different numbers of individuals
        

# Output HET data for all sampling points
# assemble into dictionary for het data
het_data = { 'pedigree_id':pedigree_id, 
            'het':het, 
            'gen':gen,
           }

# convert het to dataframe and export
het_df = pd.DataFrame(data=het_data)
het_df.to_csv(outdir+"/"+prefix+"_het.txt", sep=' ', index=True)


# export relatedness 
rel_df = pd.DataFrame(data=rel)
rel_df.to_csv(outdir+"/"+prefix+"_relatedness.txt", sep=',', index=True)
