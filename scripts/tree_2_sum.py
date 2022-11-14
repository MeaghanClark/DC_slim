#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random
import datetime # FOR DEBUGGING 
import csv
np.set_printoptions(threshold=sys.maxsize)

# uncomment these lines when running from command line
# sys.argv = ['tree_processing.py', '../slim_output_11082022/tree_nWF_5_10_89.trees','../slim_output_11082022/metaInd_nWF_5_10_89.txt', '/Users/meaghan/Desktop/DC_slim/het', '11_8_test', 1e-8, 5, 5]
# arguments: 
# [0] -- python script name
# [1] -- tree file
# [2] -- meta file
# [3] -- outdir
# [4] -- prefix
# [5] -- mu
# [6] -- gen time
# [7] -- avg age

# print some record keeping output
seed=random.randint(1,1e6)
print(f"random seed is {seed}")

treefile = sys.argv[1]
print(f"treefile is {treefile}")

metafile = sys.argv[2]

outdir = sys.argv[3]
print(f"outdir is {outdir}")

prefix = sys.argv[4]
print(f"prefix is {prefix}")

mu = sys.argv[5]
print(f"mu is {mu}")

gen_time = sys.argv[6]
print(f"gen time is {gen_time}")

avg_age = sys.argv[7]
print(f"average age is {avg_age}")

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

# load metadata

metadata = pd.read_csv(metafile, sep = "\t")

# define sampling periods
before = range(0, 450, 50)
during = range(450, 500, 5)
after = range(500, 750, 50)

sampling = [*before, *during, *after]
print(sampling)

# make cycle key
# needs to change for different lifespans
if avg_age == 2: 
    burn = 126300
elif avg_age == 5: 
    burn = 315900  # changed b/c of my slim error with % 50 
elif avg_age == 10: 
    burn = 635600
elif avg_age == 20: 
    burn = 1264600

before = range(burn-200, burn, 50)
during = range(burn, burn+50, 5)
after = range(burn+50, burn+550, 50)

cycles = [*before, *during, *after][::-1]

# make conversion df
convert_time = pd.DataFrame({'tskit_time':sampling, 'slim_time':cycles}, columns=['tskit_time', 'slim_time'])

# initalize data arrays
## individual
pedigree_id = []
het = []
gen = [] 
age = []

## pairwise
rel = []
pi = []

# loop through time points to calculate het, pi, and relatedness using tskit
for n in [*range(0, 24, 1)]: 

    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")

    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = mts.individuals_alive_at(tskit_time)

    # define pedigree ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    x = [mts.individual(i).metadata['pedigree_id'] for i in alive]   
        
    # make list of nodes for individuals sampled in slim
    ind_nodes = []
    for i in samp_pdids.to_numpy():                          # for each individual sampled in slim
        focal_ind = mts.individual(x.index(i))               # get inidvidual id by matching pedigree id to tskit id
        ind_nodes.append(focal_ind.nodes)                    # make list of nodes
    print(f"length of ind_nodes is {len(ind_nodes)}")
    
        # make vector of per-individual heterozygosities:
    ind_het = mts.diversity(ind_nodes, mode="site")
    het = np.append(het, ind_het)

    print(f"length of ind_het is {len(ind_het)}")
    
    # make array of pedigree_ids 
    ped_ids = samp_pdids.to_numpy() 
    # save het output for nth sampling point
    pedigree_id = np.append(pedigree_id, ped_ids) 
    
    # make array of ages
    ind_ages = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["age"])
    age = np.append(age, ind_ages)
    
    # make array of what cycle it is 
    gen = np.append(gen, np.repeat(tskit_time, len(ind_het))) # generation
    
    print(f"done calculating het for sampling point {n}")
    print('Timestamp: {:%H:%M:%S}'.format(datetime.datetime.now()))
        
    # define pairs
    nind = len(samp_pdids)
    pairs = [(i, j) for i in range(nind) for j in range(nind)]
        
    print(f"done defining pairs at sampling time {n}")
    print('Timestamp: {:%H:%M:%S}'.format(datetime.datetime.now()))
    print(f"there are {len(pairs)} pairs")
        
    # make matrix of genetic relatedness for this sampling point
    ind_rel = mts.divergence(ind_nodes, indexes=pairs)
    rel = np.append(rel, str(ind_rel)) # convert relatedness into a string so there aren't issue with different numbers of individuals
    print(f"there are {len(ind_rel)} relatedness values")
    print(f"done calculating relatedness for sampling time {n}")
    print('Timestamp: {:%H:%M:%S}'.format(datetime.datetime.now()))
    
    # make matrix of pairwise pi for this sampling point
    ind_pi = mts.divergence(ind_nodes, indexes=pairs)
    pi = np.append(pi, str(ind_pi)) # convert relatedness into a string so there aren't issue with different numbers of individuals
    
    print(f"done calculating pi for sampling time {n}")
    print('Timestamp: {:%H:%M:%S}'.format(datetime.datetime.now()))

# Output data for all sampling points
# assemble into dictionary for het data
het_data = { 'pedigree_id':pedigree_id, 
            'het':het, 
            'gen':gen,
            'age':age,
           }

# convert het to dataframe and export
het_df = pd.DataFrame(data=het_data)
het_df.to_csv(outdir+"/"+prefix+"_het.txt", sep=' ', index=True)

# Output Relatedness data for all sampling points
rel_df = pd.DataFrame(data=rel)
rel_df.to_csv(outdir+"/"+prefix+"_relatedness.txt", sep=',', index=True)

# Output pairwise pi data for all sampling points
pi_df = pd.DataFrame(data=pi)
pi_df.to_csv(outdir+"/"+prefix+"_pi.txt", sep=',', index=True)