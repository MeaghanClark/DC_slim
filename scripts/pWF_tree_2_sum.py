#!/usr/bin/env python

# last updates 2/14/2024

# removed age cohort analyses and added age bin permutations

import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random
import datetime # FOR DEBUGGING 
import csv
import tskit
import scipy.stats as stats
import momi
#import matplotlib.pyplot as plt

np.set_printoptions(threshold=sys.maxsize)

# define custom functions ------------------------------------------------------------------------------------------------------------------------------------------

def getNodes(ids, inds_alive, ts):
    # this function returns a list of nodes from a given list of individuals alive at a specific time from a treesequence
    x = [ts.individual(i).metadata['pedigree_id'] for i in inds_alive]   # define pedigree ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)

    nodes = []
    for i in ids.to_numpy():
        focal_ind = ts.individual(int(inds_alive[np.where(x==i)[0][0]])) # get inidvidual id by matching pedigree id to tskit id
        nodes.append(focal_ind.nodes.tolist())
    nodes = [item for sublist in nodes for item in sublist] 
    return nodes

# ------------------------------------------------------------------------------------------------------------------------------------------

# uncomment these lines when running from command line
# sys.argv = ['tree_processing.py', '../slim_output_11082022/tree_nWF_5_10_89.trees','../slim_output_11082022/metaInd_nWF_5_10_89.txt', '/Users/meaghan/Desktop/DC_slim/het', 's_sites_test', 1e-8, 5, 5]
# arguments: 
# [0] -- python script name
# [1] -- tree file
# [2] -- meta file
# [3] -- outdir
# [4] -- prefix
# [5] -- mu
# [6] -- gen time
# [7] -- avg age
# [8] -- burn in time
# [9] -- rVal

# print some record keeping output
seed = random.randint(1,1e6)
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

avg_age = int(sys.argv[7])
print(f"average age is {avg_age}")
print(f"the type of avg_age is {type(avg_age)}")

burn = int(sys.argv[8])
print(f"burn in time is {burn}")

rVal =  int(sys.argv[9])
print(f"bottleneck intensity is {rVal}")

# read in treefile 
orig_ts = tskit.load(treefile)

orig_ts = pyslim.update(orig_ts)

print(f"Loaded tree file")

# recapitate tree
rts = pyslim.recapitate(orig_ts, ancestral_Ne=(23677*float(gen_time)), recombination_rate = (1e-8/float(gen_time))) 

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
rate = float(mu)/float(gen_time)
print(f"using rate {rate}")

mts = msprime.sim_mutations(rts, rate = rate, random_seed = seed, keep = True)

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

before = range(burn-200, burn, 50)
during = range(burn, burn+50, 5)
after = range(burn+50, burn+550, 50)

cycles = [*before, *during, *after][::-1]

# make conversion df
convert_time = pd.DataFrame({'tskit_time':sampling, 'slim_time':cycles}, columns=['tskit_time', 'slim_time'])

# initalize data lists


# save some characteristic of the tree sequence for later
seq_length = mts.sequence_length
positions = mts.tables.sites.position


df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta'])
df_temporal = pd.DataFrame(columns = ['timepoint', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])
df_permut_temporal = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])

# define sample sizes for temporal sampling, calculated as the mean sample size for each R value from the perennial simulations in get_pWF_sample_sizes.R
if rVal == 2: 
    old_sizes = [111, 113, 115, 112, 115, 111, 113, 114, 112, 111, 113, 112, 112, 113, 114, 113, 114, 111, 111, 111, 112, 115, 116, 115]
    young_sizes = [201, 202, 203, 202, 203, 202, 202, 203, 204, 204, 201, 202, 202, 202, 201, 202, 202, 202, 201, 203, 203, 204, 201, 201]
if rVal == 10: 
    old_sizes = [113, 111, 116, 115, 112, 113, 114, 114, 111, 111, 114, 114, 116, 113, 113, 112, 112, 113, 113, 113, 113, 113, 114, 112]
    young_sizes = [202, 202, 202, 203, 201, 202, 204, 204, 202, 203, 200, 203, 202, 202, 202, 202, 202, 203, 202, 202, 202, 202, 203, 202]
if rVal == 100: 
    old_sizes = [105, 106, 104, 104, 104, 108, 104, 103, 104, 105, 107, 103, 105, 105, 104, 107, 106, 108, 107, 107, 113, 116, 111, 113]
    young_sizes = [185, 186, 185, 186, 186, 185, 186, 186, 186, 186, 187, 187, 185, 185, 185, 186, 185, 184, 186, 186, 203, 202, 203, 201]

old_sizes.reverse()
young_sizes.reverse()



# loop through time points to calculate stats using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 2, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store stats values that are calculated once per time point
    
    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta'])
    tp_temporal = pd.DataFrame(columns = ['timepoint', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])

    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    
    # assign timepoint to output files    
    tp_summary.loc[0, 'timepoint'] = n
    tp_temporal.loc[0,'timepoint'] = n

    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    meta = metadata[metadata["generation"] == convert_time.iloc[n][1]] # subset metadata for this timepoint

    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = pyslim.individuals_alive_at(mts, tskit_time)

    #### START CALCULATIONS ------------------------------------------------------------------------------------------------------------------------------------------
        
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
    
    all_nodes = getNodes(ids = samp_pdids, inds_alive = alive, ts = mts)

    ### Summary stats for entire sample------------------------------------------------------------------------------------------------------------------------------------------
        
    # with all_nodes
    tp_summary.loc[0, 'pi'] = mts.diversity(sample_sets = all_nodes)
    tp_summary.loc[0, 'theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    
    ### Temporal Sampling------------------------------------------------------------------------------------------------------------------------------------------
    # actual values
    # define nodes for sample
    if(n < 23):
    #     if len(meta) > young_sizes[n]:
    #         num_inds_now = young_sizes[n] 
    #     else: 
    #         num_inds_now = len(meta)
    #     if len(meta) > old_sizes[n]:
    #         num_inds_past = old_sizes[n]
    #     else: 
    #         num_inds_past = len(meta)        
        
        num_inds_now = young_sizes[n] 
        num_inds_past = old_sizes[n]

        print(f'number of individuals now: {num_inds_now}')
        print(f'length of meta: {len(meta)}')
        print(f'number of individuals past: {num_inds_past}')
        
        now_sample = meta.sample(n=num_inds_now, replace=False)
        now_nodes = getNodes(ids = now_sample["pedigree_id"], inds_alive = alive, ts = mts)
    
        # define past nodes (timepoint n + 1) 
        meta_past = metadata[metadata["generation"] == convert_time.iloc[n+1][1]] # subset metadata for this timepoint
        
        # remove individuals in now_sample
        # meta_past = meta_past[~meta_past['pedigree_id'].isin(now_sample['pedigree_id'])]
        
        past_sample = meta_past.sample(n = num_inds_past, replace = False)
        past_nodes = getNodes(ids = past_sample["pedigree_id"], inds_alive = pyslim.individuals_alive_at(mts, convert_time.iloc[n+1][0]), ts = mts)

        tp_temporal.loc[0, 'theta_now'] = mts.segregating_sites(sample_sets = now_nodes) / np.sum([1/i for i in np.arange(1,len(now_nodes))])
        tp_temporal.loc[0, 'pi_now'] = mts.diversity(sample_sets = now_nodes)
        
        tp_temporal.loc[0, 'theta_past'] = mts.segregating_sites(sample_sets = past_nodes) / np.sum([1/i for i in np.arange(1,len(past_nodes))])
        tp_temporal.loc[0, 'pi_past'] = mts.diversity(sample_sets = past_nodes)
    
        # permutations      tp_permut_temporal = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])
        temp_meta = pd.concat([meta, meta_past], ignore_index=True)
        
        temp_permuts = [] 
        
        for j in range(1, 101):
        #for j in range(1, 5):
            now_sample = temp_meta.sample(n=num_inds_now, replace=False)
            # remove individuals in now_sample
            # temp_meta_unique = temp_meta[~temp_meta['pedigree_id'].isin(now_sample['pedigree_id'])]

            past_sample = temp_meta.sample(n=num_inds_past, replace=False)
            
            # redo temporal comparisons with permuted timepoints
                    
            # now ------------------------------------------------------------------------------------------------------------------------------------------
            now_nodes = getNodes(ids = now_sample["pedigree_id"], inds_alive = np.concatenate((pyslim.individuals_alive_at(mts, convert_time.iloc[n][0]), pyslim.individuals_alive_at(mts, convert_time.iloc[n+1][0]))), ts = mts)
            
            # past ------------------------------------------------------------------------------------------------------------------------------------------
            past_nodes = getNodes(ids = past_sample["pedigree_id"], inds_alive = np.concatenate((pyslim.individuals_alive_at(mts, convert_time.iloc[n][0]), pyslim.individuals_alive_at(mts, convert_time.iloc[n+1][0]))), ts = mts)
                
            # save output ------------------------------------------------------------------------------------------------------------------------------------------
            pi_now = mts.diversity(sample_sets = now_nodes)
            pi_past = mts.diversity(sample_sets = past_nodes)
    
            theta_now = mts.segregating_sites(sample_sets = now_nodes) / np.sum([1/i for i in np.arange(1,len(now_nodes))])
            theta_past = mts.segregating_sites(sample_sets = past_nodes) / np.sum([1/i for i in np.arange(1,len(past_nodes))])
    
            permut_output = [theta_now, theta_past, pi_now, pi_past]
            permut_output.insert(0, n) 
            permut_output.insert(1, j)
            temp_permuts.append(permut_output)

        tp_permut_temporal = pd.DataFrame(temp_permuts)
        tp_permut_temporal.columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past']

        # save output ------------------------------------------------------------------------------------------------------------------------------------------
        df_temporal = pd.concat([df_temporal, tp_temporal], axis=0)   
        df_permut_temporal = pd.concat([df_permut_temporal, tp_permut_temporal], axis = 0)

    # save output ------------------------------------------------------------------------------------------------------------------------------------------
    df_summary = pd.concat([df_summary, tp_summary], axis=0)
    # end of for loop
    
df_summary.to_csv(outdir+"/"+prefix+"_summary.txt", sep=',', index=False)
df_temporal.to_csv(outdir+"/"+prefix+"_temporal.txt", sep=',', index=False)
df_permut_temporal.to_csv(outdir+"/"+prefix+"_permut_temporal.txt", sep=',', index=False)

print(f"done saving output")
