#!/usr/bin/env python

# last updates 3/7/2024

# combined functionality of pWF and nWF tree processing scripts 

import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random
import csv
import tskit

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
# sys.argv = ['tree_processing.py', '../sample_size_diagnose/tree_pWF_1_100_100.trees','../sample_size_diagnose/metaInd_pWF_1_100_100.txt', '/Users/meaghan/Desktop/DC_slim/sample_size_diagnose', 'dist_test', 1e-8, 1, 1, 100000, 100, "pWF"]
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
# [10] -- model type 

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

model_type =  sys.argv[10]
print(f"model type is is {model_type}")

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

df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta'])
df_temporal = pd.DataFrame(columns = ['timepoint', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])
df_permut_temporal = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])

if model_type == "nWF":
    df_age_bin = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_older_exp', 'pi_younger', 'pi_older', 'pi_older_exp'])
    df_permut_age_bin = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_younger', 'theta_older', 'theta_older_exp', 'pi_younger', 'pi_older', 'pi_older_exp'])
    #df_bin_sample_size = pd.DataFrame(columns = ['young_bin', 'old_bin'])

    # define age bounds, as calculated in define_age_bounds.R
    if avg_age == 2:
        if rVal == 2: 
            upper_bound = [5, 5, 5, 5, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
            upper_bound_exp = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if rVal == 10: 
            upper_bound = [5, 5, 5, 5, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
            upper_bound_exp = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if rVal == 100: 
            upper_bound = [5, 5, 5, 5, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5]
            upper_bound_exp = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if avg_age == 5:
        if rVal == 2: 
            upper_bound = [12, 12, 12, 12, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
            upper_bound_exp = [3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if rVal == 10: 
            upper_bound = [12, 12, 12, 12, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
            upper_bound_exp = [3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        if rVal == 100: 
            upper_bound = [12, 12, 12, 12, 13, 13, 13, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
            upper_bound_exp = [3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
            lower_bound = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    if avg_age == 10:
        if rVal == 2: 
            upper_bound = [23, 23, 23, 24, 25, 25, 25, 25, 25, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 23, 24, 23, 24, 24]
            upper_bound_exp = [7, 7, 7, 7, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
            lower_bound = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        if rVal == 10: 
            upper_bound = [24, 24, 24, 24, 25, 25, 25, 25, 25, 24, 23, 23, 24, 24, 24, 24, 23, 23, 24, 24, 23, 23, 24, 24]
            upper_bound_exp = [7, 7, 7, 7, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
            lower_bound = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        if rVal == 100: 
            upper_bound = [24, 24, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 23, 23, 24, 24, 24, 23, 24, 24, 24, 24, 24]
            upper_bound_exp = [7, 7, 7, 7, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7]
            lower_bound = [1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    if avg_age == 20:
        if rVal == 2: 
            upper_bound = [46, 46, 46, 47, 48, 48, 48, 47, 48, 47, 47, 47, 47, 47, 47, 47, 46, 46, 47, 47, 47, 47, 46, 47]
            upper_bound_exp = [14, 14, 14, 14, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
            lower_bound = [2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        if rVal == 10: 
            upper_bound = [47, 46, 47, 46, 48, 47, 48, 48, 48, 47, 48, 47, 47, 47, 47, 47, 47, 46, 46, 46, 46, 46, 47, 46]
            upper_bound_exp = [14, 14, 14, 14, 15, 15, 15, 14, 14, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
            lower_bound = [2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        if rVal == 100: 
            upper_bound = [47, 47, 46, 47, 47, 47, 47, 48, 47, 48, 47, 47, 48, 47, 47, 47, 47, 47, 46, 47, 47, 47, 46, 46]
            upper_bound_exp = [14, 14, 14, 14, 15, 15, 15, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14]
            lower_bound = [2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
    
    upper_bound.reverse()
    upper_bound_exp.reverse()
    lower_bound.reverse()

if model_type == "pWF":
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

    
# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
# for n in [*range(22, 24 , 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")

    # initialize data object to store stats values that are calculated once per time point

    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta']) 

    # data object to store temporal stats
    tp_temporal = pd.DataFrame(columns = ['timepoint', 'theta_now', 'theta_past', 'pi_now', 'pi_past']) # newBoot(future_nodes, now_nodes, niter = no_straps)
    tp_permut_temporal = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])
    
    # assign timepoint to output files    
    tp_summary.loc[0, 'timepoint'] = n
    tp_temporal.loc[0, 'timepoint'] = n 

    # data object to store age bin stats
    if model_type == "nWF":
        tp_age_bin = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_older_exp', 'pi_younger', 'pi_older', 'pi_older_exp'])
        tp_permut_age_bin = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_younger', 'theta_older', 'theta_older_exp', 'pi_younger', 'pi_older', 'pi_older_exp'])
        #tp_bin_sample_size = pd.DataFrame(columns = ['timepoint', 'young_bin', 'old_bin'])

        tp_age_bin.loc[0, 'timepoint'] = n
        #tp_bin_sample_size.loc[0, 'timepoint'] = n 

    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = pyslim.individuals_alive_at(mts, tskit_time)

    meta = metadata[metadata["generation"] == convert_time.iloc[n][1]] # subset metadata for this timepoint
    
    #### START CALCULATIONS ------------------------------------------------------------------------------------------------------------------------------------------
        
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
   
    all_nodes = getNodes(ids = samp_pdids, inds_alive = alive, ts = mts)
    
    ### Summary stats for entire sample------------------------------------------------------------------------------------------------------------------------------------------
        
    # with all_nodes
    tp_summary.loc[0, 'pi'] = mts.diversity(sample_sets = all_nodes)
    tp_summary.loc[0, 'theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    
    ### Age Bins------------------------------------------------------------------------------------------------------------------------------------------
    if model_type == "nWF":
        # create list of individual ages
        ages = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["age"])
        all_ages = [item for sublist in ages.to_numpy().tolist() for item in sublist] # list of individual ages for this timepoint

        # sort individuals by age
        meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
        
        # upper bin ------------------------------------------------------------------------------------------------------------------------------------------
        upper_ids = meta_sorted[meta_sorted["age"] >= upper_bound[n]]["pedigree_id"]
        upper_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)

        # lower bin ------------------------------------------------------------------------------------------------------------------------------------------
        lower_ids = meta_sorted[meta_sorted["age"] <= lower_bound[n]]["pedigree_id"][0:len(upper_ids)] # keep sample sizes the same between upper and lower bins
        lower_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
        
        # record sample sizes for pWF processing
        #tp_bin_sample_size.loc[0, 'young_bin'] = len(lower_nodes)
        #tp_bin_sample_size.loc[0, 'old_bin'] = len(upper_nodes)
    
        # expanded upper bin ------------------------------------------------------------------------------------------------------------------------------------------
        upper_exp_ids = meta_sorted[meta_sorted["age"] >= upper_bound_exp[n]]["pedigree_id"][(len(meta_sorted[meta_sorted["age"] >= upper_bound_exp[n]]) - len(upper_ids)):len(meta_sorted[meta_sorted["age"] >= upper_bound_exp[n]])]
        upper_exp_nodes = getNodes(ids = upper_exp_ids, inds_alive = alive, ts = mts)
        
        #print(f"There are {len(lower_nodes)} nodes in the young group and {len(upper_nodes)} nodes in the old group and {len(upper_exp_nodes)} in the expanded old group")
    
        tp_age_bin.loc[0, 'pi_younger'] = mts.diversity(sample_sets = lower_nodes)
        tp_age_bin.loc[0, 'pi_older'] = mts.diversity(sample_sets = upper_nodes)
        tp_age_bin.loc[0, 'pi_older_exp'] = mts.diversity(sample_sets = upper_exp_nodes)
                                                                                    
        tp_age_bin.loc[0, 'theta_younger'] = mts.segregating_sites(sample_sets = lower_nodes) / np.sum([1/i for i in np.arange(1,len(lower_nodes))])
        tp_age_bin.loc[0, 'theta_older'] = mts.segregating_sites(sample_sets = upper_nodes) / np.sum([1/i for i in np.arange(1,len(upper_nodes))])
        tp_age_bin.loc[0, 'theta_older_exp'] = mts.segregating_sites(sample_sets = upper_exp_nodes) / np.sum([1/i for i in np.arange(1,len(upper_exp_nodes))])
    
        age_permut = all_ages.copy()
        meta_permut = meta.copy()
        permuts = [] 
    
        for j in range(1, 101):
        # for j in range(1, 5):
            random.shuffle(age_permut)
            meta_permut['age_permut'] = age_permut
    
            meta_permut_sorted = meta_permut.sort_values(by='age_permut', ascending=True) # sort metadata by age
            
            # redo age bins with permuted ages             

            # upper nodes ------------------------------------------------------------------------------------------------------------------------------------------
            upper_ids = meta_permut_sorted[meta_permut_sorted["age_permut"] >= upper_bound[n]]["pedigree_id"]
            upper_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
        
            # lower nodes ------------------------------------------------------------------------------------------------------------------------------------------
            lower_ids = meta_permut_sorted[meta_permut_sorted["age_permut"] <= lower_bound[n]]["pedigree_id"][0:len(upper_ids)]
            lower_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
        
            # expanded upper nodes ------------------------------------------------------------------------------------------------------------------------------------------
            upper_exp_ids = meta_permut_sorted[meta_permut_sorted["age"] >= upper_bound_exp[n]]["pedigree_id"][(len(meta_permut_sorted[meta_permut_sorted["age"] >= upper_bound_exp[n]]) - len(upper_ids)):len(meta_sorted[meta_sorted["age"] >= upper_bound_exp[n]])]
            upper_exp_nodes = getNodes(ids = upper_exp_ids, inds_alive = alive, ts = mts)
            
            #print(f"There are {len(lower_nodes)} nodes in the young group and {len(upper_nodes)} nodes in the old group and {len(upper_exp_nodes)} in the expanded old group for permutations")
            
            # save output ------------------------------------------------------------------------------------------------------------------------------------------
            pi_younger = mts.diversity(sample_sets = lower_nodes)
            pi_older = mts.diversity(sample_sets = upper_nodes)
            pi_older_exp = mts.diversity(sample_sets = upper_exp_nodes)
    
            theta_younger = mts.segregating_sites(sample_sets = lower_nodes) / np.sum([1/i for i in np.arange(1,len(lower_nodes))])
            theta_older = mts.segregating_sites(sample_sets = upper_nodes) / np.sum([1/i for i in np.arange(1,len(upper_nodes))])
            theta_older_exp = mts.segregating_sites(sample_sets = upper_exp_nodes) / np.sum([1/i for i in np.arange(1,len(upper_exp_nodes))])
    
            permut_output = [theta_younger, theta_older, theta_older_exp, pi_younger, pi_older, pi_older_exp]
            permut_output.insert(0, n) 
            permut_output.insert(1, j)
            permuts.append(permut_output)
        
        tp_permut_age_bin = pd.DataFrame(permuts)
        tp_permut_age_bin.columns = ['timepoint', 'permutation', 'theta_younger', 'theta_older', 'theta_older_exp', 'pi_younger', 'pi_older', 'pi_older_exp']
    
    ### Temporal Sampling------------------------------------------------------------------------------------------------------------------------------------------
    # actual values 
    # define nodes for sample
    if(n < 23):
        if model_type == "nWF":
            num_inds_now = len(upper_nodes) # keep sample size consistent
            num_inds_past = len(upper_nodes)
        if model_type == "pWF": 
            num_inds_now = old_sizes[n] # keep sample size consistent
            num_inds_past = old_sizes[n]
    
        # define past nodes (timepoint n + 1) 
        meta_past = metadata[metadata["generation"] == convert_time.iloc[n+1][1]] # subset metadata for this timepoint
        # remove individuals in now_sample
        meta_past_unique = meta_past[~meta_past['pedigree_id'].isin(now_sample['pedigree_id'])]
        
        past_sample = meta_past_unique.sample(n = num_inds_past, replace = False)
        past_nodes = getNodes(ids = past_sample["pedigree_id"], inds_alive = pyslim.individuals_alive_at(mts, convert_time.iloc[n+1][0]), ts = mts)
        
        tp_temporal.loc[0, 'theta_now'] = mts.segregating_sites(sample_sets = now_nodes) / np.sum([1/i for i in np.arange(1,len(now_nodes))])
        tp_temporal.loc[0, 'pi_now'] = mts.diversity(sample_sets = now_nodes)
        
        tp_temporal.loc[0, 'theta_past'] = mts.segregating_sites(sample_sets = past_nodes) / np.sum([1/i for i in np.arange(1,len(past_nodes))])
        tp_temporal.loc[0, 'pi_past'] = mts.diversity(sample_sets = past_nodes)
    
        # permutations      tp_permut_temporal = pd.DataFrame(columns = ['timepoint', 'permutation', 'theta_now', 'theta_past', 'pi_now', 'pi_past'])
        temp_meta = pd.concat([meta, meta_past], ignore_index=True)
        # remove duplicate individuals 
        temp_meta = temp_meta.drop_duplicates(subset = ['pedigree_id'], keep = 'first')

        temp_permuts = [] 
    
        for j in range(1, 101):
        # for j in range(1, 5):
            now_sample = temp_meta.sample(n=num_inds_now, replace=False)
            # remove individuals in now_sample
            #temp_meta_unique = temp_meta[~temp_meta['pedigree_id'].isin(now_sample['pedigree_id'])]

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
    if model_type == "nWF":
        df_age_bin = pd.concat([df_age_bin, tp_age_bin], axis=0)   
        df_permut_age_bin = pd.concat([df_permut_age_bin, tp_permut_age_bin], axis = 0)
        # df_bin_sample_size = pd.concat([df_bin_sample_size, tp_bin_sample_size], axis = 0)
    # end of for loop

df_summary.to_csv(outdir+"/"+prefix+"_summary.txt", sep=',', index=False)
df_temporal.to_csv(outdir+"/"+prefix+"_temporal.txt", sep=',', index=False)
df_permut_temporal.to_csv(outdir+"/"+prefix+"_permut_temporal.txt", sep=',', index=False)

if model_type == "nWF":
    df_age_bin.to_csv(outdir+"/"+prefix+"_age_bin.txt", sep=',', index=False)
    df_permut_age_bin.to_csv(outdir+"/"+prefix+"_permut_age_bin.txt", sep=',', index=False)
    # df_bin_sample_size.to_csv(outdir+"/"+prefix+"_bin_sample_size.txt", sep=',', index=False)

print(f"done saving output")
