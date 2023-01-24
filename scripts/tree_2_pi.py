#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd
import random
import datetime # FOR DEBUGGING 
import csv
import math
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

avg_age = int(sys.argv[7])
print(f"average age is {avg_age}")
print(f"the type of avg_age is {type(avg_age)}")

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
# debug line: 

print(f"average age is {avg_age}")
print(f"the type of avg_age is {type(avg_age)}")

if avg_age == 2: 
    burn = 126300
elif avg_age == 5: 
    burn = 315900  # changed b/c of my slim error with % 50 
elif avg_age == 10: 
    burn = 635600
elif avg_age == 20: 
    burn = 1264600

print(f"Using a burn-in time of {burn}")

before = range(burn-200, burn, 50)
during = range(burn, burn+50, 5)
after = range(burn+50, burn+550, 50)

cycles = [*before, *during, *after][::-1]

# make conversion df
convert_time = pd.DataFrame({'tskit_time':sampling, 'slim_time':cycles}, columns=['tskit_time', 'slim_time'])

# initalize data lists

# object to store overall pi and upper and lower 10% pi
pi_data_bins = pd.DataFrame(columns = ['timepoint', 'overall_pi', 'lower_pi', 'upper_pi', 'overall_theta', 'lower_theta', 'upper_theta', 'theta_bootstraps'])

# object to store cohort pi
pi_by_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'cohort_pi', 'cohort_theta'])

# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 

    # initialize data object to store pi values that are calculated once per time point
    binned_pi = pd.DataFrame(columns = ['timepoint', 'overall_pi', 'lower_pi', 'upper_pi', 'overall_theta', 'lower_theta', 'upper_theta', 'theta_bootstraps'])
    
    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")

    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = mts.individuals_alive_at(tskit_time)

    # define pedigree ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    x = [mts.individual(i).metadata['pedigree_id'] for i in alive]   

    # create list of individual ages
    meta = metadata[metadata["generation"] == convert_time.iloc[n][1]] # subset metadata for this timepoint
    ages = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["age"])
    all_ages = [item for sublist in ages.to_numpy().tolist() for item in sublist] # list of individual ages for this timepoint

    #### measure overall pi at sampling point
    
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
    ind_nodes = []    
    for i in samp_pdids.to_numpy():                            # for each individual sampled in slim
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        ind_nodes.append(focal_ind.nodes.tolist())                      # make list of nodes
    print(f"length of ind_nodes is {len(ind_nodes)}")
    all_nodes = [item for sublist in ind_nodes for item in sublist]
    binned_pi.loc[0, 'overall_pi'] = mts.diversity(sample_sets = all_nodes) # calculate pi and save to dataframe
    binned_pi.loc[0, 'overall_theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])

    #### measure overall pi within age cohorts 
    
    unique_ages = list(set(all_ages))
    
    cohort_pi_df = pd.DataFrame(index=[*range(0, len(unique_ages), 1)], columns = ['timepoint', 'age', 'cohort_pi', 'cohort_theta'])
    pi_from_cohorts = [] 
    for a in [*range(0, len(unique_ages), 1)]:
        cohort_nodes = [] 
        ids = meta[meta['age'] == unique_ages[a]][["pedigree_id"]]
        for i in ids.to_numpy():
            focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
            cohort_nodes.append(focal_ind.nodes.tolist())   # make list of nodes
        cohort_nodes = [item for sublist in cohort_nodes for item in sublist] # get rid of sub-lists to get overall pi 
        cohort_pi_df.loc[a, 'timepoint'] = n
        cohort_pi_df.loc[a, 'age'] = unique_ages[a] # record focal age in dataframe
        cohort_pi_df.loc[a, 'cohort_pi'] = mts.diversity(sample_sets = cohort_nodes) # calculate pi and save to dataframe
        cohort_pi_df.loc[a, 'cohort_theta'] = mts.segregating_sites(sample_sets = cohort_nodes) / np.sum([1/i for i in np.arange(1,len(cohort_nodes))])

    # end product: cohort_pi_df

    #### measure overall pi within age bins

    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    lower_index = int(len(meta_sorted) * 0.1)
    upper_index = int(len(meta_sorted) * 0.9)

    # Get the lower 10% of data points
    lower_10_percent = meta_sorted[:lower_index]
    
    # Get the upper 10% of data points
    upper_10_percent = meta_sorted[upper_index:]
    
    # get lists of nodes
    
    # lower 10%
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_nodes = []
    for i in lower_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        lower_nodes.append(focal_ind.nodes.tolist())
    lower_nodes = [item for sublist in lower_nodes for item in sublist] 
    binned_pi.loc[0, 'lower_pi'] = mts.diversity(sample_sets = lower_nodes) # calculate pi and save to dataframe
    binned_pi.loc[0, 'lower_theta'] = mts.segregating_sites(sample_sets = lower_nodes) / np.sum([1/i for i in np.arange(1,len(lower_nodes))])
    
    # upper 10%
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_nodes = []
    for i in upper_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        upper_nodes.append(focal_ind.nodes.tolist())
    upper_nodes = [item for sublist in upper_nodes for item in sublist] 
    binned_pi.loc[0, 'upper_pi'] = mts.diversity(sample_sets = upper_nodes) # calculate pi and save to dataframe
    binned_pi.loc[0, 'upper_theta'] = mts.segregating_sites(sample_sets = upper_nodes) / np.sum([1/i for i in np.arange(1,len(upper_nodes))])

    #### bootstrapping
    bootstraps = []
    for i in [*range(0, 100, 1)]: 
        # sample nodes
        lower_sample = random.sample(lower_nodes, 4)
        upper_sample = random.sample(upper_nodes, 4)
        # calculate theta
        lower_theta = mts.segregating_sites(sample_sets = lower_sample) / np.sum([1/i for i in np.arange(1,len(lower_sample))])
        upper_theta = mts.segregating_sites(sample_sets = upper_sample) / np.sum([1/i for i in np.arange(1,len(upper_sample))])
        # take difference
        difference = lower_theta - upper_theta
        # add to data storage object
        bootstraps.append(difference)
    # add to dataframe
    binned_pi.loc[0,'theta_bootstraps'] = bootstraps
    
    #### add information from this timepoint to final dataframes
    binned_pi.loc[0, 'timepoint'] = n

    pi_data_bins = pd.concat([pi_data_bins, binned_pi], axis=0)

    ## object to store cohort pi
    pi_by_cohort = pd.concat([pi_by_cohort, cohort_pi_df], axis = 0)

# save dataframes: pi_data_bins, pi_by_cohort

pi_data_bins.to_csv(outdir+"/"+prefix+"_pi_bins.txt", sep=',', index=False)
pi_by_cohort.to_csv(outdir+"/"+prefix+"_pi_cohort.txt", sep=',', index=False)
