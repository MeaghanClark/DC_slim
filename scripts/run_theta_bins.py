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

# object to store bin distributions 

# [1] oldest 10% vs. youngest 10%
# [2] oldest 5% vs. youngest 5%
# [3] all individuals of max age vs. all individuals of min age?


theta_breakdown = pd.DataFrame(columns = ['timepoint', 'upper_ten', 'lower_ten', 'upper_five', 'lower_five', 'max_age', 'min_age'])

# object to store cohort pi
# pi_by_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'cohort_pi', 'cohort_theta'])

# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 1, 1)]: 

    # initialize data object to store pi values that are calculated once per time point
    bin_tests = pd.DataFrame(columns = ['timepoint', 'upper_ten', 'lower_ten', 'upper_five', 'lower_five', 'max_age', 'min_age'])

    
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

    #### define age bins *ENDED HERE 2/28 * 

    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    #### define windows
    num_windows = 10000

    ######################################################################################################################################################
    #### LOWER/UPPER 10 
    
    lower_index = int(len(meta_sorted) * 0.1)
    upper_index = int(len(meta_sorted) * 0.9)

    # Get the lower 10% of data points
    lower_10_percent = meta_sorted[:lower_index]
    
    # Get the upper 10% of data points
    upper_10_percent = meta_sorted[upper_index:]
    
    # get lists of nodes and generate theta distributions
    
    # lower 10%
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_nodes = []
    for i in lower_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        lower_nodes.append(focal_ind.nodes.tolist())
    lower_nodes = [item for sublist in lower_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = lower_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))
    
    ### calculate distribution of overall thetas from segregating sites
    bin_tests.loc[0, 'lower_ten'] = s_sites / np.sum([1/i for i in np.arange(1,len(lower_nodes))])

    # upper 10%
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_nodes = []
    for i in upper_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        upper_nodes.append(focal_ind.nodes.tolist())
    upper_nodes = [item for sublist in upper_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = upper_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))

    bin_tests.loc[0, 'upper_ten'] = s_sites / np.sum([1/i for i in np.arange(1,len(upper_nodes))])

    ######################################################################################################################################################
    #### LOWER/UPPER 5%
    lower_index = int(len(meta_sorted) * 0.05)
    upper_index = int(len(meta_sorted) * 0.95)

    # Get the lower 5% of data points
    lower_5_percent = meta_sorted[:lower_index]
    
    # Get the upper 5% of data points
    upper_5_percent = meta_sorted[upper_index:]
    
    # lower 5%
    lower_ids = lower_5_percent["pedigree_id"] 
    lower_nodes = []
    for i in lower_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        lower_nodes.append(focal_ind.nodes.tolist())
    lower_nodes = [item for sublist in lower_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = lower_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))

    bin_tests.loc[0, 'lower_five'] = s_sites / np.sum([1/i for i in np.arange(1,len(lower_nodes))])

    # upper 5%
    upper_ids = upper_5_percent["pedigree_id"] 
    upper_nodes = []
    for i in upper_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        upper_nodes.append(focal_ind.nodes.tolist())
    upper_nodes = [item for sublist in upper_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = upper_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))
    bin_tests.loc[0, 'upper_five'] = s_sites / np.sum([1/i for i in np.arange(1,len(upper_nodes))])

    ######################################################################################################################################################
    #### MAX MIN AGES
    meta_sorted
    lower_index = int(len(meta_sorted) * 0.05)
    upper_index = int(len(meta_sorted) * 0.95)

    # Get the lower 10% of data points
    min_age = meta_sorted[meta_sorted["age"] == min(meta_sorted["age"])]
    
    # Get the upper 10% of data points
    max_age = meta_sorted[meta_sorted["age"] == max(meta_sorted["age"])]
    
    print(f"min age is {min_age}; max age is {max_age}")
    print(f"there are {len(meta_sorted[meta_sorted['age'] == min(meta_sorted['age'])])} individuals of min age; there are {len(meta_sorted[meta_sorted['age'] == max(meta_sorted['age'])])} individuals of max age")

    # get lists of nodes and generate theta distributions
    
    # min age
    lower_ids = min_age["pedigree_id"] 
    lower_nodes = []
    for i in lower_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        lower_nodes.append(focal_ind.nodes.tolist())
    lower_nodes = [item for sublist in lower_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = lower_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))
    bin_tests.loc[0, 'min_age'] = s_sites / np.sum([1/i for i in np.arange(1,len(lower_nodes))])

    # max age
    upper_ids = max_age["pedigree_id"] 
    upper_nodes = []
    for i in upper_ids.to_numpy():
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        upper_nodes.append(focal_ind.nodes.tolist())
    upper_nodes = [item for sublist in upper_nodes for item in sublist] 
    s_sites = mts.segregating_sites(sample_sets = upper_nodes, windows = np.linspace(0, mts.sequence_length, num_windows + 1))
    bin_tests.loc[0, 'max_age'] = s_sites / np.sum([1/i for i in np.arange(1,len(upper_nodes))])
    
   ## object to store cohort pi
    theta_breakdown = pd.concat([theta_breakdown, bin_tests], axis = 0)

# save dataframes: pi_data_bins, pi_by_cohort

theta_breakdown.to_csv(outdir+"/"+prefix+"bin_tests.txt", sep=',', index=False)
