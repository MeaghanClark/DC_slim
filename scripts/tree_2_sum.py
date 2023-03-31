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

# define custom functions
def getNodes(ids, inds_alive, ts):
    nodes = []
    for i in ids.to_numpy():
        focal_ind = ts.individual(int(inds_alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        nodes.append(focal_ind.nodes.tolist())
    nodes = [item for sublist in nodes for item in sublist] 
    return nodes

def bootstrapTheta(nodes, niter, nDraws): # returns a list of theta values with length of niter, each value is Wu and Watterson's theta calculated from nDraws nodes
    samples = []
    theta = []
    for i in [*range(0, niter, 1)]: 
        sample = random.choices(nodes, k=nDraws)
        sample_theta = mts.segregating_sites(sample_sets = list(set(sample))) / np.sum([1/i for i in np.arange(1,len(sample))])
        theta.append(sample_theta)
    return theta

## def bootstrapPi(nodes, niter, nDraws): # returns a list of pi values with length of niter, each value is pi calculated from nDraws nodes
##     samples = []
##     pi = []
##     for i in [*range(0, niter, 1)]: 
##         sample = random.choices(nodes, k=nDraws)
##         sample_pi = mts.diversity(sample_sets = sample) # broken... doens't work with duplicated nodes
##         pi.append(sample_pi)
##     return pi

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

# data we are extracting: 
#    [1] pairwise pi 
#         - age cohorts 
#         - age bins x 3 (upper/lower 10%, upper/lower 5%, min/max) -- on hold for now! 
#         - random sampling
#    [2] theta
#         - age cohorts
#         - age bins x 2 (upper/lower 10%, upper/lower 5%)
#         - random sampling
#    [3] decay of LD
#         - age cohorts
#         - age bins x 2 (upper/lower 10%, upper/lower 5%)
#         - random sampling

# For now, I am just calculating theta... 

# data files: 
#    [1] age cohorts file
#         - columns = ['timepoint', 'age', 'pi', 'theta', 'LD']
#         - rows = calculations different age cohorts from 24 timepoints 
#    [2] age bins file
#         - columns = ['timepoint', 'pi_upper_ten', 'pi_lower_ten', 'pi_upper_five', 'pi_lower_five', 'pi_max_age', 'pi_min_age', 
#                      'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five', 'theta_max_age', 'theta_min_age', 
#                      'LD_upper_ten', 'LD_lower_ten', 'LD_upper_five', 'LD_lower_five', 'LD_max_age', 'LD_min_age']
#         - rows = calculations from 24 timepoints 
#    [3] random sampling file 
#         - columns = ['timepoint', 'pi', 'theta', 'LD']
#         - rows = calculations from 24 timepoints 



df_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta']) # eventually want to add some measure of LD
df_age_bins = pd.DataFrame(columns = ['timepoint', # 'pi_upper_ten', 'pi_lower_ten', 'pi_upper_five', 'pi_lower_five',
                       'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five'])
                       #'LD_upper_ten', 'LD_lower_ten', 'LD_upper_five', 'LD_lower_five'])
df_random_samp = pd.DataFrame(columns = ['timepoint', 'pi_ten', 'theta_ten', 'pi_all', 'theta_all']) # add eventually 'pi_ten', 'LD_ten',


# object to store cohort pi
# pi_by_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'cohort_pi', 'cohort_theta'])

# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 1, 1)]: 

    # initialize data object to store pi values that are calculated once per time point
    tp_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta'])
    tp_age_bins = pd.DataFrame(columns = ['timepoint', 'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five'])
    tp_random_samp = pd.DataFrame(columns = ['timepoint', 'pi_ten', 'theta_ten', 'pi_all', 'theta_all']) 
    
    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    # assign timepoint to output file
    tp_age_cohort.loc[0, 'timepoint'] = n 
    tp_age_bins.loc[0, 'timepoint'] = n 
    tp_random_samp.loc[0, 'timepoint'] = n 
    
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

    #### START CALCUATIONS
    
    # define number of bootstraps
    no_straps = 1000
        
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
    ind_nodes = []    
    for i in samp_pdids.to_numpy():                            # for each individual sampled in slim
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        ind_nodes.append(focal_ind.nodes.tolist())                      # make list of nodes
    print(f"length of ind_nodes is {len(ind_nodes)}")
    all_nodes = [item for sublist in ind_nodes for item in sublist]

    ### Random Sampling 
        # what to fill in: 'theta_ten', 'pi_all', 'theta_all'
    # with all_nodes
    tp_random_samp.loc[0, 'pi_all'] = mts.diversity(sample_sets = all_nodes)
    tp_random_samp.loc[0, 'theta_all'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    #tp_random_samp.loc[0, 'LD_all'] = ###
    
    # with ten sampled nodes
    ten_rand_nodes = random.sample(all_nodes, k = 10)

    tp_random_samp.loc[0, 'pi_ten'] = mts.diversity(sample_sets = ten_rand_nodes)
    tp_random_samp.loc[0, 'theta_ten'] = bootstrapTheta(ten_rand_nodes, no_straps, 10)
    
    
    #tp_random_samp.loc[0, 'LD_ten'] = ###
    
    print(f"done with random sampling for sampling point {n} representing tskit time {tskit_time}")

    ### Age Cohorts
            # what to fill in: 'age', 'pi', 'theta'
    unique_ages = list(set(all_ages))
    for a in [*range(0, len(unique_ages), 1)]:
        cohort_nodes = [] 
        ids = meta[meta['age'] == unique_ages[a]][["pedigree_id"]]
        for i in ids.to_numpy():
            focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
            cohort_nodes.append(focal_ind.nodes.tolist())   # make list of nodes
        cohort_nodes = [item for sublist in cohort_nodes for item in sublist] # get rid of sub-lists to get overall pi 
        tp_age_cohort.loc[0, 'age'] = unique_ages[a]
        tp_age_cohort.loc[0, 'pi'] = mts.diversity(sample_sets = cohort_nodes)
        tp_age_cohort.loc[0, 'theta'] = mts.segregating_sites(sample_sets = cohort_nodes) / np.sum([1/i for i in np.arange(1,len(cohort_nodes))])
        #tp_age_cohort.loc[0, 'LD'] =
   
    print(f"done with age cohort sampling for sampling point {n} representing tskit time {tskit_time}")

    ### Age Bins
        # what to fill in: 'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five'
        
    # sort individuals by age
    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    ## [1] upper/lower 10%
    
    # define nodes
    
    lower_index_ten = int(len(meta_sorted) * 0.1)
    upper_index_ten = int(len(meta_sorted) * 0.9)

    # Get the lower 10% of data points
    lower_10_percent = meta_sorted[:lower_index]
    
    # Get the upper 10% of data points
    upper_10_percent = meta_sorted[upper_index:]
    
    # lower 10%
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_ten_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
    tp_age_bins.loc[0, 'theta_lower_ten'] = bootstrapTheta(lower_ten_nodes, no_straps, 20)
    
    # upper 10%
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_ten_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
    tp_age_bins.loc[0, 'theta_upper_ten'] = bootstrapTheta(upper_ten_nodes, no_straps, 20)
    
    ####
    ###
    ##
    #
    # stopped here because bootstrapPi isn't working b/c of duplicate node issues
    #
    ##
    ###
    ####
    
    # save output 
    ## tp_age_bins.loc[0, 'pi_upper_ten'] = 
    ## tp_age_bins.loc[0, 'pi_lower_ten'] = 
    ## tp_age_bins.loc[0, 'LD_upper_ten'] = 
    ## tp_age_bins.loc[0, 'LD_lower_ten'] = 

    
    ## [2] upper/lower 5%
    
    lower_index_ten = int(len(meta_sorted) * 0.05)
    upper_index_ten = int(len(meta_sorted) * 0.95)

    # Get the lower 5% of data points
    lower_5_percent = meta_sorted[:lower_index]
    
    # Get the upper 5% of data points
    upper_5_percent = meta_sorted[upper_index:]
    
    # lower 5%
    lower_ids = lower_5_percent["pedigree_id"] 
    lower_five_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
    tp_age_bins.loc[0, 'theta_lower_five'] = bootstrapTheta(lower_five_nodes, no_straps, 10)
    
    # upper 10%
    upper_ids = upper_5_percent["pedigree_id"] 
    upper_five_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
    tp_age_bins.loc[0, 'theta_upper_five'] = bootstrapTheta(upper_five_nodes, no_straps, 10)

    
    # save output 
    ## tp_age_bins.loc[0, 'pi_upper_five'] = 
    ## tp_age_bins.loc[0, 'pi_lower_five'] = 
    ## tp_age_bins.loc[0, 'LD_upper_five'] = 
    ## tp_age_bins.loc[0, 'LD_lower_five'] = 

    ## [3] min/max -- leave for now 

    ### how many bootstrap replicates to do? 1000
    ### how many samples to take? sample same number of nodes with replacement, how the group looks by chance, otherwise inflating sampling variance 
    print(f"done with age bin sampling for sampling point {n} representing tskit time {tskit_time}")

    #### add information from this timepoint to final dataframes

    df_age_cohort = pd.concat([df_age_cohort, tp_age_cohort], axis=0)
    df_age_bins = pd.concat([df_age_bins, tp_age_bins], axis=0)
    df_random_samp = pd.concat([df_random_samp, tp_random_samp], axis=0)

df_age_cohort.to_csv(outdir+"/"+prefix+"_age_cohort.txt", sep=',', index=False)
df_age_bins.to_csv(outdir+"/"+prefix+"_age_bins.txt", sep=',', index=False)
df_random_samp.to_csv(outdir+"/"+prefix+"_rand_samp.txt", sep=',', index=False)