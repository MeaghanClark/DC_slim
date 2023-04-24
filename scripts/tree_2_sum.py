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
import tskit
np.set_printoptions(threshold=sys.maxsize)

# define custom functions
def getNodes(ids, inds_alive, ts):
    # this function returns a list of nodes from a given list of individuals alive at a specific time from a treesequence
    nodes = []
    for i in ids.to_numpy():
        focal_ind = ts.individual(int(inds_alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        nodes.append(focal_ind.nodes.tolist())
    nodes = [item for sublist in nodes for item in sublist] 
    return nodes

def bootstrapTheta(nodes, niter, nDraws): 
    # returns a list of theta values with length of niter, each value is Wu and Watterson's theta calculated from nDraws nodes
    samples = []
    theta = []
    for i in [*range(0, niter, 1)]: 
        sample = random.choices(nodes, k=nDraws)
        sample_theta = mts.segregating_sites(sample_sets = list(set(sample))) / np.sum([1/i for i in np.arange(1,len(sample))])
        theta.append(sample_theta)
    return theta    
    

def getAlleleCounts(gt_matrix) : 
    # This function returns a matrix of biallelic allele counts and the number of multiallelic loci from a genotype matrix 
    # get mask for multiallelic sites
    multiallelic_mask = (gt_matrix == 2).any(axis=1)
    
    # extract allele counts for multiallelic sites and biallelic sites
    m_gt = gt_matrix[multiallelic_mask]
    bi_gt = gt_matrix[~multiallelic_mask]
        
    # calculate allele counts for biallelic sites
    bi_ac = np.sum(bi_gt, axis=1)

    return(bi_ac, len(m_gt))


def getSFS(gt_matrix, sample_set): 
    # this function returns a site frequency spectrum and the number of multiallelic sites from a genotype matrix and set of sampled nodes
    # based on pyslim manual https://tskit.dev/pyslim/docs/stable/tutorial.html
    output = getAlleleCounts(gt_matrix)
    m_sites = output[1]
    
    # convert the n x 1 array of floats to a vector of integers
    freqs = output[0].flatten().astype(int)
    sfs = np.bincount(freqs, minlength=len(sample_set) + 1)
    return(sfs, m_sites)

def calcPiSFS(afs, mts, m_sites):
    # this function returns pi, calculated from the site frequency spectrum, as defined by Korneliussen et al 2013
    n = len(afs)-1
    list = []
    for i in [*range(1, n-1, 1)]:
        x = i * (n-i) * afs[i]
        list.append(x)
    pi = (sum(list)/math.comb(n, 2))/(mts.sequence_length - m_sites) # need to subtract number of biallelic sites
    return(pi)

def bootstrapPi(nodes, gt_matrix, ts, niter, nDraws): # returns a list of pi values with length of niter, each value is pi calculated from nDraws nodes
    # this function generates bootstrapped replicates of pi from the SFS by sampling nDraws numbers of nodes from a subset with replacement niter times. 
    samples = []
    pi = []
    for i in [*range(0, niter, 1)]: 
        #isolate appropraite columns of the genotype matrix
        sample = random.choices(range(0, len(nodes), 1), k=nDraws)
        sfs = getSFS(gt_matrix[:,sample], sample)
        sample_pi = calcPiSFS(sfs[0], mts, sfs[1])
        pi.append(sample_pi)
    return(pi)

def saveStraps(bootstrapped_data, timepoint, calc_id, dataframe):
    # this function appends bootstrapped_data to a dataframe with specified timepoint and calculation informatio
    rows = []
    row = {'timepoint': timepoint, 'calculation_id': calc_id}
    for i, rep in enumerate(pi_ten):
        row[f'replicate_{i+1}'] = rep
    rows.append(row)
    dataframe = pd.concat([dataframe, pd.DataFrame(rows)], ignore_index=True)
    return(dataframe)

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
orig_ts = tskit.load(treefile)

orig_ts = pyslim.update(orig_ts)

print(f"Loaded tree file")

# recapitate tree
rts = pyslim.recapitate(orig_ts, ancestral_Ne=5487, recombination_rate = 1e-8) # need to change ancestral Ne! 

orig_max_roots = max(t.num_roots for t in orig_ts.trees())
recap_max_roots = max(t.num_roots for t in rts.trees())
print(f"Maximum number of roots before recapitation: {orig_max_roots}\n"
      f"After recapitation: {recap_max_roots}")

# overlay mutations
rate = float(mu)/float(gen_time)
print(f"using rate {rate}")

# mts = pyslim.SlimTreeSequence(msprime.mutate(rts, rate=rate, random_seed = seed, keep=True)) 
mts = msprime.sim_mutations(rts, rate = rate, random_seed = seed, keep = True)
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
elif avg_age == 1: 
    burn = 100000

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
#         - age bins x 2 (upper/lower 10%, upper/lower 5%)
#         - random sampling
#    [2] theta
#         - age cohorts
#         - age bins x 2 (upper/lower 10%, upper/lower 5%)
#         - random sampling
#    [3] decay of LD
#         - age cohorts
#         - age bins x 2 (upper/lower 10%, upper/lower 5%)
#         - random sampling

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


no_straps = 1000

df_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta']) # eventually want to add some measure of LD
df_random_samp_summary = pd.DataFrame(columns = ['timepoint', 'pi_all', 'theta_all']) # add eventually 'pi_ten', 'LD_ten',
df_age_bins = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, (no_straps + 1))])
df_random_samp = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, no_straps + 1)])

# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
# for n in [*range(0, 1, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store pi values that are calculated once per time point
    tp_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta'])
    
    # tp_age_bins = pd.DataFrame(columns = ['timepoint', 'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five'])
    tp_age_bins = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, (no_straps + 1))])

    # tp_random_samp = pd.DataFrame(columns = ['timepoint', 'pi_ten', 'theta_ten', 'pi_all', 'theta_all']) 
    tp_random_samp = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, (no_straps+1))])

    tp_random_samp_summary = pd.DataFrame(columns = ['timepint', 'pi_all', 'theta_all']) # add eventually 'pi_ten', 'LD_ten',

    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    
    # assign timepoint to output file
    #tp_age_bins.loc[0, 'timepoint'] = n 
    #tp_random_samp.loc[0, 'timepoint'] = n 
    tp_random_samp_summary.loc[0, 'timepoint'] = n
    
    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = pyslim.individuals_alive_at(mts, tskit_time)

    # define pedigree ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    x = [mts.individual(i).metadata['pedigree_id'] for i in alive]   

    # create list of individual ages
    meta = metadata[metadata["generation"] == convert_time.iloc[n][1]] # subset metadata for this timepoint
    ages = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["age"])
    all_ages = [item for sublist in ages.to_numpy().tolist() for item in sublist] # list of individual ages for this timepoint

    #### START CALCULATIONS ------------------------------------------------------------------------------------------------------------------------------------------
        
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
    ind_nodes = []    
    for i in samp_pdids.to_numpy():                            # for each individual sampled in slim
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        ind_nodes.append(focal_ind.nodes.tolist())                      # make list of nodes
    print(f"length of ind_nodes is {len(ind_nodes)}")
    all_nodes = [item for sublist in ind_nodes for item in sublist]
    
    ### Random Sampling------------------------------------------------------------------------------------------------------------------------------------------
        # what to fill in: 'theta_ten', 'pi_all', 'theta_all'
        
    # with all_nodes
    tp_random_samp_summary.loc[0, 'pi_all'] = mts.diversity(sample_sets = all_nodes)
    tp_random_samp_summary.loc[0, 'theta_all'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    #tp_random_samp.loc[0, 'LD_all'] = ###
    
    # with ten sampled nodes

    ten_rand_nodes = random.sample(all_nodes, k = 10)
    
    # bootstrap pi and save to tp_random_samp
    pi_ten = bootstrapPi(nodes = ten_rand_nodes, gt_matrix = mts.genotype_matrix(samples = ten_rand_nodes), ts = mts, niter = no_straps, nDraws = 10)
    tp_random_samp = saveStraps(bootstrapped_data = pi_ten, timepoint = n, calc_id = "pi_ten", dataframe = tp_random_samp)

    # bootstrap theta and save to tp_random_samp
    theta_ten = bootstrapTheta(nodes = ten_rand_nodes, niter = no_straps, nDraws = 10)
    tp_random_samp = saveStraps(bootstrapped_data = theta_ten, timepoint = n, calc_id = "theta_ten", dataframe = tp_random_samp)

    
    ### Age Cohorts------------------------------------------------------------------------------------------------------------------------------------------
            # what to fill in: 'age', 'pi', 'theta'
    unique_ages = list(set(all_ages))
    for a in [*range(0, len(unique_ages), 1)]:
        cohort_nodes = [] 
        ids = meta[meta['age'] == unique_ages[a]][["pedigree_id"]]
        for i in ids.to_numpy():
            focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
            cohort_nodes.append(focal_ind.nodes.tolist())   # make list of nodes
        cohort_nodes = [item for sublist in cohort_nodes for item in sublist] # get rid of sub-lists to get overall pi 
        tp_age_cohort.loc[a, 'age'] = unique_ages[a]
        tp_age_cohort.loc[a, 'pi'] = mts.diversity(sample_sets = cohort_nodes)
        tp_age_cohort.loc[a, 'theta'] = mts.segregating_sites(sample_sets = cohort_nodes) / np.sum([1/i for i in np.arange(1,len(cohort_nodes))])
        tp_age_cohort.loc[a, 'timepoint'] = n 
   
    print(f"done with age cohort sampling for sampling point {n} representing tskit time {tskit_time}")

    ### Age Bins------------------------------------------------------------------------------------------------------------------------------------------
        # what to fill in: 'theta_upper_ten', 'theta_lower_ten', 'theta_upper_five', 'theta_lower_five'
        
    # sort individuals by age
    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    ## [1] upper/lower 10%
    
    # define nodes
    
    lower_index_ten = int(len(meta_sorted) * 0.1)
    upper_index_ten = int(len(meta_sorted) * 0.9)

    # Get the lower 10% of data points
    lower_10_percent = meta_sorted[:lower_index_ten]
    
    # Get the upper 10% of data points
    upper_10_percent = meta_sorted[upper_index_ten:]
    
    # lower 10%
    # define nodes 
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_ten_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_lower_ten = bootstrapTheta(nodes = lower_ten_nodes, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = theta_lower_ten, timepoint = n, calc_id = "theta_lower_ten", dataframe = tp_age_bins)

    # bootstrap pi and save to tp_age_bins
    pi_lower_ten = bootstrapPi(nodes = lower_ten_nodes, gt_matrix = mts.genotype_matrix(samples = lower_ten_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = pi_lower_ten, timepoint = n, calc_id = "pi_lower_ten", dataframe = tp_age_bins)

    # upper 10%
    # define nodes 
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_ten_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_upper_ten = bootstrapTheta(nodes = upper_ten_nodes, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = theta_upper_ten, timepoint = n, calc_id = "theta_upper_ten", dataframe = tp_age_bins)

    # bootstrap pi and save to tp_age_bins
    pi_upper_ten = bootstrapPi(nodes = upper_ten_nodes, gt_matrix = mts.genotype_matrix(samples = upper_ten_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = pi_upper_ten, timepoint = n, calc_id = "pi_upper_ten", dataframe = tp_age_bins)
    
    
    ## [2] upper/lower 5%
    
    lower_index_five = int(len(meta_sorted) * 0.05)
    upper_index_five = int(len(meta_sorted) * 0.95)

    # Get the lower 5% of data points
    lower_5_percent = meta_sorted[:lower_index_five]
    
    # Get the upper 5% of data points
    upper_5_percent = meta_sorted[upper_index_five:]
    
    # lower 5%
    # define nodes
    lower_ids = lower_5_percent["pedigree_id"] 
    lower_five_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_lower_five = bootstrapTheta(nodes = lower_five_nodes, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = theta_lower_five, timepoint = n, calc_id = "theta_lower_five", dataframe = tp_age_bins)

    # bootstrap pi and save to tp_age_bins
    pi_lower_five = bootstrapPi(nodes = lower_five_nodes, gt_matrix = mts.genotype_matrix(samples = lower_five_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = pi_lower_five, timepoint = n, calc_id = "pi_lower_five", dataframe = tp_age_bins)

    # upper 5%
    # define nodes
    upper_ids = upper_5_percent["pedigree_id"] 
    upper_five_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_upper_five = bootstrapTheta(nodes = upper_five_nodes, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = theta_upper_five, timepoint = n, calc_id = "theta_upper_five", dataframe = tp_age_bins)

    # bootstrap pi and save to tp_age_bins
    pi_upper_five = bootstrapPi(nodes = upper_five_nodes, gt_matrix = mts.genotype_matrix(samples = upper_five_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_age_bins = saveStraps(bootstrapped_data = pi_upper_five, timepoint = n, calc_id = "pi_upper_five", dataframe = tp_age_bins)

    
    # save output------------------------------------------------------------------------------------------------------------------------------------------

    print(f"done with age bin sampling for sampling point {n} representing tskit time {tskit_time}")

    #### add information from this timepoint to final dataframes

    df_age_cohort = pd.concat([df_age_cohort, tp_age_cohort], axis=0)
    df_age_bins = pd.concat([df_age_bins, tp_age_bins], axis=0)
    df_random_samp = pd.concat([df_random_samp, tp_random_samp], axis=0)

    
df_age_cohort.to_csv(outdir+"/"+prefix+"_age_cohort.txt", sep=',', index=False)
df_age_bins.to_csv(outdir+"/"+prefix+"_age_bins.txt", sep=',', index=False)
df_random_samp.to_csv(outdir+"/"+prefix+"_rand_samp.txt", sep=',', index=False)

print(f"done saving output")
