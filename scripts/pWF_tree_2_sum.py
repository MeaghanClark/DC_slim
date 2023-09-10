#!/usr/bin/env python

# last updates 08/18/2023

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
np.set_printoptions(threshold=sys.maxsize)

# define custom functions ------------------------------------------------------------------------------------------------------------------------------------------

def getNodes(ids, inds_alive, ts):
    # this function returns a list of nodes from a given list of individuals alive at a specific time from a treesequence
    nodes = []
    for i in ids.to_numpy():
        focal_ind = ts.individual(int(inds_alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        nodes.append(focal_ind.nodes.tolist())
    nodes = [item for sublist in nodes for item in sublist] 
    return nodes


def plotDists(g1_data, g2_data, g1_real, g2_real, stat, n): 
    # modified from: https://stackoverflow.com/questions/72931022/how-to-calculate-histogram-intersection

    rng = min(g1_data.min(),g2_data.min()),max(g1_data.max(),g2_data.max()) # define total range of STAT values
    ttest = stats.ttest_rel(g1_data, g2_data, alternative = "less") 
    
    # single pane
    plt.subplots(1, 1, tight_layout=True, sharex=True)
    
    n1, bins1, _ = plt.hist(g1_data, alpha=0.4, bins=100, range=rng, color = "palevioletred")
    n2, bins2, _ = plt.hist(g2_data, alpha=0.4, bins=100, range=rng, color = "lightseagreen")
    
    intersection = np.minimum(n1, n2)
    area = intersection.sum()
    # plot overlap
    plt.bar(bins1[:-1], intersection, width=bins1[1]- bins1[0], color = "mediumpurple")
    
    plt.title(f"overlap of g1 and g2, {stat} at {n}")
    
    if(ttest[1] <= 0.05):
        plt.axvline(x = g1_real, color = 'crimson')
        plt.axvline(x = g2_real, color = 'teal')
    
    print(f"area of overlap is {area} for {stat} at {n}")
    
    
def measureOverlap(g1_data, g2_data): 
    # modified from: https://stackoverflow.com/questions/72931022/how-to-calculate-histogram-intersection

    rng = min(g1_data.min(),g2_data.min()),max(g1_data.max(),g2_data.max()) # define total range of STAT values
    
    # without plotting version
    n1, bins1 = np.histogram(a = g1_data, bins=100, range=rng)
    n2, bins2 = np.histogram(a = g2_data, bins=100, range=rng)
    
    intersection = np.minimum(n1, n2)
    area = intersection.sum()
    return(area)
    
    
def newBoot(group1_nodes, group2_nodes, niter):
    # new bootstrapping function
    
    # input: group 1 nodes, group 2 nodes, number of interations
    
    ## calculate stats using all nodes
    # theta
    real_theta_g1 = mts.segregating_sites(sample_sets = group1_nodes) / np.sum([1/i for i in np.arange(1,len(group1_nodes))])
    real_theta_g2 = mts.segregating_sites(sample_sets = group2_nodes) / np.sum([1/i for i in np.arange(1,len(group2_nodes))])
    
    # pi
    real_pi_g1 = mts.diversity(sample_sets = group1_nodes)
    real_pi_g2 = mts.diversity(sample_sets = group2_nodes)
    
    ## bootstrap stats using genomic windows
    
    # define windows
    num_windows = 1000
    windows = np.linspace(0, mts.sequence_length, num_windows+1)
    window_length = np.diff(windows)
    
    # resampling with replacement is equivalent to multinomial reweighting
    # e.g. matrix-vector product where each matrix row is a draw from a uniform multinomial
    weights = np.random.multinomial(num_windows, [1/num_windows]*num_windows, size=niter)
    
    # theta
    # calculate theta for all genomic windows
    theta_g1 = mts.segregating_sites(sample_sets = group1_nodes, windows = windows, mode = 'site', span_normalise = False) / np.sum([1/i for i in np.arange(1,len(group1_nodes))])
    theta_g2 = mts.segregating_sites(sample_sets = group2_nodes, windows = windows, mode = 'site', span_normalise = False) / np.sum([1/i for i in np.arange(1,len(group2_nodes))])
    
    # matrix multiplication b/w weights and windows gives vector of bootstrap replicates;
    # do this for numerator (sum of pi) and denominator (sum of window length)
    theta_boot_g1 = (weights @ theta_g1) / (weights @ window_length)
    theta_boot_g2 = (weights @ theta_g2) / (weights @ window_length)
        
    # pi 
    # calculate pi for all genomic windows 
    pi_g1 = mts.diversity(sample_sets = group1_nodes, windows=windows, mode='site', span_normalise=False)
    pi_g2 = mts.diversity(sample_sets = group2_nodes, windows=windows, mode='site', span_normalise=False)
    
    pi_boot_g1 = (weights @ pi_g1) / (weights @ window_length)
    pi_boot_g2 = (weights @ pi_g2) / (weights @ window_length)

    # plot for troubleshooting
    # plotDists(theta_boot_g1, theta_boot_g2, real_theta_g1, real_theta_g2, "theta", n)
    # plotDists(pi_boot_g1, pi_boot_g2, real_pi_g1, real_pi_g2, "pi", n)
    
    # measure overlap in distributions 
    
    theta_OVL = measureOverlap(theta_boot_g1, theta_boot_g2)
    pi_OVL = measureOverlap(pi_boot_g1, pi_boot_g2)
    
    # compare bootstrapped distributions using t-test
    
    theta_ttest = stats.ttest_rel(theta_boot_g1, theta_boot_g2, alternative = "less") # ‘less’: the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
    pi_ttest = stats.ttest_rel(pi_boot_g1, pi_boot_g2, alternative = "less") 
    
    
    # return output
    return(real_theta_g1, real_theta_g2, theta_OVL, theta_ttest[1], theta_ttest[0], real_pi_g1, real_pi_g2, pi_OVL, pi_ttest[1], pi_ttest[0]) 

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
rts = pyslim.recapitate(orig_ts, ancestral_Ne=(10000*float(gen_time)), recombination_rate = (1e-8/float(gen_time))) 

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

if avg_age == 1: 
    burn = 100000
else: 
	print(f"hmmm... seems like you should be using a different script!")

print(f"Using a burn-in time of {burn}")

before = range(burn-200, burn, 50)
during = range(burn, burn+50, 5)
after = range(burn+50, burn+550, 50)

cycles = [*before, *during, *after][::-1]

# make conversion df
convert_time = pd.DataFrame({'tskit_time':sampling, 'slim_time':cycles}, columns=['tskit_time', 'slim_time'])

# initalize data lists

no_straps = 1000 # 100 for troubleshooting, 1000 for running

# save some characteristic of the tree sequence for later
seq_length = mts.sequence_length
positions = mts.tables.sites.position


df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta']) #'LD'

# loop through time points to calculate stats using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 2, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store stats values that are calculated once per time point
    
    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta']) # "LD"
        
    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    
    # assign timepoint to output files    
    tp_summary.loc[0, 'timepoint'] = n
    
    # define pedigree ids sampled by slim, representing individuals we have we have age information for
    samp_pdids = metadata[metadata["generation"] == convert_time.iloc[n][1]].filter(["pedigree_id"])
    
    # define tskit ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    alive = pyslim.individuals_alive_at(mts, tskit_time)

    # define pedigree ids of individuals alive at the sampling point in the tree sequences (can be longer than samp_pdids)
    x = [mts.individual(i).metadata['pedigree_id'] for i in alive]   

    #### START CALCULATIONS ------------------------------------------------------------------------------------------------------------------------------------------
        
    # define nodes to work with -- nodes from sampled individuals at this sampling point regardless of age
    ind_nodes = []    
    for i in samp_pdids.to_numpy():                            # for each individual sampled in slim
        focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        ind_nodes.append(focal_ind.nodes.tolist())                      # make list of nodes
    print(f"length of ind_nodes is {len(ind_nodes)}")
    all_nodes = [item for sublist in ind_nodes for item in sublist]
    
    ### Summary stats for entire sample------------------------------------------------------------------------------------------------------------------------------------------
        
    # with all_nodes
    tp_summary.loc[0, 'pi'] = mts.diversity(sample_sets = all_nodes)
    tp_summary.loc[0, 'theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    gt_matrix_all_nodes = mts.genotype_matrix(samples = all_nodes)
    tp_summary.loc[0, 'LD'] = getrSquaredDecayDist(all_nodes, gt_matrix_all_nodes, positions, seq_length)[3]

    
    # save output ------------------------------------------------------------------------------------------------------------------------------------------
    df_summary = pd.concat([df_summary, tp_summary], axis=0)

    # end of for loop

df_summary.to_csv(outdir+"/"+prefix+"_summary.txt", sep=',', index=False)

print(f"done saving output")