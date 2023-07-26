#!/usr/bin/env python

# last updates 07/26/2023

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
import seaborn as sns
import matplotlib.pyplot as plt
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
    for i, rep in enumerate(bootstrapped_data):
        row[f'replicate_{i+1}'] = rep
    rows.append(row)
    dataframe = pd.concat([dataframe, pd.DataFrame(rows)], ignore_index=True)
    return(dataframe)


def getBiGenoMatrix(gt_matrix):
    # This function returns a matrix of biallelic allele counts from a genotype matrix

    # Filter for biallelic loci
    mask = ((gt_matrix != 2).all(axis=1)) & (np.sum(gt_matrix, axis=1) > 0) & (np.sum(gt_matrix, axis=1) < np.shape(gt_matrix)[1])
    bi_gt = gt_matrix[mask]

    return bi_gt


def getVariantPositions(positions, gt_matrix):
    # This function returns the positions corresponding to biallelic allele counts from a genotype matrix

    # Filter for biallelic loci
    mask = ((gt_matrix != 2).all(axis=1)) & (np.sum(gt_matrix, axis=1) > 0) & (np.sum(gt_matrix, axis=1) < np.shape(gt_matrix)[1])
    bi_gt = gt_matrix[mask]
    positions_flt = positions[mask]

    if len(positions_flt) != len(bi_gt):
        print("Length of position object is not the same as the length of the genotype object")
    else:
        return positions_flt
    
    
def getrSquared(bi_gt, A_index, B_index):
    # consider two LOCI at a time 
    # in bi_gt, rows are loci, columns are genotypes

    A_gt = bi_gt[A_index,:]
    B_gt = bi_gt[B_index,:]
    
    genotype_counts = {
        '00': 0,
        '01': 0,
        '10': 0,
        '11': 0
    }
        
    genotypes = np.char.add(A_gt.astype(str), B_gt.astype(str)) # generate genotypes
    
    unique_genotypes, counts = np.unique(genotypes, return_counts=True) # get counts of genotypes 
    
    for genotype, count in zip(unique_genotypes, counts):   # update genotype_counts dict
        genotype_counts[genotype] = count
    
    total = sum(genotype_counts.values()) # number of haplotypes
    
    
    # Calculate allele frequencies
    fA = (genotype_counts["11"] + genotype_counts["10"]) / total
    fB = (genotype_counts["11"] + genotype_counts["01"]) / total
    fa = 1 - fA
    fb = 1 - fB
    F_AB = genotype_counts["11"]/ total
    
    # Calculate DAB
    DAB = F_AB - (fA * fB)
    
    # Calculate r-squared
    r_squared = (DAB**2) / (fA * fa * fB * fb)
    
    if np.isnan(r_squared):
        print("Hmmm.. D is 0? Are you sure you filtered out monomorphic loci?")
    
    return r_squared


def getrSquaredDecayDist(nodes, gt_matrix):

    # get genotype matrix for subset of individuals
    # gt_matrix = mts.genotype_matrix(samples = nodes) # genotype_matrix retains order of nodes passed to it

    # retain only biallelic sites
    bi_gt = getBiGenoMatrix(gt_matrix)

    # downsample loci, randomly for now
    target_loci = random.sample([*range(0, np.shape(bi_gt)[0], 1)], 100)

    # filter bi_gt by target loci
    bi_gt_target = bi_gt[target_loci,] # biallelic loci filtered to retain only target_loci

    # initialize r_squared matrix
    rSquared_matrix = np.empty(shape=(len(target_loci),len(target_loci)))
    rSquared_matrix.fill(-9) # fill with -9s 

    # populate r^2 matrix
    # for each locus.. 
        # calculate r^2 between all pairs of loci
    for i in range(len(target_loci)):
        for j in range(i, len(target_loci)):
            r_squared = getrSquared(bi_gt_target, i, j)
            rSquared_matrix[i,j] = r_squared
            rSquared_matrix[j,i] = r_squared

    # get distance between loci         
    # following code adapted from JesseGarcia562 on tskit github (https://github.com/tskit-dev/tskit/discussions/1118)
    # record positions of loci
    positions = mts.tables.sites.position # this is all positions, need to retain just positions biallelic in sample
    
    # filter positions to match r^2 
    filtered_positions = getVariantPositions(positions, gt_matrix) # filter to retain biallelic sites
    filtered_positions_target = filtered_positions[target_loci] # filter to retain target loci
    
    # Annotate rSquared_matrix
    df=pd.DataFrame(rSquared_matrix)
    df.index=filtered_positions_target
    df.columns=filtered_positions_target

    # Turn into long dataframe
    long_format_distance_df=df.unstack().reset_index()
    long_format_distance_df.columns = ['position_1', 'position_2', 'r_2']
    long_format_distance_df["distance"] = np.fabs(long_format_distance_df["position_1"] - long_format_distance_df["position_2"])
    long_format_distance_df = long_format_distance_df[long_format_distance_df.distance != 0] # remove distance = 0 
    
    # ## Get coordinates for plot/fitting and plot
    # x=long_format_distance_df["distance"]
    # y=long_format_distance_df["r_2"]
    # 
    # plt.scatter(x, y, alpha=0.005)
    
    # define number of bins and bin edges
    num_bins = 10000
    bin_edges = np.linspace(0, mts.sequence_length, num_bins + 1)
    
    # use bins to add distance class column to data frame
    long_format_distance_df['distance_class'] = pd.cut(long_format_distance_df['distance'], bins=bin_edges, right=True, labels = bin_edges[0:len(bin_edges)-1]) #.apply(lambda x: x.left)
    # long_format_distance_df.head()
    
    # calculate mean r2 value for each distance class bin
    mean_r_2 = long_format_distance_df.groupby('distance_class')['r_2'].mean()
    # plot
    #plt.scatter(x = mean_r_2.index, y = mean_r_2, alpha = 0.15)
    
    # calculate mean and standard error of r^2 for last ~20% of bins 
    
    mean_r_2_df = mean_r_2.to_frame().reset_index()
    asym_mean = mean_r_2_df[mean_r_2_df['distance_class'].astype(int) > 80000000]['r_2'].mean()
    
    asym_se = mean_r_2_df[mean_r_2_df['distance_class'].astype(int) > 80000000]['r_2'].sem()
    
    # find distance class at which average r^2 goes within 2 SE of the mean r^2
    
    limit = asym_mean + (asym_se * 2)
    
    distance = mean_r_2_df[mean_r_2_df['r_2'] < limit].iloc[0]['distance_class']
    
    # record distance class, mean r^2 for last 20% and SE for last 20%    
    
    #print(mean_r_2_df.head())
    #print(mean_r_2_df['r_2'].isna().sum())

    return(asym_mean, asym_se, limit, distance)
    
    
def bootstrapLD(ids, nodes, gt_matrix, niter, nDraws): # returns a list of pi values with length of niter, each value is pi calculated from nDraws nodes
    # this function generates bootstrapped replicates of pi from the SFS by sampling nDraws numbers of nodes from a subset with replacement niter times. 
    samples = []
    LD = []
    for i in [*range(0, niter, 1)]: 
        # define sample for bootstrapping
        sample = random.choices(range(0, len(ids), 1), k=nDraws)
        
        # filter genotype matrix to retain only sampled individuals 
        sample_nodes = getNodes(ids = ids.iloc[sample], inds_alive = alive, ts = mts)  # what are the node ids associated with the sample? 
        matching_indices = [nodes.index(node) for node in sample_nodes]  # where to sample nodes match lower_ten_nodes?
        gt_matrix_sample = gt_matrix[:,matching_indices]        # subset gt_matrix

        # calculate LD decay distance 
        sample_LD = getrSquaredDecayDist(sample_nodes, gt_matrix_sample)[3]
        
        LD.append(sample_LD)
    return(LD)


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

no_straps = 1000

df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta', 'LD']) # add eventually 'pi_ten', 'LD_ten',
df_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta']) # eventually want to add some measure of LD
df_sub_samp = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, (no_straps + 1))])


# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 2, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store stats values that are calculated once per time point
    
    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta', "LD"])
    
    # data object to store summary stats for age cohorts
    tp_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'pi', 'theta'])
    
    # data object to store bootstrapped replicates of summary stats for sub-sample and age bins
    tp_sub_samp = pd.DataFrame(columns=['timepoint', 'calculation_id'] + [f'replicate_{i}' for i in range(1, (no_straps + 1))])

    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    
    # assign timepoint to output files    
    tp_summary.loc[0, 'timepoint'] = n
    tp_sub_samp.loc[0, 'timepoint'] = n
    
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
    
    ### Summary stats for entire sample------------------------------------------------------------------------------------------------------------------------------------------
        
    # with all_nodes
    tp_summary.loc[0, 'pi'] = mts.diversity(sample_sets = all_nodes)
    tp_summary.loc[0, 'theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    gt_matrix_all_nodes = mts.genotype_matrix(samples = all_nodes)
    tp_summary.loc[0, 'LD'] = getrSquaredDecayDist(all_nodes, gt_matrix_all_nodes)[3]
    
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
        
    # sort individuals by age
    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    ## [1] upper/lower 10%
    # define nodes from upper and lower 10% of individuals, calculate theta, pi, and LD of entire subsample, then calculate theta, pi and LD of age bins
    
    # define nodes
    
    lower_index_ten = int(len(meta_sorted) * 0.1)
    upper_index_ten = int(len(meta_sorted) * 0.9)

    # Get the lower 10% of individuals by age
    lower_10_percent = meta_sorted[:lower_index_ten]
    
    # Get the upper 10% of individuals by age
    upper_10_percent = meta_sorted[upper_index_ten:]
    
    # lower 10% ------------------------------------------------------------------------------------------------------------------------------------------
    # define nodes 
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_ten_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_lower_ten = bootstrapTheta(nodes = lower_ten_nodes, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = theta_lower_ten, timepoint = n, calc_id = "theta_lower_ten", dataframe = tp_sub_samp)

    # bootstrap pi and save to tp_age_bins
    pi_lower_ten = bootstrapPi(nodes = lower_ten_nodes, gt_matrix = mts.genotype_matrix(samples = lower_ten_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = pi_lower_ten, timepoint = n, calc_id = "pi_lower_ten", dataframe = tp_sub_samp)

    # bootstrap LD and save to tp_age_bins
    # get gt_matrix for all lower nodes
    gt_matrix_lower_ten = mts.genotype_matrix(samples = lower_ten_nodes)
    LD_lower_ten = bootstrapLD(lower_ids, lower_ten_nodes, gt_matrix_lower_ten, niter = no_straps, nDraws = 10) 
    tp_sub_samp= saveStraps(bootstrapped_data = LD_lower_ten, timepoint = n, calc_id = "LD_lower_ten", dataframe = tp_sub_samp)
    

    # upper 10% ------------------------------------------------------------------------------------------------------------------------------------------
    # define nodes 
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_ten_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)
    
    # bootstrap theta and save to tp_age_bins
    theta_upper_ten = bootstrapTheta(nodes = upper_ten_nodes, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = theta_upper_ten, timepoint = n, calc_id = "theta_upper_ten", dataframe = tp_sub_samp)

    # bootstrap pi and save to tp_age_bins
    pi_upper_ten = bootstrapPi(nodes = upper_ten_nodes, gt_matrix = mts.genotype_matrix(samples = upper_ten_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = pi_upper_ten, timepoint = n, calc_id = "pi_upper_ten", dataframe = tp_sub_samp)
    
    # bootstrap LD and save to tp_age_bins
    # get gt_matrix for all lower nodes
    gt_matrix_upper_ten = mts.genotype_matrix(samples = upper_ten_nodes)
    LD_upper_ten = bootstrapLD(upper_ids, upper_ten_nodes, gt_matrix_upper_ten, niter = no_straps, nDraws = 10) 
    tp_sub_samp= saveStraps(bootstrapped_data = LD_upper_ten, timepoint = n, calc_id = "LD_upper_ten", dataframe = tp_sub_samp)
   

    ### No age info------------------------------------------------------------------------------------------------------------------------------------------
    subsamp_ids = pd.concat([lower_ids, upper_ids])
    subsamp_nodes = getNodes(ids = subsamp_ids, inds_alive = alive, ts = mts)

    # bootstrap theta and save to tp_age_bins
    theta_subsamp = bootstrapTheta(nodes = subsamp_nodes, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = theta_subsamp, timepoint = n, calc_id = "theta_subsamp", dataframe = tp_sub_samp)

    # bootstrap pi and save to tp_age_bins
    pi_subsamp = bootstrapPi(nodes = subsamp_nodes, gt_matrix = mts.genotype_matrix(samples = subsamp_nodes), ts = mts, niter = no_straps, nDraws = 20)
    tp_sub_samp = saveStraps(bootstrapped_data = pi_subsamp, timepoint = n, calc_id = "pi_subsamp", dataframe = tp_sub_samp)
    
    # bootstrap LD and save to tp_age_bins
    # get gt_matrix for all lower nodes
    gt_matrix_subsamp= mts.genotype_matrix(samples = subsamp_nodes)
    LD_subsamp = bootstrapLD(subsamp_ids, subsamp_nodes, gt_matrix_subsamp, niter = no_straps, nDraws = 10) 
    tp_sub_samp= saveStraps(bootstrapped_data = LD_subsamp, timepoint = n, calc_id = "LD_subsamp", dataframe = tp_sub_samp)
   

    # save output------------------------------------------------------------------------------------------------------------------------------------------

    print(f"done with age bin sampling for sampling point {n} representing tskit time {tskit_time}")

    #### add information from this timepoint to final dataframes

    df_summary = pd.concat([df_summary, tp_summary], axis=0)
    df_age_cohort = pd.concat([df_age_cohort, tp_age_cohort], axis=0)
    df_sub_samp = pd.concat([df_sub_samp, tp_sub_samp], axis=0)    
    
df_summary.to_csv(outdir+"/"+prefix+"_summary.txt", sep=',', index=False)
df_age_cohort.to_csv(outdir+"/"+prefix+"_age_cohort.txt", sep=',', index=False)
df_sub_samp.to_csv(outdir+"/"+prefix+"_sub_samp.txt", sep=',', index=False)

print(f"done saving output")
