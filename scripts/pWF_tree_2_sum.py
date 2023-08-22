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
from numba import jit
np.set_printoptions(threshold=sys.maxsize)

# define custom functions ------------------------------------------------------------------------------------------------------------------------------------------

@jit(nopython=True)    
def factorial(x):
    return math.exp(math.lgamma((x+1))

                    
@jit(nopython=True)    
def binomialCoeff(n, k): 
    coeff = factorial(n) / (factorial(k) * factorial(n-k))
    return(coeff)


def getNodes(ids, inds_alive, ts):
    # this function returns a list of nodes from a given list of individuals alive at a specific time from a treesequence
    nodes = []
    for i in ids.to_numpy():
        focal_ind = ts.individual(int(inds_alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
        nodes.append(focal_ind.nodes.tolist())
    nodes = [item for sublist in nodes for item in sublist] 
    return nodes


@jit(nopython=True)
def getAlleleCounts(gt_matrix) : 
    # This function returns a matrix of biallelic allele counts and the number of multiallelic loci from a genotype matrix 
    # get mask for multiallelic sites
    multiallelic_mask = (np.sum(gt_matrix == 2, axis = 1) > 0)
    
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


@jit(nopython=True)
def calcPiSFS(afs, m_sites, seq_length):
    # this function returns pi, calculated from the site frequency spectrum, as defined by Korneliussen et al 2013
    n = len(afs)-1
    list = []
    for i in [*range(1, n-1, 1)]:
        x = i * (n-i) * afs[i]
        list.append(x)
    pi = (sum(list)/binomialCoeff(n, 2))/(seq_length - m_sites) # need to subtract number of biallelic sites
    return(pi)


@jit(nopython=True)
def getBiGenoMatrix(gt_matrix):
    # This function returns a matrix of biallelic allele counts from a genotype matrix

    # Filter for biallelic loci
    mask = ((np.sum(gt_matrix !=2, axis=1) == np.shape(gt_matrix)[1]) & (np.sum(gt_matrix, axis=1) > 0) & (np.sum(gt_matrix, axis=1) < np.shape(gt_matrix)[1]))
    bi_gt = gt_matrix[mask]

    return bi_gt


@jit(nopython=True)
def getVariantPositions(positions, gt_matrix):
    # This function returns the positions corresponding to biallelic allele counts from a genotype matrix

    # Filter for biallelic loci
    mask = ((np.sum(gt_matrix !=2, axis =1) == np.shape(gt_matrix)[1]) & (np.sum(gt_matrix, axis=1) > 0) & (np.sum(gt_matrix, axis=1) < np.shape(gt_matrix)[1]))
    bi_gt = gt_matrix[mask]
    positions_flt = positions[mask]

    if len(positions_flt) != len(bi_gt):
        print("Length of position object is not the same as the length of the genotype object")
    else:
        return positions_flt
    
    
@jit(nopython=True)    
def getrSquared(bi_gt, A_index, B_index):
    # This function calculates r^2 for two loci, A and B
    # consider two LOCI at a time 
    # in bi_gt, rows are loci, columns are genotypes

    A_gt = bi_gt[A_index,:]
    B_gt = bi_gt[B_index,:]
    
    A_gt = bi_gt[A_index, :]
    B_gt = bi_gt[B_index, :]
    
    total = len(A_gt)
    
    # Calculate allele counts
    count_11 = np.sum(A_gt & B_gt)
    count_10 = np.sum(A_gt & ~B_gt)
    count_01 = np.sum(~A_gt & B_gt)
    count_00 = np.sum(~A_gt & ~B_gt)
    
    fA = (count_11 + count_10) / total
    fB = (count_11 + count_01) / total
    fa = 1 - fA
    fb = 1 - fB
    F_AB = count_11 / total        
    
    # Calculate DAB
    DAB = F_AB - (fA * fB)
    
    # Calculate r-squared
    r_squared = (DAB**2) / (fA * fa * fB * fb)
    
    if np.isnan(r_squared):
        print("Hmmm.. D is 0? Are you sure you filtered out monomorphic loci?")
    
    return r_squared


def getrSquaredDecayDist(nodes, gt_matrix, positions, seq_length): #positions = mts.tables.sites.position
    # This function calculates the distance class at which r^2 decays to near its minimum value (within 0.2) 
    # get genotype matrix for subset of individuals
    # gt_matrix = mts.genotype_matrix(samples = nodes) # genotype_matrix retains order of nodes passed to it

    # retain only biallelic sites
    bi_gt = getBiGenoMatrix(gt_matrix)

    # downsample loci, randomly for now
    target_loci = random.sample([*range(0, np.shape(bi_gt)[0], 1)], 1000) 

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
        
    # define number of bins and bin edges
    num_bins = 10000
    bin_edges = np.linspace(0, seq_length, num_bins + 1)
    
    # use bins to add distance class column to data frame
    long_format_distance_df['distance_class'] = pd.cut(long_format_distance_df['distance'], bins=bin_edges, right=True, labels = bin_edges[0:len(bin_edges)-1]) #.apply(lambda x: x.left)
    # long_format_distance_df.head()
    
    # calculate mean r2 value for each distance class bin
    mean_r_2 = long_format_distance_df.groupby('distance_class')['r_2'].mean()
    
    # calculate mean and standard error of r^2 for last ~20% of bins 
    
    mean_r_2_df = mean_r_2.to_frame().reset_index()
    asym_mean = mean_r_2_df[mean_r_2_df['distance_class'].astype(int) > 80000000]['r_2'].mean()
    
    asym_se = mean_r_2_df[mean_r_2_df['distance_class'].astype(int) > 80000000]['r_2'].sem()
    
    # find distance class at which average r^2 goes within 2 SE of the mean r^2
    
    limit = asym_mean + (asym_se * 2)
    
    distance = mean_r_2_df[mean_r_2_df['r_2'] < limit].iloc[0]['distance_class']
    
    # return distance class, mean r^2 for last 20% and SE for last 20%    
    return(asym_mean, asym_se, limit, distance)
    
    
def bootstrap(group1_nodes, group2_nodes, g1_len, g2_len, niter):
    # This function does a series of bootstrap hypothesis tests, comparing pi, theta, and LD distance decay for the two groups of nodes given 
    
    # get list of all nodes
    node_combo = group1_nodes + group2_nodes
    node_no = len(node_combo)

    #get genotype matrix for (all) samples, accounting for the possibility of group 1 and group 2 having duplicate nodes between them during temporal sampling
    if len(node_combo) == len(np.unique(node_combo)):
        gt_matrix = mts.genotype_matrix(samples = node_combo) # group 1 will be 0-9 and group 2 will be 10-21?
    else: 
        print(f"duplicates in group 1 nodes and group 2 nodes!")
        group1_matrix = mts.genotype_matrix(samples = group1_nodes)
        group2_matrix = mts.genotype_matrix(samples = group2_nodes)
        gt_matrix = np.hstack((group1_matrix, group2_matrix))    # group 1 and 2 matrices should be: 
    
    # define indicies that correspond to different groups
    group1 = [*range(0, g1_len, 1)]
    group2 = [*range(g1_len, g1_len + g2_len, 1)]
    
    # get SFS from each sample set
    sfs1 = getSFS(gt_matrix[:,group1], group1)
    sfs2 = getSFS(gt_matrix[:,group2], group2)
    
    # calculate pi for each sample set
    pi1 = calcPiSFS(sfs1[0], sfs1[1], seq_length)
    pi2 = calcPiSFS(sfs2[0], sfs2[1], seq_length)
    
    # calculate theta for each sample set
    theta1 = mts.segregating_sites(sample_sets = list(set(group1_nodes))) / np.sum([1/i for i in np.arange(1,len(group1_nodes))])
    theta2 = mts.segregating_sites(sample_sets = list(set(group2_nodes))) / np.sum([1/i for i in np.arange(1,len(group2_nodes))])
    
    # calculate LD for each sample set
    LD1 = getrSquaredDecayDist(group1_nodes, gt_matrix[:,group1], positions, seq_length)[3]
    LD2 = getrSquaredDecayDist(group2_nodes, gt_matrix[:,group2], positions, seq_length)[3]
    
    # find differences
    obs_dif_pi = np.float64(pi2 - pi1)
    obs_dif_theta = theta2 - theta1
    obs_dif_LD = LD2 - LD1

    # generate bootstrapped differences
    pi_strap_difs = []
    theta_strap_difs = []
    LD_strap_difs = []
    
    for i in range(niter):     #repeat for no_straps
    
        # resample from all nodes
        group1_rand_samp = np.random.choice(range(0, node_no), size = g1_len, replace = True)
        group2_rand_samp = np.random.choice(range(0, node_no), size = g2_len, replace = True)
                                
        # calculate stats from artificial groups and find difference
    
        # pi
        # get SFS from each sample set
        rand_sfs1 = getSFS(gt_matrix[:,group1_rand_samp], group1_rand_samp)
        rand_sfs2 = getSFS(gt_matrix[:,group2_rand_samp], group2_rand_samp)
        
        # calculate pi for each sample set
        rand_pi1 = calcPiSFS(rand_sfs1[0], rand_sfs1[1], seq_length)
        rand_pi2 = calcPiSFS(rand_sfs2[0], rand_sfs2[1], seq_length)
    
        rand_dif_pi = rand_pi2 - rand_pi1
        
        # theta
        rand_theta1 = mts.segregating_sites(sample_sets = list(set([node_combo[i] for i in group1_rand_samp]))) / np.sum([1/i for i in np.arange(1,len(group1_rand_samp))])
        rand_theta2 = mts.segregating_sites(sample_sets = list(set([node_combo[i] for i in group2_rand_samp]))) / np.sum([1/i for i in np.arange(1,len(group2_rand_samp))])
        
        rand_dif_theta = rand_theta2 - rand_theta1
            
        # LD 
        rand_LD1 = getrSquaredDecayDist([node_combo[i] for i in group1_rand_samp], gt_matrix[:,group1_rand_samp], positions, seq_length)[3]
        rand_LD2 = getrSquaredDecayDist([node_combo[i] for i in group2_rand_samp], gt_matrix[:,group2_rand_samp], positions, seq_length)[3]    
        
        rand_dif_LD = rand_LD2 - rand_LD1
        
        pi_strap_difs.append(rand_dif_pi)
        theta_strap_difs.append(rand_dif_theta)
        LD_strap_difs.append(rand_dif_LD)


    # caclulate p-value
    p_val_pi = sum([i for i in pi_strap_difs >= obs_dif_pi]) / niter
    p_val_theta = sum([i for i in theta_strap_difs >= obs_dif_theta]) / niter
    p_val_LD = sum([i for i in LD_strap_difs >= obs_dif_LD]) / niter
    
    p_vals = [pi2, pi1, p_val_pi, theta2, theta1, p_val_theta, LD1, LD2, p_val_LD]
    return(p_vals)


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


df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta', 'LD']) # add eventually 'pi_ten', 'LD_ten',

# loop through time points to calculate stats using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(0, 2, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store stats values that are calculated once per time point
    
    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta', "LD"])
        
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