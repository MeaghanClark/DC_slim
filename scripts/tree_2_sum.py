#!/usr/bin/env python

# last updates 10/17/2023

# added ARG segment comparison

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
    

def getDoubleSingleRatio(afs):
    X = afs[1]/afs[2]
    return(X)


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
    
    # proporation of differences between bootstrapped regions that are greater than 0 
    
    theta_prop = sum((theta_boot_g2 - theta_boot_g1) > 0.0)/len(theta_boot_g1) 
    pi_prop = sum((pi_boot_g2 - pi_boot_g1) > 0.0)/len(pi_boot_g1) 
    
    # compare bootstrapped distributions using t-test
    
    theta_ttest = stats.ttest_rel(theta_boot_g1, theta_boot_g2, alternative = "less") # ‘less’: the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
    pi_ttest = stats.ttest_rel(pi_boot_g1, pi_boot_g2, alternative = "less") 
    
    
    # return output
    return(real_theta_g1, real_theta_g2, theta_prop, theta_ttest[1], theta_ttest[0], real_pi_g1, real_pi_g2, pi_prop, pi_ttest[1], pi_ttest[0]) 


#def getGeneticRel(ped_inds, ts, inds_alive, ts_inds):
#    nodes = []
#    for i in ped_inds.to_numpy():
#        focal_ind = ts.individual(int(inds_alive[np.where(ts_inds==i)])) # get inidvidual id by matching pedigree id to tskit id
#        nodes.append(focal_ind.nodes.tolist())
#    pairs = [(i, j) for i in range(len(nodes)) for j in range(len(nodes)) if i != j] # range (# of groups)
#    plt.hist(mts.genetic_relatedness(nodes, indexes= pairs, mode = 'site'))
#    return(np.mean(mts.genetic_relatedness(nodes, indexes= pairs, mode = 'site')))
#    

def getAIC(log_like, K):
    AIC = (2*K) - (2*log_like)
    return(AIC)


def mod2AIC(model):
    model_like = model.log_likelihood()
    no_params = len(model.get_params())
    AIC = getAIC(model_like, no_params)
    return(AIC)


def getMomiSFS(nodes, ts, pop_name):
    afs = ts.allele_frequency_spectrum(sample_sets = [nodes], span_normalise = False, polarised = False)
    configs = [[len(nodes) - i, i] for i in range(1, len(nodes))]
    afs_dict = {}
    for config, value in zip(configs, afs[1:200]):
        key = (tuple(config),)  # Convert the list to a tuple and create a tuple containing it
        afs_dict[key] = value

    afs_dict_flt = {key: value for key, value in afs_dict.items() if value != 0}

    afs_momi = momi.site_freq_spectrum([pop_name], [afs_dict_flt], length = ts.sequence_length)    
    
    return(afs_momi)
    
    
def runMomiModels(SFS, gen_time):

    # define model –– no population size change
    mod_constant = momi.DemographicModel(N_e = 10000, gen_time = float(gen_time), muts_per_gen=float(mu))

    # add data
    mod_constant.set_data(SFS)

    # define parameters to infer
    mod_constant.add_size_param("N_c", lower = 10, upper = 1e7)

    # add samples
    mod_constant.add_leaf("pop", N="N_c")

    # infer parameters

    mod_constant.optimize(method="TNC")


    # define model –– bottleneck
    mod_bn = momi.DemographicModel(N_e = 10000, gen_time = float(gen_time), muts_per_gen=float(mu))

    # add data
    mod_bn.set_data(SFS)

    # define parameters to infer
    mod_bn.add_size_param("N_pre", lower = 10, upper = 1e7)
    mod_bn.add_size_param("N_post", lower = 10, upper = 1e7)
    mod_bn.add_size_param("T_bn", lower = 0, upper = 1e3)


    # add samples
    mod_bn.add_leaf("pop", t=0, N="N_post")
    mod_bn.set_size("pop", t="T_bn", N="N_pre")

    # infer parameters

    mod_bn.optimize(method="TNC")

    # compare models
    # "verdict" is whether a bottleneck was detected (True) or not (False) 
    if(mod2AIC(mod_bn) < mod2AIC(mod_constant)):
        if(mod_bn.get_params()['N_pre'] > mod_bn.get_params()['N_post']):
            verdict = True
        else: 
            verdict = False
    elif(mod2AIC(mod_bn) >= mod2AIC(mod_constant)):
        verdict = False
    
    list = [verdict, mod2AIC(mod_constant), mod_constant.log_likelihood(), mod_constant.get_params()['N_c'], 
              mod2AIC(mod_bn), mod_bn.log_likelihood(), mod_bn.get_params()['N_pre'], mod_bn.get_params()['N_post'], mod_bn.get_params()['T_bn']]

    return(list)

def ARGtest(group1_nodes, group2_nodes, ts):
    # new function that calculates theta and pi for each segment in the ARG instead of resampling/bootstrapping
    
    # input: group 1 nodes, group 2 nodes, number of interations
    
    ## calculate stats using all nodes
    # theta
    real_theta_g1 = mts.segregating_sites(sample_sets = group1_nodes) / np.sum([1/i for i in np.arange(1,len(group1_nodes))])
    real_theta_g2 = mts.segregating_sites(sample_sets = group2_nodes) / np.sum([1/i for i in np.arange(1,len(group2_nodes))])
    
    # pi
    real_pi_g1 = mts.diversity(sample_sets = group1_nodes)
    real_pi_g2 = mts.diversity(sample_sets = group2_nodes)
    
    ## define ARG segment breakpoints
    
    tree_coords = [i for i in ts.breakpoints()] # Thanks to A. Lewanski for this code! 
    
    segment_theta_g1 = ts.segregating_sites(sample_sets = group1_nodes, windows = tree_coords, mode = 'site', span_normalise = False) / np.sum([1/i for i in np.arange(1,len(group1_nodes))])
    segment_theta_g2 = ts.segregating_sites(sample_sets = group2_nodes, windows = tree_coords, mode = 'site', span_normalise = False) / np.sum([1/i for i in np.arange(1,len(group2_nodes))])    

    segment_pi_g1 = ts.diversity(sample_sets = group1_nodes, windows=tree_coords, mode='site', span_normalise=False)
    segment_pi_g2 = ts.diversity(sample_sets = group2_nodes, windows=tree_coords, mode='site', span_normalise=False)
    
    # plot for troubleshooting
    # plotDists(theta_boot_g1, theta_boot_g2, real_theta_g1, real_theta_g2, "theta", n)
    # plotDists(pi_boot_g1, pi_boot_g2, real_pi_g1, real_pi_g2, "pi", n)
    
    # proporation of differences between bootstrapped regions that are greater than 0 
    
    theta_prop = sum((segment_theta_g2 - segment_theta_g1) > 0.0)/len(segment_theta_g1) 
    pi_prop = sum((segment_pi_g2 - segment_pi_g1) > 0.0)/len(segment_pi_g1) 
    
    # compare bootstrapped distributions using t-test
    
    theta_ttest = stats.ttest_rel(segment_theta_g1, segment_theta_g2, alternative = "less") # ‘less’: the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
    pi_ttest = stats.ttest_rel(segment_pi_g1, segment_pi_g2, alternative = "less") 
    
    
    # return output
    return(real_theta_g1, real_theta_g2, theta_prop, theta_ttest[1], theta_ttest[0], real_pi_g1, real_pi_g2, pi_prop, pi_ttest[1], pi_ttest[0]) 


# ------------------------------------------------------------------------------------------------------------------------------------------

# uncomment these lines when running from command line
# sys.argv = ['tree_processing.py', '../troubleshooting/tree_nWF_20_2_100.trees','../troubleshooting/metaInd_nWF_20_2_100.txt', '/Users/meaghan/Desktop/DC_slim/het', 'hpcc_trouble', 1e-8, 21, 20]
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
# needs to change for different lifespans
# debug line: 

print(f"average age is {avg_age}")
print(f"the type of avg_age is {type(avg_age)}")

if avg_age == 2: 
    burn = 126300
elif avg_age == 5: 
    burn = 315900
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

no_straps = 100 # 100 for troubleshooting, 1000 for running

# seq_length = mts.sequence_length
# positions = mts.tables.sites.position

df_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta'])
df_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'N', 'pi', 'theta']) 

df_age_bin_boot = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_prop', 'theta_pval', 'theta_T', 'pi_younger', 'pi_older', 'pi_prop', 'pi_pval', 'pi_T'])
df_temporal_boot = pd.DataFrame(columns = ['timepoint', 'theta_future', 'theta_now', 'theta_prop', 'theta_pval', 'theta_T', 'pi_future', 'pi_now', 'pi_prop', 'pi_pval', 'pi_T'])
df_age_bin_test = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_prop', 'theta_pval', 'theta_T', 'pi_younger', 'pi_older', 'pi_prop', 'pi_pval', 'pi_T'])
df_temporal_test = pd.DataFrame(columns = ['timepoint', 'theta_future', 'theta_now', 'theta_prop', 'theta_pval', 'theta_T', 'pi_future', 'pi_now', 'pi_prop', 'pi_pval', 'pi_T'])

# momi2 
df_demo_params_all = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])
df_demo_params_young = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])
df_demo_params_yoy = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])

# loop through time points to calculate pi using tskit

for n in [*range(0, 24, 1)]: 
#for n in [*range(16, 17, 1)]: # ------------------------------------------------------------------------------------------------------------------------------------------

    # initialize data object to store stats values that are calculated once per time point
    
    # data object to store summary stats calculated from all nodes
    tp_summary = pd.DataFrame(columns = ['timepoint', 'pi', 'theta']) # "LD"
    
    # data object to store summary stats for age cohorts
    tp_age_cohort = pd.DataFrame(columns = ['timepoint', 'age', 'N', 'pi', 'theta']) 
    
    # data objects to store bootstrapped replicates of summary stats for age bins and temporal comparison
    tp_age_bin_test = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_prop', 'theta_pval', 'theta_T', 'pi_younger', 'pi_older', 'pi_prop', 'pi_pval', 'pi_T'])
    tp_temporal_test = pd.DataFrame(columns = ['timepoint', 'theta_future', 'theta_now', 'theta_prop', 'theta_pval', 'theta_T', 'pi_future', 'pi_now', 'pi_prop', 'pi_pval', 'pi_T']) # newBoot(future_nodes, now_nodes, niter = no_straps)
    
    tp_age_bin_boot = pd.DataFrame(columns = ['timepoint', 'theta_younger', 'theta_older', 'theta_prop', 'theta_pval', 'theta_T', 'pi_younger', 'pi_older', 'pi_prop', 'pi_pval', 'pi_T'])
    tp_temporal_boot = pd.DataFrame(columns = ['timepoint', 'theta_future', 'theta_now', 'theta_prop', 'theta_pval', 'theta_T', 'pi_future', 'pi_now', 'pi_prop', 'pi_pval', 'pi_T']) # newBoot(future_nodes, now_nodes, niter = no_straps)

    # data object to store results from demographic inference
    tp_demo_params_all = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])
    tp_demo_params_young = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])
    tp_demo_params_yoy = pd.DataFrame(columns = ['timepoint', 'verdict', 'con_AIC', 'con_llike', 'con_N_c', 'bn_AIC', 'bn_llike', 'bn_N_pre', 'bn_N_post', 'bn_T_bn'])

    # define tskit time
    tskit_time = convert_time.iloc[n][0]
    print(f"processing sampling point {n} representing tskit time {tskit_time}")
    
    # assign timepoint to output files    
    tp_summary.loc[0, 'timepoint'] = n
    tp_age_bin_boot.loc[0, 'timepoint'] = n
    tp_temporal_boot.loc[0, 'timepoint'] = n 
    tp_age_bin_test.loc[0, 'timepoint'] = n
    tp_temporal_test.loc[0, 'timepoint'] = n 
    tp_demo_params_all.loc[0, 'timepoint'] = n 
    tp_demo_params_young.loc[0, 'timepoint'] = n 
    tp_demo_params_yoy.loc[0, 'timepoint'] = n 
    

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
    all_nodes = [item for sublist in ind_nodes for item in sublist]
    # print(f"length of all_nodes is {len(all_nodes)}")

    ### Summary stats for entire sample------------------------------------------------------------------------------------------------------------------------------------------
        
    # with all_nodes
    tp_summary.loc[0, 'pi'] = mts.diversity(sample_sets = all_nodes)
    tp_summary.loc[0, 'theta'] = mts.segregating_sites(sample_sets = all_nodes) / np.sum([1/i for i in np.arange(1,len(all_nodes))])
    #gt_matrix_all_nodes = mts.genotype_matrix(samples = all_nodes)
    # tp_summary.loc[0, 'LD'] = getrSquaredDecayDist(all_nodes, gt_matrix_all_nodes, positions, seq_length)[3] 
    
    ### Age Cohorts------------------------------------------------------------------------------------------------------------------------------------------
            # what to fill in: 'age', 'pi', 'theta'
        
    # # calculate individual stats for samples of each age cohort
    # for i in [*range(0,len(meta['pedigree_id']))]:
    #     ind = list(meta['pedigree_id'])[i]
    #     focal_ind = mts.individual(int(alive[x.index(ind)]))# get inidvidual id by matching pedigree id to tskit id
    #     nodes = focal_ind.nodes.tolist()
    #     tp_age_cohort.loc[i, 'age'] = list(meta['age'])[i]
    #     tp_age_cohort.loc[i, 'pi'] = mts.diversity(sample_sets = nodes)
    #     tp_age_cohort.loc[i, 'theta'] = mts.segregating_sites(sample_sets = nodes) / np.sum([1/i for i in np.arange(1,len(nodes))])
    #     tp_age_cohort.loc[i, 'timepoint'] = n 
    #     tp_age_cohort.loc[i, 'pedigree_id'] = ind
    
    # calculate cohort stats 
    unique_ages = list(set(all_ages))
    noSegSites = mts.segregating_sites(sample_sets = all_nodes, span_normalise = False)
    for a in [*range(0, len(unique_ages), 1)]:
        cohort_nodes = [] 
        ids = meta[meta['age'] == unique_ages[a]][["pedigree_id"]]
        for i in ids.to_numpy():
            focal_ind = mts.individual(int(alive[np.where(x==i)])) # get inidvidual id by matching pedigree id to tskit id
            cohort_nodes.append(focal_ind.nodes.tolist())   # make list of nodes
        cohort_nodes = [item for sublist in cohort_nodes for item in sublist] # get rid of sub-lists to get overall pi 
        tp_age_cohort.loc[a, 'timepoint'] = n 
        tp_age_cohort.loc[a, 'age'] = unique_ages[a]
        tp_age_cohort.loc[a, 'N'] = len(cohort_nodes)
        tp_age_cohort.loc[a, 'pi'] = mts.diversity(sample_sets = cohort_nodes)
        tp_age_cohort.loc[a, 'theta'] = mts.segregating_sites(sample_sets = cohort_nodes) / np.sum([1/i for i in np.arange(1,len(cohort_nodes))])
        #cohort_SFS = mts.allele_frequency_spectrum(sample_sets = [cohort_nodes], span_normalise = False, polarised = False)
        #tp_age_cohort.loc[a, 'X_ratio'] = cohort_SFS[1]/cohort_SFS[2]
        #tp_age_cohort.loc[a, 'prop_SS'] = (mts.segregating_sites(sample_sets = cohort_nodes)/mts.segregating_sites(sample_sets = all_nodes))/np.sum([1/i for i in np.arange(1,len(cohort_nodes))])
        #tp_age_cohort.loc[a, 'prop_SS'] = cohort_SFS[1]/noSegSites

    # print(f"done with age cohort sampling for sampling point {n} representing tskit time {tskit_time}")
  
    ### Age Bins------------------------------------------------------------------------------------------------------------------------------------------
        
    # sort individuals by age
    meta_sorted = meta.sort_values(by='age', ascending=True) # sort metadata by age
    
    ## [1] upper/lower 10%
    # define nodes from upper and lower 10% of individuals, calculate theta, pi, and LD of entire subsample, then calculate theta, pi and LD of age bins
    
    lower_index_ten = int(len(meta_sorted) * 0.1)
    upper_index_ten = int(len(meta_sorted) * 0.9)
    
    # Get the lower 10% of individuals by age
    lower_10_percent = meta_sorted[:lower_index_ten]
    
    # Get the upper 10% of individuals by age
    upper_10_percent = meta_sorted[upper_index_ten:]
            
    # lower 10% ------------------------------------------------------------------------------------------------------------------------------------------
    lower_ids = lower_10_percent["pedigree_id"] 
    lower_ten_nodes = getNodes(ids = lower_ids, inds_alive = alive, ts = mts)

    # upper 10% ------------------------------------------------------------------------------------------------------------------------------------------
    upper_ids = upper_10_percent["pedigree_id"] 
    upper_ten_nodes = getNodes(ids = upper_ids, inds_alive = alive, ts = mts)


    # bootstrap ------------------------------------------------------------------------------------------------------------------------------------------    
    bin_vals = newBoot(lower_ten_nodes, upper_ten_nodes, niter = no_straps)
    ARG_test = ARGtest(lower_ten_nodes, upper_ten_nodes, ts = mts)
    
    # save output ------------------------------------------------------------------------------------------------------------------------------------------
    tp_age_bin_boot.iloc[:, 1:] = bin_vals
    tp_age_bin_test.iloc[:, 1:] = ARG_test

    # before nodes 
    # bootstrap ------------------------------------------------------------------------------------------------------------------------------------------
    now_nodes = lower_ten_nodes + upper_ten_nodes

    if 'future_nodes' in locals():         
        temp_vals = newBoot(future_nodes, now_nodes, niter = no_straps)
        temp_ARG_test = ARGtest(future_nodes, now_nodes, ts = mts)
       
    # save output ------------------------------------------------------------------------------------------------------------------------------------------
        tp_temporal_boot.iloc[:, 1:] = temp_vals
        tp_temporal_test.iloc[:, 1:] = temp_ARG_test

    # momi model selection
    SFS = getMomiSFS(now_nodes, mts, "pop")
    mod_summary = runMomiModels(SFS, gen_time)
    
    # momi inference on 20 young individuals only
    young_nodes = getNodes(ids = meta_sorted[0:20]["pedigree_id"], inds_alive = alive, ts = mts) # twenty youngest individuals = 40 nodes total
    young_SFS = getMomiSFS(young_nodes, mts, "pop")
    young_mod_summary = runMomiModels(young_SFS, gen_time)
    
    # momi inference on individuals younger than avg_age/10 only
    yoy_nodes = getNodes(ids = meta_sorted[meta_sorted["age"] <= round(avg_age/10)]["pedigree_id"], inds_alive = alive, ts = mts) # all individuals under avg_age/10 (rounded)
    yoy_SFS = getMomiSFS(young_nodes, mts, "pop")
    yoy_mod_summary = runMomiModels(young_SFS, gen_time)


    # save output to data object
    tp_demo_params_all.iloc[:, 1:] = mod_summary
    tp_demo_params_young.iloc[:, 1:] = young_mod_summary
    tp_demo_params_yoy.iloc[:,1:] = yoy_mod_summary

    # save output ------------------------------------------------------------------------------------------------------------------------------------------
    df_summary = pd.concat([df_summary, tp_summary], axis=0)
    df_age_cohort = pd.concat([df_age_cohort, tp_age_cohort], axis=0)
    #---
    df_age_bin_boot = pd.concat([df_age_bin_boot, tp_age_bin_boot], axis=0)   
    df_temporal_boot = pd.concat([df_temporal_boot, tp_temporal_boot], axis=0)  
    
    df_age_bin_test = pd.concat([df_age_bin_test, tp_age_bin_test], axis=0)   
    df_temporal_test = pd.concat([df_temporal_test, tp_temporal_test], axis=0)   
    #---
    df_demo_params_all = pd.concat([df_demo_params_all, tp_demo_params_all], axis = 0)
    df_demo_params_young = pd.concat([df_demo_params_young, tp_demo_params_young], axis = 0)
    df_demo_params_yoy = pd.concat([df_demo_params_yoy, tp_demo_params_yoy], axis = 0)


    # end of for loop
    future_nodes = lower_ten_nodes + upper_ten_nodes

df_summary.to_csv(outdir+"/"+prefix+"_summary.txt", sep=',', index=False)
df_age_cohort.to_csv(outdir+"/"+prefix+"_age_cohort.txt", sep=',', index=False)

df_age_bin_boot.to_csv(outdir+"/"+prefix+"_age_bin_boot.txt", sep=',', index=False)
df_temporal_boot.to_csv(outdir+"/"+prefix+"_temporal_boot.txt", sep=',', index=False)
df_age_bin_test.to_csv(outdir+"/"+prefix+"_age_bin_test.txt", sep=',', index=False)
df_temporal_test.to_csv(outdir+"/"+prefix+"_temporal_test.txt", sep=',', index=False)

df_demo_params_all.to_csv(outdir+"/"+prefix+"_demo_params_all.txt", sep=',', index=False)
df_demo_params_young.to_csv(outdir+"/"+prefix+"_demo_params_young.txt", sep=',', index=False)
df_demo_params_yoy.to_csv(outdir+"/"+prefix+"_demo_params_yoy.txt", sep=',', index=False)


print(f"done saving output")