#!/usr/bin/env python
import sys
import msprime
import pyslim
import numpy as np
import pandas as pd

# use: python tree_processing.py input.trees prefix
# need to run in slim_3.4 environemnt (source activate slim_3.4)


# In[3]:

# uncomment these lines when running from command line
# sys.argv = ['tree_processing.py', './path/to/tree/file', './outdir', 'big_test']

treefile = sys.argv[1]
prefix = sys.argv[3]
outdir = sys.argv[2]

print(f"treefile is {treefile}")
print(f"prefix is {prefix}")

# set up logfile
#logfile = prefix + ".log"
#logfile = open(logfile, "w+")


# In[9]:

print(f"Loading tree file", sys.argv[1])

# read in treefile 
orig_ts = pyslim.load(treefile)


# recapitate tree
rts = orig_ts.recapitate(recombination_rate = 1e-8, Ne=5000) #may what to change this Ne


# In[111]:

# subset out individuals over the 1000 year time period of interest

# In[166]:
ind_subsets = {}
ind_label = {}
for x in range(121):
    ind_subsets[x] = rts.individuals_alive_at(x*10).tolist()
    y = []
    for i in ind_subsets[x]: 
        y.append(rts.individual(i).metadata['pedigree_id'])
    ind_label[x] = y
    
    file_name = (f"./labels/"+prefix+"/"+prefix+"_labs_"+str(x)+".txt")
    with open(file_name, 'w') as f:
        for item in y:
            f.write("%s\n" % item)

