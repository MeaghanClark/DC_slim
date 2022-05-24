

alive = sts.individuals_alive_at(0) # modern age

groups = {
    'alive' : np.random.choice(alive, size=100, replace=False), # could add more groups to this
}

group_order = ['alive'] # add in groups
sampled_nodes = [[] for _ in groups]
for j, k in enumerate(group_order):
   for ind in groups[k]:
      sampled_nodes[j].extend(ts.individual(ind).nodes)


#group divergence - change range() to the number of groups compared, gives pairwise distance between group, diagonal is distance within group
pairs = [(i, j) for i in range(5) for j in range(5)] # range (# of groups)
group_div = ts.divergence(sampled_nodes, indexes=pairs).reshape((5, 5))

#ind divergence
ind_nodes = []
ind_group = []
ind_ids = []
for j, group in enumerate(group_order):
   for ind in groups[group]:
      ind_ids.append(ind)
      ind_nodes.append(ts.individual(ind).nodes)
      ind_group.append(group_order[j])

nind = len(ind_ids)
pairs = [(i, j) for i in range(nind) for j in range(nind)]
ind_div = ts.divergence(ind_nodes, indexes=pairs) # this is what gives pi, diagonal is mean pi within population,
# individual level divergence (heterozygosity) 

# save output
x = []
for i in ind_ids:
   ind = mts.individual(i)
   label = f"tsk_{ind.id}"
   x.append(label)

b = np.reshape(ind_div, (92,92))
panda_df = pd.DataFrame(data = b, columns = x)
panda_df.to_csv("../full_output/"+prefix+"pi.txt", sep=' ', index=True)

