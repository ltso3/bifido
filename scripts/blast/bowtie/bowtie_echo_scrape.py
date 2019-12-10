import pandas as pd 
import re
import math
import numpy as np
from os import listdir
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

# read in all bowtie files, one for each sample
bowtie_output = [f for f in listdir("/Users/laurentso/Desktop/repos/bifido/scripts/blast/bowtie/counts/")]

# compute blon mapping and lengths here
blons = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/blon_map.txt", "r")
blon_dict = {}
blon_names = []
for line in blons.readlines():
   blon, nc = line.split()
   blon_dict[nc.split(":")[1]] = [blon]
   blon_names.append(blon)

lengths = open("/Users/laurentso/Desktop/repos/bifido/blast_broad/query/hmo_genes_lengths.txt", "r")
for line in lengths.readlines():
    if re.search("^NC", line):
       gene = line.split()[0].split(":")[1]
    else:
        blon_dict[gene].append(line)

match_dict = {}
for filename in bowtie_output:
    match_dict[filename] = {}
    df = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/blast/bowtie/counts/{}'.format(filename), sep='\t')
    
    for gene in df["ContigName"]:
        hmo = blon_dict[gene][0]
        if hmo not in match_dict[filename]:
            match_dict[filename][hmo] = int(df["totalreadsgoodmapped"][df["ContigName"] == gene])
        else:
            match_dict[filename][hmo] += int(df["totalreadsgoodmapped"][df["ContigName"] == gene])

# ---------------------------------------------------------------------------------------------------------

meta = pd.DataFrame(np.nan, index = [i.split('_S')[0] for i in match_dict.keys()], columns = blon_dict.keys()) 
for key in match_dict.keys():
    norms = [int(num)/float(len[1]) for num, len in zip(pd.Series(match_dict[key]), blon_dict.values())]
    meta.loc[key.split('_S')[0]] = [math.log(num+0.00000000001) for num in norms]
meta.columns = blon_names

# print(meta)
# print(match_dict.keys()) // M0769-2E-1A_S1_counts.txt

plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(meta.astype(int), yticklabels=False)
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_log_norm_adj.png")

mapping = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/master_fecal_sample_12.csv')

# need to get a dictionary of sample id to subject id
ids = list(zip(mapping["SampleID"], mapping["SubjectID"]))
id_dict = {}
for sample, subject in ids:
    if sample.startswith("M"):
        id_dict[sample] = str(subject)+"_m"
    else:
        id_dict[sample] = subject

meta_sort = meta.copy()
samples = [re.sub('-', '_', i.split('_S')[0]) for i in match_dict.keys()]
subjects = [str(id_dict[sample]) for sample in samples]

# create an adjusted index such that children sorted on top, then mothers sorted on bottom
adj_index = []
for i in subjects:
    if len(i.split("_")) > 1:
        adj_index.append(int(i.split("_")[0]) + 9999) # add 9999 to the adjusted index of mothers
    else:
        adj_index.append(int(i))

meta_sort.index = subjects

# assign adj_index to column, sort df by it, then drop it
meta_sort["adj_index"] = adj_index
meta_sort.sort_values(by=['adj_index'], inplace=True)
kid_df = meta_sort.copy()
meta_sort = meta_sort.drop('adj_index', 1)

kid_df = kid_df[kid_df["adj_index"] < 9999]
kid_df = kid_df.drop('adj_index', 1)

plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(meta_sort.astype(int), yticklabels=False) #, yticklabels=True, figsize=(100, figure_height)) #, cmap="YlGnBu")
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_sorted.png")

plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(kid_df.astype(int), yticklabels=False) #, yticklabels=True, figsize=(100, figure_height)) #, cmap="YlGnBu")
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_kids.png")

# plot only the subjects with blon_2355
df_2355 = meta.copy()
df_2355["samples"] = samples
df_2355 = df_2355[df_2355['Blon_2355'] > -12]
lst_2355 = list(df_2355["samples"])
df_2355 = df_2355.drop('samples', 1)

plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df_2355.astype(int)) #, cmap="YlGnBu")
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_2355.png")

# ---------------------------------------------------------------------------------------------------------

# need to create dataframe with ids as subject ids and breastfeeding (metadata)
metadata = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/bf.csv')

children = [i for i in match_dict.keys() if not i.startswith("M")]
samples = [re.sub('-', '_', i.split('_S')[0]) for i in children]
subjects = [id_dict[sample] for sample in samples]
meta_df = pd.DataFrame(np.nan, index = subjects, columns = ["bf"])
meta_df.sort_index(inplace=True)

for key in subjects:
    if key in list(metadata["subject"]):
        index = metadata.index[metadata['subject'] == key]
        # convert T/F to 1/0
        if(metadata.iloc[index[0]]['bf']):
            meta_df.loc[key] = 1
        else:
            meta_df.loc[key] = 0    
    else:
        meta_df.loc[key] = np.nan

cmap = mpl.colors.ListedColormap(['w', 'b'])
plt.subplots(figsize=(5,15))
heatmap = sns.heatmap(meta_df, cmap=cmap)
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_bf.png")

from matplotlib.gridspec import GridSpec
fig = plt.figure()

gs = GridSpec(2,2, height_ratios=[150,1], width_ratios=[30,600])

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax1.autoscale_view('tight')

plt.subplots(figsize=(20,15))
sns.heatmap(meta_df, cmap=cmap, ax=ax1, cbar = False)
sns.heatmap(kid_df, ax=ax2)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.xaxis.set_tick_params(labelsize=5)
ax1.yaxis.set_tick_params(labelsize=5)

fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_both.png")

# ---------------------------------------------------------------------------------------------------------

fig = plt.figure()

gs = GridSpec(2,2, height_ratios=[150,1], width_ratios=[30,600])

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax1.autoscale_view('tight')
ax2.autoscale_view('tight')

meta_df["samples"] = samples
meta_df["true"] = meta_df.samples.apply(lambda x: True if x in lst_2355 else False)

filtered = meta_df[meta_df["true"] == True]
filtered = filtered.drop('samples', 1)
filtered = filtered.drop('true', 1)

plt.subplots(figsize=(20,15))
sns.heatmap(filtered, cmap=cmap, ax=ax1, cbar = False)
sns.heatmap(df_2355, ax=ax2)

plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax2.get_yticklabels(), visible=False)
ax2.xaxis.set_tick_params(labelsize=5)
ax1.yaxis.set_tick_params(labelsize=5)

fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_2355_both.png")
