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

hmos = pd.DataFrame(np.nan, index = [i.split('_S')[0] for i in match_dict.keys()], columns = blon_dict.keys()) 
for key in match_dict.keys():
    norms = [int(num)/float(len[1]) for num, len in zip(pd.Series(match_dict[key]), blon_dict.values())]
    hmos.loc[key.split('_S')[0]] = [math.log(num+0.00000000001) for num in norms]
hmos.columns = blon_names

# need to add 13 missing sequences that have no hits
missing = ["C0388_7F_1A", "C0461_4F_1A", "C0698_2F_1A", "C0828_4F_1A", "C0839_4F_1A", "C0851_4F_1A", \
            "C0936_1F_1A", "C1026_1F_1A", "C1062_3F_1A", "C1075_1F_1A", "C1089_1F_1A", "C1090_1F_1A", \
            "C1102_1F_1A", "C2018_4F_1A"]
for miss in missing:
    hmos.loc[miss] = [-25] * 31

hmos["sample"] = hmos.index
hmos['sample'].replace('-', '_', inplace = True, regex = True)
hmos = hmos[hmos['sample'].str.contains("E") == False]
hmos.set_index('sample', inplace = True)
hmos.sort_index(inplace = True)
hmos.to_csv("~/Downloads/hmo.csv")

# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(meta.astype(int), yticklabels=False)
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_log_norm_adj.png")

# mapping = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/master_fecal_sample_12.csv')


# # need to get a dictionary of sample id to subject id
# ids = list(zip(mapping["SampleID"], mapping["SubjectID"]))
# id_dict = {}
# for sample, subject in ids:
#     if "E" in sample:
#         id_dict[sample] = str(subject)+"_E"
#     elif sample.startswith("M"):
#         id_dict[sample] = str(subject)+"_m"
#     else:
#         id_dict[sample] = subject


# meta_sort = meta.copy()
# samples = [re.sub('-', '_', i.split('_S')[0]) for i in match_dict.keys()]
# samples = samples + missing
# subjects = [str(id_dict[sample]) for sample in samples]

# # create an adjusted index such that children sorted on top, then mothers sorted on bottom
# adj_index = []
# for i in subjects:
#     if len(i.split("_")) > 1:
#         adj_index.append(int(i.split("_")[0]) + 9999) # add 9999 to the adjusted index of mothers
#     else:
#         adj_index.append(int(i))
# meta_sort.index = subjects

# # assign adj_index to column, sort df by it, then drop it
# meta_sort["adj_index"] = adj_index
# meta_sort.sort_values(by=['adj_index'], inplace=True)
# kid_df = meta_sort.copy()
# meta_sort = meta_sort.drop('adj_index', 1)

# kid_df = kid_df[kid_df["adj_index"] < 9999]
# kid_df = kid_df.drop('adj_index', 1)

# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(meta_sort.astype(int), yticklabels=False) #, yticklabels=True, figsize=(100, figure_height)) #, cmap="YlGnBu")
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_sorted.png")

# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(kid_df.astype(int), yticklabels=False) #, yticklabels=True, figsize=(100, figure_height)) #, cmap="YlGnBu")
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_kids.png")

# # plot only the subjects with blon_2355
# df_2355 = meta.copy()
# df_2355["samples"] = samples
# df_2355 = df_2355[df_2355['Blon_2355'] > -12]
# lst_2355 = list(df_2355["samples"])
# df_2355 = df_2355.drop('samples', 1)

# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(df_2355.astype(int)) #, cmap="YlGnBu")
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_2355.png")

# ---------------------------------------------------------------------------------------------------------

# need to create dataframe with ids as subject ids and breastfeeding (metadata)
# interested in sample, breastFedPercent, birthType, correctedAgeDays
metadata = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/metadatawide.csv')
metadata = metadata[['sample','birthType', 'correctedAgeDays', 'breastFedPercent']]
metadata = metadata[metadata['sample'].str.contains("E") == False] # need to filter out ethanol samples

# need to also add mgx and pcr data
# interested in Presence_mgx, Presence_gel
pcrmgx = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/infantis_pcr.csv')
pcrmgx["mgx"] = pcrmgx['Presence_mgx'].str.contains("Bifidobacterium")
pcrmgx["pcr"] = pcrmgx['Presence_gel'].str.contains("Present")
pcrmgx = pcrmgx[['SampleID_PCR', 'mgx', 'pcr']]
pcrmgx.fillna(False, inplace = True)

# drop rows that are not in the pcr data --> 922 samples
metadata = metadata[metadata['sample'].str.contains("E") == False]
metadata = metadata[metadata['sample'].isin(pcrmgx["SampleID_PCR"])]
hmos = hmos[(hmos.index).isin(pcrmgx["SampleID_PCR"])]

# order both in terms of sample
metadata = metadata.sort_values(by=['sample'])
pcrmgx = pcrmgx.sort_values(by=['SampleID_PCR'])

big = hmos.copy()
big["birthType"] = metadata["birthType"]
big["correctedAgeDays"] = metadata["correctedAgeDays"]
big["breastFedPercent"] = metadata["breastFedPercent"]
big["mgx"] = pcrmgx["mgx"]
big["pcr"] = pcrmgx["pcr"]

print(len(hmos))
hmos.to_csv("~/Downloads/hmos.csv")

big.to_csv("~/Downloads/big.csv")

# big = pd.read_csv("big.csv")

# cmap = mpl.colors.ListedColormap(['w', 'b'])
# plt.subplots(figsize=(5,15))
# heatmap = sns.heatmap(meta_df, cmap=cmap)
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_bf.png")

# from matplotlib.gridspec import GridSpec
# fig = plt.figure()

# gs = GridSpec(1,4, height_ratios=[300], width_ratios=[30,30,30,1200])

# ax1 = fig.add_subplot(gs[0,0])
# ax2 = fig.add_subplot(gs[0,1])
# ax3 = fig.add_subplot(gs[0,2])
# ax4 = fig.add_subplot(gs[0,3])
# ax1.autoscale_view('tight')

# # fragment big dataframe into breastfeeding and gene data
# sorted_kids = big.drop(["bf", "age", "mgx", "pcr"], axis = 1)
# sorted_bf = big[["bf"]]
# sorted_mgx = big[["mgx"]]
# sorted_pcr = big[["pcr"]]

# sorted_kids.to_csv("sorted_kids.csv")

# cmap1 = mpl.colors.ListedColormap(['w', 'g'])
# cmap2 = mpl.colors.ListedColormap(['w', 'r'])

# plt.subplots(figsize=(20,15))
# sns.heatmap(sorted_bf, cmap=cmap, ax=ax1, cbar = False)
# sns.heatmap(sorted_mgx, cmap=cmap1, ax=ax2, cbar = False)
# sns.heatmap(sorted_pcr, cmap=cmap2, ax=ax3, cbar = False)
# sns.heatmap(sorted_kids, ax=ax4, xticklabels=True)
# sns.set(font_scale=5)

# ax1.set_ylabel('')    
# ax2.set_ylabel('')    
# ax3.set_ylabel('')  
# ax4.set_ylabel('')    

# plt.setp(ax1.get_xticklabels(), visible=False)
# plt.setp(ax1.get_yticklabels(), visible=False)

# plt.setp(ax2.get_xticklabels(), visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)

# plt.setp(ax3.get_xticklabels(), visible=False)
# plt.setp(ax3.get_yticklabels(), visible=False)

# plt.setp(ax4.get_yticklabels(), visible=False)
# ax4.xaxis.set_tick_params(labelsize=5)

# ax1.tick_params(axis='both', which='both', length=0)
# ax2.tick_params(axis='both', which='both', length=0)
# ax3.tick_params(axis='both', which='both', length=0)
# ax4.tick_params(axis='y', which='both', length=0)

# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_both.png")

# ---------------------------------------------------------------------------------------------------------

# fig = plt.figure()

# gs = GridSpec(2,2, height_ratios=[150,1], width_ratios=[30,600])

# ax1 = fig.add_subplot(gs[0,0])
# ax2 = fig.add_subplot(gs[0,1])
# ax1.autoscale_view('tight')
# ax2.autoscale_view('tight')

# meta_df["samples"] = samples
# meta_df["true"] = meta_df.samples.apply(lambda x: True if x in lst_2355 else False)

# filtered = meta_df[meta_df["true"] == True]
# filtered = filtered.drop('samples', 1)
# filtered = filtered.drop('true', 1)

# plt.subplots(figsize=(20,15))
# sns.heatmap(filtered, cmap=cmap, ax=ax1, cbar = False)
# sns.heatmap(df_2355, ax=ax2)

# plt.setp(ax1.get_xticklabels(), visible=False)
# plt.setp(ax2.get_yticklabels(), visible=False)
# ax2.xaxis.set_tick_params(labelsize=5)
# ax1.yaxis.set_tick_params(labelsize=5)

# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/bowtie_2355_both.png")
