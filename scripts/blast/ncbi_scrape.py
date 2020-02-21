import re
import math
import pandas as pd 
import numpy as np
from os import listdir
import seaborn as sns
import matplotlib.pyplot as plt

# read in all blast files, one for each sample
blast_output = [f for f in listdir("/Users/laurentso/Desktop/repos/bifido/figure2/ncbi_out")]

match_dict = {}
for filename in blast_output:
    match_dict[filename] = {}
    f = open('/Users/laurentso/Desktop/repos/bifido/figure2/ncbi_out/{}'.format(filename))
    for line in f.readlines():
        if re.search("^# Query:", line): # if line with the HMO query
            query = "NC_011593.1:"+(line.split(":")[2].split()[0].strip())

        if re.search("hits found$", line): # if line with the number of hits
            num_hits = line.split()[1].strip()
            if query not in match_dict[filename]:
                match_dict[filename][query] = int(num_hits)
            else:
                match_dict[filename][query] += int(num_hits)

        # if alignment length is less than 90bp, consider number of hits to be 0
        if re.search("^NC", line):
            query, genome, identity, length, mismatches, gaps, \
            q_start, q_end, s_start, s_end, evalue, bit_score = line.split()
            if int(length) < 90:
                match_dict[filename][query] -= 1

blons = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/blon_map.txt", "r")
queries = []
blon_names = []
for line in blons.readlines():
   blon, nc = line.split()
   queries.append(nc)
   blon_names.append(blon)

lengths = open("/Users/laurentso/Desktop/repos/bifido/blast_broad/query/hmo_genes_lengths.txt", "r")
len_dict = {}
for line in lengths.readlines():
    if re.search("^NC", line):
        blon = line.split()[0]
    else:
        len_dict[blon] = line

df = pd.DataFrame(np.nan, index = [i for i in match_dict.keys()], columns = queries)
for key in match_dict.keys():
    df.loc[key] = pd.Series(match_dict[key])
df.columns = blon_names

normalized_df = df.copy()
for key in match_dict.keys():
    norms = [int(num)/int(len) for num, len in zip(pd.Series(match_dict[key]), len_dict.values())]
    normalized_df.loc[key] = [math.log(num+0.00000000001) for num in norms]

normalized_df.to_csv("output/heatmap_all.csv", sep=',')

# # normal results
# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(df.astype(int))#, vmin=0.0, vmax = 2)
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/blongum.png")

# # logged and normalized results
# plt.subplots(figsize=(20,15))
# heatmap = sns.heatmap(normalized_df.astype(int)) #, cmap="YlGnBu")
# fig = heatmap.get_figure()
# fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/blongum_log_norm.png")


