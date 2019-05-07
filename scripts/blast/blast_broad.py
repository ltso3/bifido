import re
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cols = ["query", "genome", "identity", "alignment_length", "mismatches", "gaps", \
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

# read in blast output in a tabular format
f = open("/Users/laurentso/Desktop/repos/bifido/blast_copy/output/blast_output_fmt", "r")

matchDict = {}
for line in f.readlines():
    # need to reformat as if starts with NC, split the line
    # first is the query, next are the according columns
    if re.search("^NC", line):
        query, genome, identity, length, mismatches, gaps, \
        q_start, q_end, s_start, s_end, evalue, bit_score = line.split()
        if query not in matchDict:
            # inner value cannot be a dictionary bc keys are replaced, must be a list
            matchDict[query] = [genome.split("!")[0]]
        else:
            matchDict[query].append(genome.split("!")[0])

# read in all genome names
g = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/genomes.txt", "r")
all_genomes = []
for line in g.readlines():
    all_genomes.append(line.strip())

# need to map query names to blon numbers for easier visualization
queries = set(matchDict.keys())

map = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/blon_map.txt", "r")
map_dict = {}
for line in map.readlines():
   blon, nc = line.split()
   map_dict[nc] = blon
queries_map = [map_dict[query] for query in queries]

# need to set genomes as index and queries as columns
df = pd.DataFrame(1, index = range(0,len(all_genomes)), columns = queries_map)

# set genomes as the index
df["genomes"] = all_genomes
df.set_index("genomes", inplace=True, drop=True)

# fill in dataframe such that genomes with a hit are marked with a 1, without with a 0
for query in matchDict.keys():
    map_query = map_dict[query]
    for genome in matchDict[query]:
        df.loc[str(genome), str(map_query)] = 0

fig, ax = plt.subplots(figsize=(10,10))

df.fillna(value=0, inplace=True)
plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df.astype(int))
fig = heatmap.get_figure()
fig.savefig("output/existing_broad.png")