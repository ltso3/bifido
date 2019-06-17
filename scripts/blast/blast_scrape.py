
import re
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cols = ["query", "genome", "identity", "alignment_length", "mismatches", "gaps", \
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

# read in blast output in a tabular format
# f = open("/Users/laurentso/Desktop/repos/bifido/blast/output/blast_output_fmt", "r")
f = open("/Users/laurentso/Desktop/repos/bifido/blast/output/blast_output_bifidum", "r")

matchDict = {}
for line in f.readlines():
    # need to reformat as if starts with NC, split the line
    # first is the query, next are the according columns
    if re.search("^NC", line):
        query, genome, identity, length, mismatches, gaps, \
        q_start, q_end, s_start, s_end, evalue, bit_score = line.split()
        if query not in matchDict:
            # inner value cannot be a dictionary bc keys are replaced, must be a list
            matchDict[query] = [[genome, {"identity": identity, "length": length, "mismatches": mismatches, \
                                         "gaps": gaps, "q_start": q_start, "q_end": q_end, "s_start": s_end, \
                                         "evalue": evalue, "bit_score": bit_score}]]
        else:
            matchDict[query].append([genome, {"identity": identity, "length": length, "mismatches": mismatches, \
                                         "gaps": gaps, "q_start": q_start, "q_end": q_end, "s_start": s_end, \
                                         "evalue": evalue, "bit_score": bit_score}])

# need to find all unique genomes and average length if more than one
filterDict = {}
for query in matchDict.keys():
        for genome in matchDict[query]:
                if genome[0]+"!"+query not in filterDict:
                        filterDict[genome[0]+"!"+query] = int(genome[1]["length"])

queries = [key.split("!")[1] for key in filterDict.keys()]
genomes = [key.split("!")[0] for key in filterDict.keys()]
lengths = [length for length in filterDict.values()]

# need to map query names to blon numbers for easier visualization
map = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/blon_map.txt", "r")
map_dict = {}
for line in map.readlines():
   blon, nc = line.split()
   map_dict[nc] = blon
queries_map = [map_dict[query] for query in queries]

df = pd.DataFrame(np.nan, index = range(0,len(queries)), columns = ["query", "genome", "length"])
df["query"] = queries_map
df["genome"] = genomes
df["length"] = lengths

df_pivot = df.pivot_table(index="genome", columns="query", values="length")
print(df_pivot)
fig, ax = plt.subplots(figsize=(10,10))
# black = present, beige = no hits (missing data)
sns.set(font_scale=0.4)
heatmap = sns.heatmap(df_pivot.isnull(), cbar=False, linewidths=.5)
fig = heatmap.get_figure()
# fig.savefig("output/missing.png")
fig.savefig("output/missing_bifidum.png")

df_pivot.fillna(value=0, inplace=True)
plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df_pivot.astype(int))
fig = heatmap.get_figure()
# fig.savefig("output/existing.png")
fig.savefig("output/existing_bifidum.png")