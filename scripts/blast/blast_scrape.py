import re
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

cols = ["query", "genome", "identity", "alignment_length", "mismatches", "gaps", \
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

# read in blast output in a tabular format
f = open("/Users/laurentso/Desktop/repos/bifido/blast/output/blast_output_fmt", "r")

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

# queries = []
# lengthDict = {}
# lengths are not matching up to genomes and queries
# for key in filterDict.keys():
#         queries.append(key.split("!")[1])
#         if key.split("!")[0] not in lengthDict:
#                 lengthDict[key.split("!")[0]] = [filterDict[key.split("!")[0]+"!"+key.split("!")[1]]]
#         else:
#                 lengthDict[key.split("!")[0]].append(filterDict[key.split("!")[0]+"!"+key.split("!")[1]])

queries = [key.split("!")[1] for key in filterDict.keys()]
genomes = [key.split("!")[0] for key in filterDict.keys()]
lengths = [length for length in filterDict.values()]

print(len(queries))
print(len(genomes))
print(len(lengths))
print(queries[:5])
print(genomes[:5])
print(lengths[:5])

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
df_pivot.to_csv("test.csv", sep=',')

fig, ax = plt.subplots(figsize=(10,10))
# black = present, beige = no hits (missing data)
sns.set(font_scale=0.4)
heatmap = sns.heatmap(df_pivot.isnull(), cbar=False, linewidths=.5)
fig = heatmap.get_figure()
fig.savefig("missing.png")

df_pivot.fillna(value=0, inplace=True)
plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df_pivot.astype(int))
fig = heatmap.get_figure()
fig.savefig("existing1.png")

df_pivot1 = df.pivot_table(index="query", columns="genome", values="length")

df_pivot1.fillna(value=0, inplace=True)
heatmap = sns.heatmap(df_pivot1.astype(int))
fig = heatmap.get_figure()
fig.savefig("existing2.png")
