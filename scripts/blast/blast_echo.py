import re
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# read in blast output in a tabular format
f = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/blast_output_fmt", "r")

cols = ["query", "genome", "identity", "alignment_length", "mismatches", "gaps", \
        "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

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

broadDict = {}
for genome_query in filterDict.keys():
        if genome_query.split("!")[0].split("_")[0]+"!"+genome_query.split("!")[1] not in broadDict:
                broadDict[genome_query.split("!")[0].split("_")[0]+"!"+genome_query.split("!")[1]] = \
                        filterDict[genome_query]

queries = [key.split("!")[1] for key in broadDict.keys()]
genomes = [key.split("!")[0] for key in broadDict.keys()]
lengths = [length for length in broadDict.values()]

df = pd.DataFrame(np.nan, index = range(0,len(queries)), columns = ["query", "genome", "length"])
df["query"] = queries
df["genome"] = genomes
df["length"] = lengths

df_pivot = df.pivot_table(index="genome", columns="query", values="length")

df_pivot.to_csv("output/df_pivot.csv", setp = '\t')

df_pivot.fillna(value=0, inplace=True)
plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df_pivot.astype(int))
fig = heatmap.get_figure()
fig.savefig("output/existing_echo.png")