import re
import pandas as pd 
import numpy as np
import seaborn as sns

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
        
queries = [query for query in matchDict.keys() for i in range(len(matchDict[query]))]
genomes = [genome[0]+"_"+genome[1]["s_start"] for query in matchDict.keys() for genome in matchDict[query]]
lengths = [genome[1]["length"] for query in matchDict.keys() for genome in matchDict[query]]

# need to map query names to blon numbers for easier visualization

df = pd.DataFrame(np.nan, index = range(0,len(queries)), columns = ["query", "genome", "length"])
df["query"] = queries
df["genome"] = genomes
df["length"] = lengths

df_pivot = df.pivot("genome", "query", "length")

# black = present, beige = no hits (missing data)
heatmap = sns.heatmap(df_pivot.isnull(), cbar=False)
heatmap.set(yticks=[])
heatmap.set(xticks=[])
fig = heatmap.get_figure()
fig.savefig("missing.png")

df_pivot.fillna(value=0, inplace=True)
heatmap = sns.heatmap(df_pivot.astype(int))
heatmap.set(yticks=[])
heatmap.set(xticks=[])
fig = heatmap.get_figure()
fig.savefig("existing.png")