import re
import pandas as pd 
import numpy as np
from os import listdir
import seaborn as sns
import matplotlib.pyplot as plt

# read in all blast files, one for each sample
blast_output = [f for f in listdir("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/blast_broadecho_output")]

match_dict = {}
for filename in blast_output:
    match_dict[filename] = {}
    f = open('/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/blast_broadecho_output/{}'.format(filename))
    for line in f.readlines():
        if re.search("^# Query:", line): # if line with the HMO query
            query = "NC_011593.1:"+(line.split(":")[2].split()[0].strip())
        if re.search("hits found$", line): # if line with the number of hits
            num_hits = line.split()[1].strip()
            if query not in match_dict[filename]:
                match_dict[filename][query] = int(num_hits)
            else:
                match_dict[filename][query] += int(num_hits)

blons = open("/Users/laurentso/Desktop/repos/bifido/scripts/blast/blon_map.txt", "r")
queries = []
blon_names = []
for line in blons.readlines():
   blon, nc = line.split()
   queries.append(nc)
   blon_names.append(blon)

df = pd.DataFrame(np.nan, index = [i.split('_L001')[0] for i in match_dict.keys()], columns = queries)
for key in match_dict.keys():
    df.loc[key.split('_L001')[0]] = pd.Series(match_dict[key])
df.columns = blon_names

plt.subplots(figsize=(20,15))
heatmap = sns.heatmap(df.astype(int))
fig = heatmap.get_figure()
fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/existing_broadecho.png")