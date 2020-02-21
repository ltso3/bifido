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

# need to add 13 child and 20 mother missing sequences that have no hits
missing = ["C0388_7F_1A", "C0461_4F_1A", "C0698_2F_1A", "C0828_4F_1A", "C0839_4F_1A", "C0851_4F_1A", \
            "C0936_1F_1A", "C1026_1F_1A", "C1062_3F_1A", "C1075_1F_1A", "C1089_1F_1A", "C1090_1F_1A", \
            "C1102_1F_1A", "C2018_4F_1A", "M1008_3F_1A", "M1022_2F_1A", "M1056_2F_1A", "M1057_1F_1A", \
            "M1059_1F_1A", "M1061_2F_1A", "M1064_2F_1A", "M1067_1F_1A", "M1078_1F_1A", \
            "M1082_1F_1A", "M1083_1F_1A", "M1088_1F_1A", "M1095_1F_1A", "M1096_1F_1A", "M1099_1F_1A", \
            "M1109_1F_1A", "M1110_1F_1A", "M1113_1F_1A", "M1115_1F_1A", "M1116_1F_1A"]
for miss in missing:
    hmos.loc[miss] = [-25] * 31

hmos["sample"] = hmos.index
hmos['sample'].replace('-', '_', inplace = True, regex = True)
hmos = hmos[hmos['sample'].str.contains("E") == False]
hmos.set_index('sample', inplace = True)
hmos.sort_index(inplace = True)

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
big["birthType"] = list(metadata["birthType"])
big["correctedAgeDays"] = list(metadata["correctedAgeDays"])
big["breastFedPercent"] = list(metadata["breastFedPercent"])
big["mgx"] = list(pcrmgx["mgx"])
big["pcr"] = list(pcrmgx["pcr"])

hmos.to_csv("hmos.csv")
big.to_csv("big.csv")

# ---------------------------------------------------------------------------------------------------------

# drop mothers except for the four found by PCR
mothers = ["M0886_3F_1A", "M0951_2F_1A", "M1183_1F_1A", "M1153_1F_1A"]
big["sample"] = big.index
big = big[(big["sample"].str.contains("C")) | (big.index.isin(mothers))]
big = big.drop("sample", axis = 1)

# sort by breastfeeding percentage
big = big.sort_values(by=["breastFedPercent"])

# only plot samples found to harbor infantis by HMOs or PCR
# small_big = big[(big["pcr"] == True) | (big.index.isin(["C1217_1F_1A", "C0611_1F_1A"]))]

from matplotlib.gridspec import GridSpec
fig = plt.figure(figsize = (300, 200))

gs = GridSpec(1,3, height_ratios=[300], width_ratios=[30,30,1200])

ax1 = fig.add_subplot(gs[0,0])
ax2 = fig.add_subplot(gs[0,1])
ax3 = fig.add_subplot(gs[0,2])
ax1.autoscale_view('tight')
ax2.autoscale_view('tight')

# fragment big dataframe into breastfeeding and gene data
hmos = big.drop(["breastFedPercent", "correctedAgeDays", "mgx", "pcr", "birthType"], axis = 1)
sorted_bf = big[["breastFedPercent"]]
sorted_mgx = big[["mgx"]]
sorted_pcr = big[["pcr"]]

cmap1 = mpl.colors.ListedColormap(['w', 'b'])
cmap2 = mpl.colors.ListedColormap(['w', 'r'])

sns.heatmap(sorted_mgx, cmap=cmap1, ax=ax1, cbar = False)
sns.heatmap(sorted_pcr, cmap=cmap2, ax=ax2, cbar = False)
sns.heatmap(hmos, ax=ax3, xticklabels=True)
sns.set(font_scale=5)

ax1.set_ylabel('')  
ax1.set_xlabel('mgx')

ax2.set_ylabel('')    
ax2.set_xlabel('pcr')

ax3.set_ylabel('')  

plt.setp(ax1.get_yticklabels(), visible=False)

plt.setp(ax2.get_yticklabels(), visible=False)

plt.setp(ax3.get_yticklabels(), visible=False)
ax3.xaxis.set_tick_params(labelsize=5)

ax1.tick_params(axis='both', which='both', length=0)
ax2.tick_params(axis='both', which='both', length=0)
ax3.tick_params(axis='y', which='both', length=0)

fig.savefig("/Users/laurentso/Desktop/repos/bifido/scripts/blast/output/figure4_sorted.png")


