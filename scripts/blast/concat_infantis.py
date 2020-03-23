import pandas as pd
import re 
import math

# note that for core gene analysis, blons actually refer to core genes
# with open('core_genes_cluster_dict.txt', 'r') as file:
#     dic = file.read().replace('\n', '')

# dic = eval(dic)

# df = pd.DataFrame.from_dict({(i,j): dic[i][j] 
#                             for i in dic.keys()  
#                             for j in dic[i].keys()},
#                             orient='index').reset_index()
# df = df.iloc[:,:4]
# df.columns = ["genome_blon", "position", "accession", "sequence"]

# genomes = []
# blons = []
# for row in df["genome_blon"]:
#     genome, blon = tuple(row)
#     genomes.append(genome)
#     blons.append(blon)
# df.drop('genome_blon', axis=1, inplace=True)

# df["genome"] = genomes
# df["blon"] = blons

# starts = []
# ends = []
# for row in df["position"]:
#     if row != None:
#         first, second = tuple(row)
#         mini = min(int(first), int(second))
#         maxi = max(int(first), int(second))
#         starts.append(mini)
#         ends.append(maxi)
#     else:
#         starts.append(None)
#         ends.append(None)

# df["start"] = starts
# df["end"] = ends

# sequences = []
# for row in df["sequence"]:
#     sequence = re.sub('>.*?sequence', '', str(row))
#     sequence1 = re.sub('>.*?genome', '', str(sequence))
#     sequence2 = re.sub(' assembly.*?1', '', str(sequence1))
#     sequences.append(sequence2)
# df["sequence"] = sequences

# cols = ["genome", "blon", "position", "start", "end", "accession", "sequence"]
# df = df[cols] 

# df.to_csv("core_genes_cluster_dict.csv")
df = pd.read_csv("~/Desktop/repos/bifido/figure2/core_genes_cluster_dict.csv")
mapping = pd.read_csv("~/Desktop/repos/bifido/figure2/genome_map.csv")

# letters = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", \
#            "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", \
#            "ab", "ac", "ad", "ae", "af"]

# for genome in list(set(df['genome'])):
#     i = 0
#     with open ('cytobands/{}.txt'.format(genome),'a') as f:
#         f.write("chr\tstart\tend\tname\tgieStain\n")
#         for index, row in df.loc[df['genome'] == genome].iterrows():
#             # print(row["start"], type(row['start']))
#             if not math.isnan(row["start"]):
#                 f.write("infantis\t{}\t{}\t{}\t{}\n".format(row["start"], row["end"], row["blon"], letters[i]))
#                 i += 1    

# will need to adapt for core genes - filter out core genes that don't exist across all genomes
missing = []
for index, row in df.iterrows():
    if row['sequence'] == "None":
        missing.append(row['blon'])
missing = set(missing)
print(len(missing), len(df))
# print(missing)

# check how many genomes each gene is missing in
gene_count = {}
genomes = {}
for index, row in df.iterrows():
    if row['sequence'] == "None":
        if row['blon'] not in gene_count:
            gene_count[row['blon']] = 1
            genomes[row['blon']] = [row['genome']]
        else:
            gene_count[row['blon']] += 1
            genomes[row['blon']].append(row['genome'])

counts = pd.DataFrame(gene_count.items(), columns=['blon', 'num_genomes'])

# add the names of genomes that the genes are missing from
counts["genomes"] = genomes.values()

counts.to_csv("~/Downloads/test.csv")

# create a separate dataframe showing which how many genes each genome has
gene_count = {}
for index, row in df.iterrows():
    if row["sequence"] != "None":
        if row['genome'] not in gene_count:
            gene_count[row['genome']] = 1
        else:
            gene_count[row['genome']] += 1

counts = pd.DataFrame(gene_count.items(), columns=['genome', 'num_genes'])
# counts.to_csv("~/Downloads/test.csv")

for gene in missing:
    df = df[df.blon != gene]

print(len(df))
print(len(set(list(df['genome']))))

# # add a column mapping genome GCF... names to species inf_1 names
# map_dict = dict(zip(list(mapping.gcfs), list(mapping.species)))

# # blon genes that exist across all genomes
# # blons = ["blon_2331", "blon_2332", "blon_2334", "blon_2335", "blon_2336", "blon_2337", "blon_2338" \
# #          "blon_2339", "blon_2340", "blon_2348", "blon_2349", "blon_2354", "blon_2355"]

# # need to concatenate across blon genes for genomes
# # only certain genes that are present for all of the genomes
# for genome in set(df["genome"]):
#     species = map_dict[genome.split(".fna")[0]]
#     with open('move/{}.txt'.format(species), 'w') as f:
#         print(">{}".format(species), file=f)
#         for index, row in df.loc[df['genome'] == genome].iterrows():
#             # if row['blon'] in blons:
#             print(row['sequence'], file=f, end='')
#         print("\n", file=f, end='')

# for index, row in df.loc[df['genome'] == "GCF_000226175.1_ASM22617v2.fna"].iterrows():
#     print(row["sequence"])

    