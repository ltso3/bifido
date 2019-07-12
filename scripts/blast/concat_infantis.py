import pandas as pd
import re 

with open('cluster_dict.txt', 'r') as file:
    dic = file.read().replace('\n', '')

dic = eval(dic)

df = pd.DataFrame.from_dict({(i,j): dic[i][j] 
                            for i in dic.keys()  
                            for j in dic[i].keys()},
                            orient='index').reset_index()
df = df.iloc[:,:4]
df.columns = ["genome_blon", "position", "accession", "sequence"]

genomes = []
blons = []
for row in df["genome_blon"]:
    genome, blon = tuple(row)
    genomes.append(genome)
    blons.append(blon)
df.drop('genome_blon', axis=1, inplace=True)

df["genome"] = genomes
df["blon"] = blons

starts = []
ends = []
for row in df["position"]:
    if row != None:
        first, second = tuple(row)
        mini = min(int(first), int(second))
        maxi = max(int(first), int(second))
        starts.append(mini)
        ends.append(maxi)
    else:
        starts.append(None)
        ends.append(None)

df["start"] = starts
df["end"] = ends

sequences = []
for row in df["sequence"]:
    sequence = re.sub('>.*?sequence', '', str(row))
    sequence1 = re.sub('>.*?genome', '', str(sequence))
    sequence2 = re.sub(' assembly.*?1', '', str(sequence1))
    sequences.append(sequence2)
df["sequence"] = sequences

cols = ["genome", "blon", "position", "start", "end", "accession", "sequence"]
df = df[cols] 

df.to_csv("cluster_dict.csv")

# blons = ["blon_2331", "blon_2332", "blon_2334", "blon_2335", "blon_2336", "blon_2337", "blon_2338" \
#          "blon_2339", "blon_2340", "blon_2348", "blon_2349", "blon_2354", "blon_2355"]

# # need to concatenate across blon genes for genomes
# # only certain genes that are present for all of the genomes
# for genome in set(df["genome"]):
#     with open('{}.txt'.format(genome), 'w') as f:
#         print(">{}".format(genome), file=f)
#         for index, row in df.loc[df['genome'] == genome].iterrows():
#             if row['blon'] in blons:
#                 print(row['sequence'], file=f, end='')
    