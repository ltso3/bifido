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

sequences = []
for row in df["sequence"]:
    sequences.append(re.sub('>.*?sequence', '', str(row)))
df["sequence"] = sequences

cols = ["genome", "blon", "position", "accession", "sequence"]
df = df[cols] 

df.to_csv("cluster_dict.csv")