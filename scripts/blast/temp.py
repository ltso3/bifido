import pandas as pd
import numpy as np

df = pd.read_csv("/Users/laurentso/Desktop/repos/bifido/scripts/metadata/genome_status.csv")
blongum = open("/Users/laurentso/Desktop/repos/bifido/figure2/ftps/blongum_ftps.txt", "r")
breve = open("/Users/laurentso/Desktop/repos/bifido/figure2/ftps/breve_ftps.txt", "r")
infantis = open("/Users/laurentso/Desktop/repos/bifido/figure2/ftps/infantis_ftps.txt", "r")
longum = open("/Users/laurentso/Desktop/repos/bifido/figure2/ftps/longum_ftps.txt", "r")
suis = open("/Users/laurentso/Desktop/repos/bifido/figure2/ftps/suis_ftps.txt", "r")

blongums = []
for line in blongum.readlines():
    blongums.append(line.split("/")[-1].strip()+".fna")

breves = []
for line in breve.readlines():
    breves.append(line.split("/")[-1].strip()+".fna")

infantises = []
for line in infantis.readlines():
    infantises.append(line.split("/")[-1].strip()+".fna")

longums = []
for line in longum.readlines():
    longums.append(line.split("/")[-1].strip()+".fna")

suises = []
for line in suis.readlines():
    suises.append(line.split("/")[-1].strip()+".fna")

df["species"] = np.repeat("blongum", len(df))

print(breves)

for index, row in df.iterrows():
    if row["genome"] in breves:
        df.loc[index, "species"] = "breve"
    elif row["genome"] in infantises:
        df.loc[index, "species"] = "infantis"
    elif row["genome"] in longums:
        df.loc[index, "species"] = "longum"
    elif row["genome"] in suises:
        df.loc[index, "species"] = "suis"

df.to_csv("~/Downloads/test.csv")
print(df)