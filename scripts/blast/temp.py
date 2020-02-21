import pandas as pd

df = pd.read_csv("output/heatmap_all.csv")
mapping = pd.read_csv("~/Desktop/repos/bifido/figure2/genome_map.csv")
map_dict = dict(zip(list(mapping.gcfs), list(mapping.species)))

# need to convert index to be blo, inf, lon
names = []
for name in df["Unnamed: 0"]:
    genome = str(name.split(".fna")[0])
    names.append(map_dict[genome])

df["Unnamed: 0"] = names
df.set_index("Unnamed: 0")

df.to_csv("output/heatmap_all_species.csv")
