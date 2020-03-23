import pandas as pd

df = pd.read_csv("~/Downloads/genome_status.csv")

split = list(df["genome	representation	level	coverage"])
genome = [i.split("\t", 1)[0] for i in split]
representation = [i.split("\t", 2)[1] for i in split]
level = [i.split("\t", 3)[2] for i in split]
coverage = [i.split("\t", 4)[3] for i in split]

new = pd.DataFrame(columns = ['genome', 'representation', 'level', 'coverage']) 
new['genome'] = genome
new['representation'] = representation
new['level'] = level
new['coverage'] = coverage

new.to_csv("../metadata/genome_status.csv")