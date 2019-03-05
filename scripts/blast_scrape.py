import re
import pandas as pd 
import numpy as np

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
            matchDict[query] = {genome: {"identity": identity, "length": length, "mismatches": mismatches, \
                                         "gaps": gaps, "q_start": q_start, "q_end": q_end, "s_start": s_end, \
                                         "evalue": evalue, "bit_score": bit_score}}
        else:
            matchDict[query][genome] = {"identity": identity, "length": length, "mismatches": mismatches, \
                                        "gaps": gaps, "q_start": q_start, "q_end": q_end, "s_start": s_end, \
                                         "evalue": evalue, "bit_score": bit_score}

# dictionary of (hmo_gene, list of (species, e value)
# plot as a heatmap
# pull out query, genome, evalue
queries = [query for query in matchDict.keys() for i in range(len(matchDict[query]))]
genomes = [genome for genome in matchDict[query] for query in matchDict.keys()]

df = pd.DataFrame(np.nan, index = range(0,len(queries)), columns = cols)

print(len(queries))
print(len(genomes))
print(len(matchDict["NC_011593.1:c2617928-2616333"]))
print(matchDict["NC_011593.1:c2617928-2616333"])
# evals = [genome[1]["evalue"] for genome in matchDict.values()]

listy = ['NC_011593.1_cds_WP_012578586.1_2437', {'identity': '100.000', 'length': '1119', 'mismatches': '0', 'gaps': '0', 'q_start': '1', 'q_end': '1119', 's_start': '1119', 'evalue': '0.0', 'bit_score': '2067'}]
# print(listy[0])
# print(len(genomes))
# print(queries[:10])