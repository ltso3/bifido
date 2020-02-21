import pandas as pd

# read in pcr data from sophie and nisreen
big = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/PCR Plates-B. infantis PCR.csv')

# switch samples ids to subject ids
mapping = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/scripts/metadata/master_fecal_sample_12.csv')
# need to get a dictionary of sample id to subject id
ids = list(zip(mapping["SampleID"], mapping["SubjectID"]))
id_dict = {}
for sample, subject in ids:
    if sample.startswith("M"):
        id_dict[sample] = str(subject)+"_m"
    else:
        id_dict[sample] = subject

subjects = [str(id_dict[sample]) for sample in big["SampleID_PCR"]]
big["SampleID_PCR"] = subjects

# need to drop mother samples
big = big[big["SampleID_PCR"].str.contains("_m") == False]

# segment into mgx and pcr dataframes
# first, mgx data for infantis presence
mgx = pd.DataFrame({"subject": big["SampleID_PCR"], 
                    "mgx": big["Presence_mgx"].str.contains("B")})
mgx.fillna(False, inplace=True)
mgx["subject"] = list(map(int, mgx["subject"]))
mgx.sort_values(by=['subject'], inplace=True)

# next, pcr data for weak or strong brand
pcr = pd.DataFrame({"subject": big["SampleID_PCR"], 
                    "pcr": big["Presence_gel"].str.contains("Present")})
pcr["subject"] = list(map(int, pcr["subject"]))    
pcr.sort_values(by=['subject'], inplace=True)

# save both as csv's
mgx.to_csv("mgx.csv")
pcr.to_csv("pcr.csv")