using Microbiome
using DataFrames
using FileIO
using CSVFiles

# load available taxonomic profiles (2 is a copy of 1?)
batch1 = FileIO.load("tax_profiles/batch001_metaphlan2_taxonomic_profiles.tsv") |> DataFrame
# batch2 = FileIO.load("tax_profiles/batch002_metaphlan2_taxonomic_profiles.tsv") |> DataFrame
batch3 = FileIO.load("tax_profiles/batch003_metaphlan2_taxonomic_profiles.tsv") |> DataFrame
batch4 = FileIO.load("tax_profiles/batch004_metaphlan2_taxonomic_profiles.tsv") |> DataFrame

# identify species counts
species1 = taxfilter(batch1, 7)
# species2 = taxfilter(batch2, 7)
species3 = taxfilter(batch3, 7)
species4 = taxfilter(batch4, 7)

FileIO.save("tax_profiles/species1.csv", species1)
FileIO.save("tax_profiles/species3.csv", species3)
FileIO.save("tax_profiles/species4.csv", species4)

