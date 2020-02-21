from Bio import Entrez
import pandas as pd

Entrez.email = "ltso@wellesley.edu"

coregenes = pd.read_csv('/Users/laurentso/Desktop/repos/bifido/figure2/journal.pone.0117912.s007.csv')
genes = list(coregenes["Gene ID**"])

with open('core_gene_proteins.faa', 'w') as f:
    for gene in genes:
        gi = Entrez.esearch(db="protein", term=gene, retmode="txt")
        seq = Entrez.efetch(db="protein", id=Entrez.read(gi)['IdList'][0], rettype="fasta", retmode="text")
        seq_str = seq.read().split("\n", 1)
        print(seq_str[0].replace("\n", ""), file=f, end='\n')
        print(seq_str[1].replace("\n", ""), file=f, end='\n')

f.close()