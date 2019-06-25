from Bio.Blast import NCBIXML
from Bio import Entrez
from urllib.error import HTTPError
import time

Entrez.email = "ltso@wellesley.edu"

handle = open("/Users/laurentso/Desktop/repos/bifido/blast/output/blast_output_infantis_protein.xml")
records = NCBIXML.parse(handle)
with open('protein_hits.faa', 'w') as f:
        for record in records:
                query_id = record.query
                for alignment in record.alignments:
                        hit_id = alignment.hit_id.split("|")[1]
                        try:
                                time.sleep(20)
                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        except HTTPError:
                                time.sleep(20)
                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        print(seq.read(), file=f, end='')