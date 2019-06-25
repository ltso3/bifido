from Bio.Blast import NCBIXML
from Bio import Entrez
from urllib.error import HTTPError
import time

Entrez.email = "ltso@wellesley.edu"

handle = open("/Users/laurentso/Desktop/repos/bifido/blast/output/blast_output_infantis_protein.xml")
records = NCBIXML.parse(handle)
map_dict = {}
with open('protein_hits.faa', 'w') as f:
        for record in records:
                query_id = record.query
                for alignment in record.alignments:
                        hit_id = alignment.hit_id.split("|")[1]
                        map_dict[query_id] = hit_id
                        try:
                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        except HTTPError as err:
                                time.sleep(20)
                                if err.code == 400:
                                        gis = Entrez.esearch(db="protein", term=hit_id, retmode="txt")
                                        length = 0
                                        for gi in gis['IdList']:
                                                acc = Entrez.efetch(db="protein", id=gi, rettype="acc")
                                                tmp = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                                                if len(tmp.read().split("\n", 1)[1]) > length:
                                                        seq = Entrez.efetch(db="protein", id=acc, rettype="fasta", retmode="text")
                                else:
                                        seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        print(seq.read(), file=f, end='')

f = open("protein_mapping.txt","w")
f.write(str(map_dict))
f.close()