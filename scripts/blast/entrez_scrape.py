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
                        try:
                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        except HTTPError as err:
                                time.sleep(20)
                                try:
                                        if err.code == 400:
                                                gis = Entrez.esearch(db="protein", term=hit_id, retmode="txt")
                                                for gi in Entrez.read(gis)['IdList']:
                                                        acc = Entrez.efetch(db="protein", id=gi, rettype="acc")
                                                        acc_str = acc.read().strip()
                                                        tmp = Entrez.efetch(db="protein", id=acc_str, rettype="fasta", retmode="text")
                                                        if len(tmp.read().split("\n", 1)[1]) > 100:
                                                                seq = Entrez.efetch(db="protein", id=acc_str, rettype="fasta", retmode="text")
                                                                hit_id = Entrez.efetch(db="protein", id=gi, rettype="acc").read().strip()
                                        else:
                                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                                except
                                        time.sleep(200)
                                        if err.code == 400:
                                                gis = Entrez.esearch(db="protein", term=hit_id, retmode="txt")
                                                for gi in Entrez.read(gis)['IdList']:
                                                        acc = Entrez.efetch(db="protein", id=gi, rettype="acc")
                                                        acc_str = acc.read().strip()
                                                        tmp = Entrez.efetch(db="protein", id=acc_str, rettype="fasta", retmode="text")
                                                        if len(tmp.read().split("\n", 1)[1]) > 100:
                                                                seq = Entrez.efetch(db="protein", id=acc_str, rettype="fasta", retmode="text")
                                                                hit_id = Entrez.efetch(db="protein", id=gi, rettype="acc").read().strip()
                                        else:
                                                seq = Entrez.efetch(db="protein", id=hit_id, rettype="fasta", retmode="text")
                        map_dict[hit_id] = query_id
                        print(seq.read(), file=f, end='')

f = open("protein_mapping.txt","w")
f.write(str(map_dict))
f.close()