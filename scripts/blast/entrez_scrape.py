from Bio import SeqIO
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from urllib.error import HTTPError
from os import listdir
import time

Entrez.email = "ltso@wellesley.edu"

# for every file in the list of files (each file is a sample = row)
output = [f for f in listdir("/Users/laurentso/Desktop/repos/bifido/blast_ncbi/out")]

genome_dict = {}
for filename in output:
        print(filename)
        handle = open("/Users/laurentso/Desktop/repos/bifido/blast_ncbi/out/{}".format(filename))
        records = NCBIXML.parse(handle)
        genome_dict[filename] = {}
        for record in records:
                query_id = record.query.split()[0]
                query_len = record.query_letters
                genome_dict[filename][query_id] = []
                for alignment in record.alignments:
                        hit_id = alignment.hit_id
                        hsp = alignment.hsps[0]
                        hit_pos = (hsp.sbjct_start, hsp.sbjct_end)
                        hit_len = hsp.align_length / query_len * 100
                        if hit_len > 80:
                                print(filename, query_id, hit_id)
                                genome_dict[filename][query_id].append(hit_pos)
                                genome_dict[filename][query_id].append(hit_id)
                                try:
                                        start, end = hit_pos
                                        if start < end:
                                                seq = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", strand=1, \
                                                                seq_start=min(hit_pos), seq_stop=max(hit_pos))
                                                genome_dict[filename][query_id].append(seq.read().strip().replace('\n', ''))
                                        else:
                                                seq = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", strand=1, \
                                                                seq_start=min(hit_pos), seq_stop=max(hit_pos))
                                                for record in SeqIO.parse(seq, "fasta"):
                                                        seq = Seq(str(record.seq))
                                                seq = seq.reverse_complement()
                                                genome_dict[filename][query_id].append(str(seq).strip().replace('\n', ''))
                                except HTTPError as err:
                                        time.sleep(20)
                                        start, end = hit_pos
                                        if start < end:
                                                seq = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", strand=1, \
                                                                seq_start=min(hit_pos), seq_stop=max(hit_pos))
                                                genome_dict[filename][query_id].append(seq.read().strip().replace('\n', ''))
                                        else:
                                                seq = Entrez.efetch(db="nucleotide", id=hit_id, rettype="fasta", strand=1, \
                                                                seq_start=min(hit_pos), seq_stop=max(hit_pos))
                                                for record in SeqIO.parse(seq, "fasta"):
                                                        seq = Seq(str(record.seq))
                                                seq = seq.reverse_complement()
                                                genome_dict[filename][query_id].append(str(seq).strip().replace('\n', ''))

# f = open("cluster_dict.txt","w")
# f.write(str(genome_dict))
# f.close()