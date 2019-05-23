#usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import random


path = '/path/to/circRNA.fa'

circ = SeqIO.parse(path, 'fasta')
seq_list = []
for seq in circ:
    seq_list.append(str(seq.seq))
seq = ''.join(seq_list)

circ_list = []

for i in range(1000):
    tmp = random.sample(seq, k=200)
    seq = ''
    seq = ''.join(tmp)
    print(seq)
    circ_list.append(SeqRecord(Seq(seq), id=str(i + 1), description='random'))
    print(circ_list[i])

SeqIO.write(circ_list, 'random_seq.fa', 'fasta')
