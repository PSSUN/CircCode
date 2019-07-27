from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import random

# The path to known translated circRNA
path = '/path/to/translated_circRNA.fa'
# The number of test circRNA
number = 1000
# Output path
output = ''

circ = SeqIO.parse(path, 'fasta')

seq_list = []

for seq in circ:
    for i in range(0, len(seq)-3, 3):
        seq_list.append(str(seq.seq[i:i+3]))

circ_list = []

for i in range(number):
    tmp = random.sample(seq_list, k=100)
    seq = ''.join(tmp)
    print(seq)
    circ_list.append(SeqRecord(Seq(seq), id=str(i + 1), description='random_seq'))
    print(circ_list[i])

SeqIO.write(circ_list, output, 'fasta')
