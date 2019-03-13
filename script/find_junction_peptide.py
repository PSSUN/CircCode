from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import pickle
import re
import rpy2.robjects as robjects


r_script = '''
library(BASiNET)
library(seqinr)
lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
classification(mRNA='/home/sun/sps/02_data/tair10_mRNA.fa',
               lncRNA='/home/sun/sps/02_data/tair10_lncRNA.fa', save="/home/sun/sps/At_circcode/tmp_file/tair10")
lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
a <- classification('/home/sun/sps/02_data/use_circCode_find_RCRJ/RCRJ.fa',lncRNA,load="/home/sun/sps/02_data/use_circCode_find_RCRJ/human.dat")
number <- which(a=='mRNA')
junction_seq=read.fasta('/home/sun/sps/02_data/use_circCode_find_RCRJ/RCRJ.fa')
mRNA_seq <- junction_seq[number]
name <- getName(mRNA_seq)
write.fasta(mRNA_seq,name,file.out = '/home/sun/sps/02_data/use_circCode_find_RCRJ/human_translated_circ.fa')
'''.format()
robjects.r(r_script)


junction2 = []
junction = pickle.load(open('/home/sun/PycharmProjects/circ-pro/junction', 'rb'))
seqs = SeqIO.parse('/home/sun/sps/02_data/use_circCode_find_RCRJ/RCRJ_translated_3.fa', 'fasta')

# position:
position = []
# stop_coden: {junction_1:[position_1 of '*' in sequence, position_2 of '*' in sequence...]...}
stop_coden = {}
# junction_dic: {junction_1:[start position, end position]...}
junction_dic = {}
# id_dic: {junction_1:[seq.id]...}
id_dic = {}

length_dic ={}
tmp = []
n=1
for seq in seqs:
    for i in re.finditer('\*', str(seq.seq)):
        tmp.append(i.start())
    # Link the junction to the start and end position with dic
    for jun in junction:
        # Get the start position and end position from seq.id
        start = seq.id.split(':')[-1].split('-')[0]
        end = seq.id.split(':')[-1].split('-')[1]
        if int(start) < int(jun) < int(end):
            length_dic[seq.id] = len(seq.seq)
            id_dic[jun] = str(seq.id)
            junction_dic[jun] = [int(start), int(end)]
            # Store location information for stop codons in each sequence
            stop_coden[jun] = tmp
            # Empty the temporary list
            tmp = []
            continue
    print('step1:', n)
    n+=1
print('-'*100)
# peptide_position:{junction_1:[peptide_start_position,peptide_end_position]...}
peptide_position = {}
tmp_2 = []
print(junction_dic)
print(stop_coden)
n=1
for i in junction_dic:
    p = int((int(i)-int(junction_dic[i][0])+1)/3)
    print('p:', p)
    if len(stop_coden[i]) == 0:
        peptide_position[i] = tmp_2
        tmp_2 = []
    elif len(stop_coden[i]) == 1:
        if p < stop_coden[i][0]:
            tmp_2.append('front')
            tmp_2.append(stop_coden[i][0])
            peptide_position[i] = tmp_2
            tmp_2 = []
        elif stop_coden[i][0] < p:

            tmp_2.append('back')
            tmp_2.append(stop_coden[i][0])
            peptide_position[i] = tmp_2
            tmp_2 = []
        elif stop_coden[i][0] == p:
            tmp_2.append('drop')
            peptide_position[i] = tmp_2
            tmp_2 = []

    else:
        for index, j in enumerate(stop_coden[i]):
            try:
                if stop_coden[i][index] < p < stop_coden[i][index+1]+1:
                    tmp_2.append(stop_coden[i][index])
                    tmp_2.append(stop_coden[i][index+1])
                    peptide_position[i] = tmp_2
                    tmp_2 = []
                    break
                elif p < stop_coden[i][0]:
                    tmp_2.append('front')
                    tmp_2.append(stop_coden[i][index])
                    peptide_position[i] = tmp_2
                    tmp_2 = []
                    break
                elif stop_coden[i][-1] < p:
                    print(stop_coden[i][-1])
                    tmp_2.append('back')
                    tmp_2.append(stop_coden[i][-1])
                    print(tmp_2)
                    peptide_position[i] = tmp_2
                    tmp_2 = []
                    break
                elif stop_coden[i][index] == p:
                    tmp_2.append('drop')
                    peptide_position[i] = tmp_2
                    tmp_2 = []
                    break
                else:
                    print('unknow error')
            except:
                print('Error in index')
        print('pep:', peptide_position)
        tmp_2 = []

    print('step2:', n)
    n += 1
final_dic = {}
print(peptide_position)
print('-'*100)
for i in peptide_position:
    for j in id_dic:
        if str(i) == str(j):
            final_dic[id_dic[i]] = peptide_position[i]

seq_list = []
print('step3')
seqs = SeqIO.parse('/home/sun/sps/02_data/use_circCode_find_RCRJ/RCRJ_translated_3.fa', 'fasta')
print(final_dic)
print('#'*3)

for seq in seqs:
    print(seq.id)
    if len(final_dic[seq.id]) == 0:
        seq_list.append(SeqRecord(Seq(str(seq.seq)),
                                  id=str(seq.id),
                                  description='peptides:' + str(len(seq.seq)) + '/' + str(length_dic[seq.id])))
    else:
        if final_dic[seq.id][0] == 'front':
            seq_list.append(SeqRecord(Seq(str(seq.seq[1:final_dic[seq.id][1]])),
                                      id=str(seq.id),
                                      description='peptides:'+str(len(seq.seq[1:final_dic[seq.id][1]]))+'/'+str(length_dic[seq.id])))
        elif final_dic[seq.id][0] == 'back':
            seq_list.append(SeqRecord(Seq(str(seq.seq[final_dic[seq.id][-1]+1:])),
                                      id=str(seq.id),
                                      description='peptides:'+str(len(seq.seq[final_dic[seq.id][1]+1:]))+'/'+str(length_dic[seq.id])))
        elif final_dic[seq.id][0] == 'drop':
            pass

        elif final_dic[seq.id] == []:
            seq_list.append(SeqRecord(Seq(str(seq.seq)),
                                      id=str(seq.id),
                                      description='peptides:'+str(len(seq.seq))+'/'+str(length_dic[seq.id])))
        else:
            seq_list.append(SeqRecord(Seq(str(seq.seq[final_dic[seq.id][0]+1:final_dic[seq.id][1]])),
                                      id=str(seq.id),
                                      description='peptides:' + str(len(seq.seq[final_dic[seq.id][0]+1:final_dic[seq.id][1]]))+'/'+str(length_dic[seq.id])))

print(seq_list)
SeqIO.write(seq_list, 'translate_3.pep.fa', 'fasta')