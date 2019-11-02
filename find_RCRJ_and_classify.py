#!usr/bin/python3

import yaml
import pickle
import re
import argparse
import subprocess
import rpy2.robjects as robjects

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import time
class Translate(object):

    def __init__(self, file, tmp_file_location):
        self.file = file
        self.tmp_file_location = tmp_file_location

    def translate(self):
        seqs = SeqIO.parse(self.file, 'fasta')
        translated_seqs_1 = []
        translated_seqs_2 = []
        translated_seqs_3 = []
        for seq in seqs:

            seq_record_1 = SeqRecord(Seq(str(seq.translate().seq)), id=seq.id,description=seq.description)
            seq_record_2 = SeqRecord(Seq(str(seq[1:].translate().seq)), id=seq.id, description=seq.description)
            seq_record_3 = SeqRecord(Seq(str(seq[2:].translate().seq)), id=seq.id, description=seq.description)

            translated_seqs_1.append(seq_record_1)
            translated_seqs_2.append(seq_record_2)
            translated_seqs_3.append(seq_record_3)

        SeqIO.write(translated_seqs_1, '{}/RCRJ_translated_1.fa'.format(self.tmp_file_location), 'fasta')
        SeqIO.write(translated_seqs_2, '{}/RCRJ_translated_2.fa'.format(self.tmp_file_location), 'fasta')
        SeqIO.write(translated_seqs_3, '{}/RCRJ_translated_3.fa'.format(self.tmp_file_location), 'fasta')

# Use R language to call BASiNET to classify sequences.
# MODEL = RF
# CALL = R LANGUAGE

def classify(coding_seq, non_coding_seq, tmp_file_location,name,coverage_counts):
    tmp = tmp_file_location+'/'+'tmp_file'
    circ = tmp_file_location+'/'+name+'.fa'
    print(circ)
    RCRJ = tmp_file_location+'/'+'RCRJ.fa'
    
    subprocess.call('''awk '$4>{} {}' {} > {}'''.format(coverage_counts,'{print $0}',tmp_file_location+'/'+'junction_result',tmp_file_location+'/'+'junction_filter_result'), shell=True)
    #getfasta
    subprocess.call('bedtools getfasta -s -fi {} -bed {}  -split -name | fold -w 60 > {}'
    .format(circ,tmp_file_location+'/'+'junction_filter_result',RCRJ),shell=True)
    
    r_script = '''
    library(Biostrings)
    library(BASiNET)
    
    lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
    classification(mRNA='{}',
                   lncRNA='{}', save='{}')
    lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
    print('Make model successfully!')
    a <- classification('{}',lncRNA,load="{}")
    
    junction_seq=readDNAStringSet('{}')
    number <- which(a[1:length(junction_seq)]=='mRNA')
    mRNA_seq <- junction_seq[number]
    
    writeXStringSet(mRNA_seq,filepath = '{}')
    print(length(mRNA_seq))
    '''.format(coding_seq, non_coding_seq, tmp, RCRJ, tmp+'.dat', RCRJ, tmp_file_location+'/'+'translated_circ.fa')
    print(r_script)
    f=open('r_script.r','w')
    f.write(r_script)
    f.close()
    subprocess.call('Rscript r_script.r', shell=True)
    print('Classify successfully!')
    subprocess.call('''sed -i 's/()//g' {}'''.format(tmp_file_location+'/'+'translated_circ.fa'),shell=True)
    #robjects.r(r_script)tmp_file_location

# Find the longest peptide that can be translated from the three reading frames
def find_longest(tmp_file_location, raw_read, result_file_location, number):
    trans_seq = tmp_file_location+'/'+'RCRJ_translated_{}.fa'.format(number)
    junction = pickle.load(open('{}'.format(tmp_file_location+'/'+'junction'), 'rb'))
    seqs = SeqIO.parse('{}'.format(trans_seq), 'fasta')
    
    # The format are showed follow:
    # position:
    position = []

    # stop_coden: {junction_1:[position_1 of '*' in sequence, position_2 of '*' in sequence...]...}
    stop_coden = {}

    # junction_dic: {junction_1:[start position, end position]...}
    junction_dic = {}

    # id_dic: {junction_1:[seq.id]...}
    id_dic = {}
    
    #length_dic is a dic to storage the length of each circRNA.
    length_dic = {}
    
    # tmp is a list to storage the position.
    tmp = []
    # n is the counter of this step.
    n = 1
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
        #print('step1:', n)
        n += 1
    print('step1 finashed!')
    print('-' * 100)
    # peptide_position:{junction_1:[peptide_start_position,peptide_end_position]...}
    peptide_position = {}
    # same as list tmp
    tmp_2 = []
    
    #print(junction_dic)
    #print(stop_coden)
    
    n = 1
    for i in junction_dic:
        p = int((int(i) - int(junction_dic[i][0]) + 1) / 3)
        #print('p:', p)
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
                    if stop_coden[i][index] < p < stop_coden[i][index + 1] + 1:
                        tmp_2.append(stop_coden[i][index])
                        tmp_2.append(stop_coden[i][index + 1])
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
                        #print(stop_coden[i][-1])
                        tmp_2.append('back')
                        tmp_2.append(stop_coden[i][-1])
                        #print(tmp_2)
                        peptide_position[i] = tmp_2
                        tmp_2 = []
                        break
                    elif stop_coden[i][index] == p:
                        tmp_2.append('drop')
                        peptide_position[i] = tmp_2
                        tmp_2 = []
                        break
                    else:
                        pass
                except:
                    print('WARNNING')
            #print('pep:', peptide_position)
            tmp_2 = []

        #print('step2:', n)
        n += 1
    final_dic = {}
    #print(peptide_position)
    print('-' * 100)
    for i in peptide_position:
        for j in id_dic:
            if str(i) == str(j):
                final_dic[id_dic[i]] = peptide_position[i]
    # ----------------------------------------------------pickle.dump(self.junction, open('junction', 'wb'))
    # pickle.dump(final_dic, open('final_dic', 'wb'))
    # pickle.dump(length_dic, open('length_dic', 'wb'))
    #
    # final_dic = pickle.load(open('final_dic', 'rb'))
    # length_dic = pickle.load(open('length_dic', 'rb'))
    # final_genome = SeqRecord(Seq(str(self.genome)), id='1_CircularRNA', description='DoubleSeqWith50N')
    seq_list = []
    #print('step3')
    seqs = SeqIO.parse('{}'.format(trans_seq), 'fasta')
    #print(final_dic)

    for seq in seqs:
        print(seq.id)
        print(len(final_dic[seq.id]))
     
        if len(final_dic[seq.id]) == 0:
            seq_list.append(SeqRecord(Seq(str(seq.seq)),
                                      id=str(seq.id),
                                      description='peptides:' + str(len(seq.seq)) + '/' + str(length_dic[seq.id])))
        else:
            if final_dic[seq.id][0] == 'front':
                seq_list.append(SeqRecord(Seq(str(seq.seq[1:final_dic[seq.id][1]])),
                                          id=str(seq.id),
                                          description='peptides:' + str(
                                              len(seq.seq[1:final_dic[seq.id][1]])) + '/' + str(length_dic[seq.id])))
            elif final_dic[seq.id][0] == 'back':
                seq_list.append(SeqRecord(Seq(str(seq.seq[final_dic[seq.id][-1] + 1:])),
                                          id=str(seq.id),
                                          description='peptides:' + str(
                                              len(seq.seq[final_dic[seq.id][1] + 1:])) + '/' + str(length_dic[seq.id])))
            elif final_dic[seq.id][0] == 'drop':
                pass
            elif final_dic[seq.id] == []:
                seq_list.append(SeqRecord(Seq(str(seq.seq)),
                                          id=str(seq.id),
                                          description='peptides:' + str(len(seq.seq)) + '/' + str(length_dic[seq.id])))
            else:
                seq_list.append(SeqRecord(Seq(str(seq.seq[final_dic[seq.id][0] + 1:final_dic[seq.id][1]])),
                                          id=str(seq.id),
                                          description='peptides:' + str(
                                              len(seq.seq[final_dic[seq.id][0] + 1:final_dic[seq.id][1]])) + '/' + str(
                                              length_dic[seq.id])))

    #print(seq_list)


    # Output finally result:
    reads_jun = list(peptide_position.keys())
    print(len(reads_jun))
    pickle.dump(reads_jun, open('{}/reads_jun'.format(tmp_file_location),'wb'))
    SeqIO.write(seq_list, '{}'.format(tmp_file_location+'/'+'result_pep_{}.fa'.format(number)), 'fasta')


def main():
    parse = argparse.ArgumentParser(description='This script helps to clean reads and map to genome')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()

    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.load(file)
    name = fileload['genome_name']
    tmp_file_location = fileload['tmp_file_location']
    coding_seq = fileload['coding_seq']
    non_coding_seq = fileload['non_coding_seq']
    raw_read = fileload['raw_reads']
    result_file_location = fileload['result_file_location']
    coverage_counts = fileload['coverage_counts']
    classify(coding_seq, non_coding_seq, tmp_file_location,name,coverage_counts)
#    subprocess.call('''sed -i 's/()//g' {}'''.format(tmp_file_location+'/'+'translated_circ.fa'),shell=True)
    # Translate
    raw_read = raw_read[0]
    ribo_name = raw_read.split('/')[-1].split('.')[0]
    file = tmp_file_location+'/'+'translated_circ.fa'
    print(file)
    handle = Translate(file, tmp_file_location)
    handle.translate()
    
    # Find longest peptide

    for number in range(1,4):
        find_longest(tmp_file_location, ribo_name, result_file_location, number)
    

if __name__ == '__main__':
    main()
