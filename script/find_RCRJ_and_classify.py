import yaml
import pickle
import re
import argparse
import subprocess
import rpy2.robjects as robjects

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class Translate(object):

    def __init__(self,  file):
        self.file = file

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

        SeqIO.write(translated_seqs_1, '/home/sun/sps/At_circcode/tmp_file/RCRJ_translated_1.fa', 'fasta')
        SeqIO.write(translated_seqs_2, '/home/sun/sps/At_circcode/tmp_file/RCRJ_translated_2.fa', 'fasta')
        SeqIO.write(translated_seqs_3, '/home/sun/sps/At_circcode/tmp_file/RCRJ_translated_3.fa', 'fasta')

#Use R language to call BASiNET to classify sequences.
def classify(coding_seq, non_coding_seq, tmp_file_location, raw_read, name):
    tmp = tmp_file_location+'/'+'tmp_file'
    raw_read = raw_read[0]
    ribo_name = raw_read.split('/')[-1].split('.')[-2]
    circ = tmp_file_location+'/'+name+'.fa'

    #getfasta
    subprocess.call('bedtools getfasta -s -fi {} -bed {}  -split -name | fold -w 60 > {}'
    .format(tmp_file_location,tmp_file_location+ribo_name+'_junction_result',tmp_file_location+'RCRJ.fa'),shell=True)
    r_script = '''
    library(seqinr)
    library(BASiNET)
    
    lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
    classification(mRNA='{}',
                   lncRNA='{}', save='{}')
    lncRNA <- system.file("extdata", "sequences.fasta", package = "BASiNET")
    
    a <- classification('{}',lncRNA,load="{}")
    number <- which(a=='mRNA')
    junction_seq=read.fasta('{}')
    mRNA_seq <- junction_seq[number]
    name <- getName(mRNA_seq)
    write.fasta(mRNA_seq,name,file.out = '{}')
    '''.format(coding_seq, non_coding_seq, tmp, circ, tmp, circ, tmp_file_location+raw_read+'_translated_circ.fa')
    print(r_script)
    f=open('r_script.r','w')
    f.write(r_script)
    f.close()
    subprocess.call('Rscript r_script.r', shell=True)
    #robjects.r(r_script)

# Find the longest peptide that can be translated from the three reading frames
def find_longest(tmp_file_location, raw_read, result_file_location):
    
    trans_seq = tmp_file_location+'/'+raw_read+'_translated_circ.fa'
    junction = pickle.load(open('{}'.format(tmp_file_location+'/'+'junction'), 'rb'))
    seqs = SeqIO.parse('{}'.format(trans_seq), 'fasta')
    
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
    n = 1
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
    seqs = SeqIO.parse('{}'.format(trans_seq), 'fasta')
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
    
    # print(seq_list)
    
    # Output finally result:
    SeqIO.write(seq_list, '{}'.format(result_file_location+'/'+'result_pep.fa'), 'fasta')


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
    classify(coding_seq, non_coding_seq, tmp_file_location, raw_read, name)
    
    # Translate
    file = tmp_file_location+raw_read[0]+'_translated_circ.fa'
    handle = Translate(file)
    handle.translate()
    
    # Find longest peptide
    find_longest(tmp_file_location, raw_read, result_file_location)
    

if __name__ == '__main__':
    main()
