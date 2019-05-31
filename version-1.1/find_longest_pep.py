#/usr/bin/python3

import argparse
import yaml
import time

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def main():
    parse = argparse.ArgumentParser(description='This script helps to find longest peptide.')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()
    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.load(file)
    tmp_file_location = fileload['tmp_file_location']
    raw_read = fileload['raw_reads']
    raw_read = raw_read[0].split('.')[0]
    ribo_name = raw_read.split('/')[-1].split('.')[0]

    path1 = tmp_file_location + '/result_pep_1.fa'
    path2 = tmp_file_location + '/result_pep_2.fa'
    path3 = tmp_file_location + '/result_pep_3.fa'
    
    tran_1 = SeqIO.parse(path1, 'fasta')
    tran_2 = SeqIO.parse(path2, 'fasta')
    tran_3 = SeqIO.parse(path3, 'fasta')

    len_dic_1 = []
    len_dic_2 = []
    len_dic_3 = []
    id_list = []

    len_dic = {}

    for i in tran_1:
        tmp = {'number': 1, 'id': i.id, 'length': len(i.seq)}
        len_dic_1.append(tmp)
        id_list.append(i.id)
    for i in tran_2:
        id_list.append(i.id)
        tmp = {'number': 2, 'id': i.id, 'length': len(i.seq)}
        len_dic_1.append(tmp)
    for i in tran_3:
        id_list.append(i.id)
        tmp = {'number': 3, 'id': i.id, 'length': len(i.seq)}
        len_dic_1.append(tmp)
#    print(tran_1)
#    for i in len_dic_1:
#        print(i)


    tran_1 = SeqIO.parse(path1, 'fasta')
    tran_2 = SeqIO.parse(path2, 'fasta')
    tran_3 = SeqIO.parse(path3, 'fasta')

    for pep in tran_1:
        len_dic[pep.id] = []
        len_dic[pep.id].append({1: len(pep)})
    for pep in tran_2:
        if pep.id in len_dic.keys():
            len_dic[pep.id].append({2: len(pep)})
        else:
            len_dic[pep.id] = []
            len_dic[pep.id].append({2: len(pep)})
    for pep in tran_3:
        if pep.id in len_dic.keys():
            len_dic[pep.id].append({3: len(pep)})
        else:
            len_dic[pep.id] = []
            len_dic[pep.id].append({3: len(pep)})

    tmp_list = len_dic_1 + len_dic_2 + len_dic_3

    tmp = []
    final = []
    for j in id_list:
        for i in tmp_list:
            if i['id'] == j:
                tmp.append(i)
        sort = sorted(tmp, key=lambda x: x['length'])
        final.append(sort[-1])
        tmp = []

    fi_list = []
    se_list = []
    th_list = []
    for i in final:
        if i['number'] == 1:
            fi_list.append(i)
        elif i['number'] == 2:
            se_list.append(i)
        else:
            th_list.append(i)
#    print(len(fi_list), len(se_list), len(th_list))

    fi_id = [x['id'] for x in fi_list]
    se_id = [x['id'] for x in se_list]
    th_id = [x['id'] for x in th_list]

    tran_1 = SeqIO.parse(path1, 'fasta')
    tran_2 = SeqIO.parse(path2, 'fasta')
    tran_3 = SeqIO.parse(path3, 'fasta')

    longest_seq_1 = []
    longest_seq_2 = []
    longest_seq_3 = []
    for seq in tran_1:
        if seq.id in fi_id:
            longest_seq_1.append(seq)
    for seq in tran_2:
        if seq.id in se_id:
            longest_seq_2.append(seq)
    for seq in tran_3:
        if seq.id in th_id:
            longest_seq_3.append(seq)
#    print(longest_seq_1)
    longest_seq = longest_seq_1 + longest_seq_2 + longest_seq_3
    print('Finished!')
#    print(longest_seq)
    SeqIO.write(longest_seq, tmp_file_location + '/longest_pep.fa', 'fasta')


if __name__ == '__main__':
    main()
