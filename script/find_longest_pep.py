import argparse
import yaml


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


class Translate(object):

    def __init__(self,  file, tmp_file_location):
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

        SeqIO.write(translated_seqs_1, self.tmp_file_location+'/RCRJ_translated_1.fa', 'fasta')
        SeqIO.write(translated_seqs_2, self.tmp_file_location+'/RCRJ_translated_2.fa', 'fasta')
        SeqIO.write(translated_seqs_3, self.tmp_file_location+'/RCRJ_translated_3.fa', 'fasta')


def main():
    parse = argparse.ArgumentParser(description='This script helps to find longest peptide.')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()
    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.load(file)
    tmp_file_location = fileload['tmp_file_location']
    path1 = tmp_file_location + '/translate_1.pep.fa'
    path2 = tmp_file_location + '/translate_2.pep.fa'
    path3 = tmp_file_location + '/translate_3.pep.fa'
    tran_1 = SeqIO.parse(path1, 'fasta')
    tran_2 = SeqIO.parse(path2, 'fasta')
    tran_3 = SeqIO.parse(path3, 'fasta')

    len_dic_1 = []
    len_dic_2 = []
    len_dic_3 = []
    id_list = []

    len_dic = {}
    # Translate sequence:
    file = tmp_file_location + '/translated_circ.fa'
    handle = Translate(file, tmp_file_location)
    handle.translate()

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

    id_list = list(set(id_list))

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
    print(len(fi_list), len(se_list), len(th_list))

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
    longest_seq = longest_seq_1 + longest_seq_2 + longest_seq_3
    print('Finished!')
    SeqIO.write(longest_seq, tmp_file_location + '/longest_pep.fa', 'fasta')


if __name__ == '__main__':
    main()