#!/usr/bin/env python3

import yaml
import pickle
import re
import argparse
import subprocess
import os

from Bio import SeqIO


# Use R language to call BASiNET to classify sequences.
# MODEL = RF
# CALL = R LANGUAGE


def classify(
    coding_seq,
    non_coding_seq,
    tmp_file_location,
    name,
    coverage_counts,
    rcrj,
    merge,
    circrnas,
    result_file_location,
    use_classify,
):

    tmp = os.path.join(tmp_file_location, 'tmp_file')
    circ = os.path.join(tmp_file_location, f'{name}.fa')
    print(circ)
    RCRJ = os.path.join(tmp_file_location, f'{rcrj}_RCRJ.fa')

    subprocess.call(
        '''awk '$4>{} {}' {} > {}'''.format(
            coverage_counts,
            '{print $0}',
            os.path.join(tmp_file_location, rcrj),
            os.path.join(tmp_file_location, f'{rcrj}.junction_filter_result'),
        ),
        shell=True,
    )

    # getfasta
    subprocess.call(
        'bedtools getfasta -s -fi {} -bed {} -split | fold -w 60 > {}'.format(
            circ,
            os.path.join(tmp_file_location, f'{rcrj}.junction_filter_result'),
            RCRJ,
        ),
        shell=True,
    )

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
    '''.format(
        coding_seq,
        non_coding_seq,
        tmp,
        RCRJ,
        tmp + '.dat',
        RCRJ,
        os.path.join(tmp_file_location, f'{rcrj}_translated_circ.fa'),
    )

    if use_classify == 'T':
        f = open(os.path.join(tmp_file_location, f'{rcrj}_rscript.r'), 'w')
        f.write(r_script)
        f.close()
        subprocess.call(
            'Rscript ' + os.path.join(tmp_file_location, f'{rcrj}_rscript.r'),
            shell=True,
        )
        print('Classify successfully!')
    else:
        subprocess.call(
            'mv {} {}'.format(
                RCRJ,
                os.path.join(tmp_file_location, f'{rcrj}_translated_circ.fa'),
            ),
            shell=True,
        )

    subprocess.call(
        '''sed -i 's/()//g' {}'''.format(
            os.path.join(tmp_file_location, f'{rcrj}_translated_circ.fa')
        ),
        shell=True,
    )
    final_trans_file = os.path.join(
        tmp_file_location, f'{rcrj}_translated_circ.fa'
    )
    seqs = SeqIO.parse(final_trans_file, 'fasta')
    id_dic = pickle.load(
        open(os.path.join(tmp_file_location, 'junction_name_dic'), 'rb')
    )
    translated_circ_id_list = []
    for seq in seqs:
        for i in id_dic:
            start = seq.id.split(':')[-1].split('-')[0]
            end = seq.id.split(':')[-1].split('-')[1]
            if int(start) < int(i) < int(end):
                translated_circ_id_list.append(id_dic[i])
                break
    total_circ = SeqIO.parse(circrnas, 'fasta')
    final_trans_circ_seq = []
    for i in total_circ:
        if i.id in translated_circ_id_list:
            final_trans_circ_seq.append(i)
            print(i)
    base_name = rcrj.replace('.RCRJ_result.csv', '')
    suffix = 'Aligned.sortedByCoord.out.bam.merge_result'
    if base_name.endswith(suffix):
        fastqid = base_name[: -len(suffix)]
    else:
        fastqid = base_name

    final_trans_circ_seq_name = os.path.join(
        result_file_location, f'{fastqid}.fa'
    )
    subprocess.call('mkdir -p {}'.format(result_file_location), shell=True)
    SeqIO.write(final_trans_circ_seq, final_trans_circ_seq_name, 'fasta')


def fgs(result_file_location):
    fa_list = [f for f in os.listdir(result_file_location) if f.endswith('.fa')]
    for i in fa_list:
        subprocess.call(
            './requiredSoft/FragGeneScan -s {} -o {} -w 0 -t complete'.format(
                os.path.join(result_file_location, i),
                os.path.join(
                    result_file_location,
                    i.split('.')[0] + '_translated_peptides.fa',
                ),
            ),
            shell=True,
        )


def main():
    parse = argparse.ArgumentParser(
        description='This script helps to clean reads and map to genome'
    )
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()

    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.safe_load(file)
    name = fileload['genome_name']
    tmp_file_location = fileload['tmp_file_location']
    coding_seq = fileload['coding_seq']
    non_coding_seq = fileload['non_coding_seq']
    raw_read = fileload['raw_reads']
    result_file_location = fileload['result_file_location']
    coverage_counts = fileload['coverage_counts']
    merge = fileload['merge']
    circrnas = fileload['circrnas']
    use_classify = fileload['classify']

    if merge == 'T':
        rcrj = "merge_result.RCRJ_result.csv"
        classify(
            coding_seq,
            non_coding_seq,
            tmp_file_location,
            name,
            coverage_counts,
            rcrj,
            merge,
            circrnas,
            result_file_location,
            use_classify,
        )
    else:
        rcrj_results = [
            f for f in os.listdir(tmp_file_location) if f.endswith('RCRJ_result.csv')
        ]
        for rcrj in rcrj_results:
            classify(
                coding_seq,
                non_coding_seq,
                tmp_file_location,
                name,
                coverage_counts,
                rcrj,
                merge,
                circrnas,
                result_file_location,
                use_classify,
            )

    fgs(result_file_location)


if __name__ == '__main__':
    main()

