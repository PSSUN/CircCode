#/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import subprocess
import pickle
import argparse
import yaml
import os


class Genome(object):

    def __init__(self, circ_rnas, tmp_file_location):
        self.circ_rnas = circ_rnas
        # store pieces of genome sequence in a list to avoid
        # expensive string concatenation inside the loop
        self.genome_parts = []
        self.newGenome = []
        self.circ_name = []
        self.start = []
        self.end = []
        self.increase_length = 0
        self.junction = []
        self.gene_type = []
        self.junction_name_dic = {}
        self.tmp_file_location = tmp_file_location

    def make_genome(self):
        N_number = 100
        polyN = 'N' * N_number
        n = 1
        for circrna in self.circ_rnas:

            # The same annotation description is required three times:
            # gene, transcript and exon.  Build it once and reuse.
            feature_desc = (
                'gene_id "{0}";'
                ' transcript_id "{0}";'
                ' exon_number "1";'
                ' gene_name "{0}";'
                ' gene_source "araport11";'
                ' gene_biotype "protein_coding";'
                ' transcript_source "araport11";'
                ' protein_id "{0}";'
                ' protein_version "1";'
            ).format(circrna.id)
            self.circ_name.extend([feature_desc] * 3)

            # This step aims to cut candidate sequences which are longer than 500bp
            print('No.', n, 'length:', len(circrna))

            if len(circrna) <= 500:
                seq_piece = str(circrna.seq * 2 + polyN)
                self.increase_length += len(circrna) * 2 + N_number
            else:
                seq_piece = str(circrna.seq[-500:] + circrna.seq[:500] + polyN)
                self.increase_length += 1000 + N_number
            self.genome_parts.append(seq_piece)

            if n == 1:
                start_position = 1
                end_position = self.increase_length - N_number
                next_start = self.increase_length + 1
            else:
                start_position = next_start
                end_position = self.increase_length - N_number
                next_start = self.increase_length + 1
            n += 1

            if len(circrna) <= 500:
                jun = int(start_position + len(circrna))
            else:
                jun = int(start_position + 500)
            self.junction.append(jun)
            self.junction_name_dic[jun] = str(circrna.id)

            self.start.extend([start_position] * 3)
            self.end.extend([end_position] * 3)
            self.gene_type.extend(['gene', 'transcript', 'exon'])

        self.genome = ''.join(self.genome_parts)
        final_genome = SeqRecord(Seq(self.genome), id='1_CircularRNA', description='DoubleSeqWith50N')

        self.newGenome.append(final_genome)

        # Dump junction into temp file
        pickle.dump(self.junction, open('{}/junction'.format(self.tmp_file_location), 'wb'))
        pickle.dump(self.junction_name_dic, open('{}/junction_name_dic'.format(self.tmp_file_location), 'wb'))

    def make_gff_file(self):
        self.gff = pd.DataFrame()
        self.gff['start'] = self.start
        self.gff['end'] = self.end
        self.gff['chr'] = 1
        self.gff['araport11'] = 'araport11'
        self.gff['type'] = self.gene_type
        self.gff['.1'] = '.'
        self.gff['+/-'] = '+'
        self.gff['.2'] = '.'
        self.gff['des'] = self.circ_name

        # adjust the order of column of gff
        self.gff = self.gff[['chr', 'araport11', 'type', 'start', 'end', '.1', '+/-', '.2', 'des']]

    def to_fasta(self, name):
        SeqIO.write(self.newGenome, name, 'fasta')

    def to_gff(self, name):
        newlist = [j for j in self.junction for _ in range(3)]
        self.gff['junction'] = newlist
        self.gff.to_csv(name, sep='\t', header=False, index=False, escapechar='"', doublequote=0)
        subprocess.call('sed -i \'s/""/"/g\' {}'.format(name), shell=True)


def main():
    parse = argparse.ArgumentParser(description='This script helps to clean reads and map to genome')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()

    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.safe_load(file)
    tmp_file_location = fileload['tmp_file_location']
    circ_rnas_file = fileload['circrnas']
    name = fileload['genome_name']
    try:
        os.makedirs(tmp_file_location, exist_ok=True)
    except Exception:
        print('Error during make new dir.')
    circ_rnas = SeqIO.parse(circ_rnas_file, 'fasta')
    info = Genome(circ_rnas, tmp_file_location)
    info.make_genome()
    info.make_gff_file()
    info.to_gff('{}/{}.gff'.format(tmp_file_location, name))
    info.to_fasta('{}/{}.fa'.format(tmp_file_location, name))
    print('finished!')


if __name__ == '__main__':
    main()
