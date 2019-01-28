from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import subprocess
import pickle
import argparse
import yaml


class Genome(object):

    def __init__(self, circrnas, tmp_file_location):
        self.circrnas = circrnas
        self.genome = ''
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
        polyN = 'N' * 50
        n = 1
        for circrna in self.circrnas:
            # try:
            self.circ_name.append(
                'gene_id "{}"; transcript_id "{}"; exon_number "1"; gene_name "{}"; gene_source "araport11"; gene_biotype "protein_coding"; transcript_source "araport11"; protein_id "{}"; protein_version "1";'
                .format(circrna.id, circrna.id, circrna.id, circrna.id))
            self.circ_name.append(
                'gene_id "{}"; transcript_id "{}"; exon_number "1"; gene_name "{}"; gene_source "araport11"; gene_biotype "protein_coding"; transcript_source "araport11"; protein_id "{}"; protein_version "1";'
                    .format(circrna.id, circrna.id, circrna.id, circrna.id))
            self.circ_name.append(
                'gene_id "{}"; transcript_id "{}"; exon_number "1"; gene_name "{}"; gene_source "araport11"; gene_biotype "protein_coding"; transcript_source "araport11"; protein_id "{}"; protein_version "1";'
                    .format(circrna.id, circrna.id, circrna.id, circrna.id))

            # This step aims to cut candidate sequences which are longer than 500bp
            if len(circrna) <= 500:
                self.genome += (circrna.seq * 2 + polyN)
                self.increase_length += len(circrna) * 2 + 100
            else:
                self.genome += (circrna.seq[-500:] + circrna.seq[:500] + polyN)
                self.increase_length += 1050
            # except:
            #     print('Ignored:', circrna.id)
            #     continue
            if n == 1:
                start_position = 1
                end_position = self.increase_length - 100
                next_start = self.increase_length + 1
            else:
                start_position = next_start
                end_position = self.increase_length - 100
                next_start = self.increase_length + 1
            n += 1

            if len(circrna) <= 500:
                jun = int(start_position + len(circrna))
                self.junction.append(start_position + len(circrna))
                self.junction_name_dic[jun] = str(circrna.id)
            else:
                jun = int(start_position + 1000)
                self.junction.append(start_position + 1000)
                self.junction_name_dic[jun] = str(circrna.id)

            self.start.append(start_position)
            self.start.append(start_position)
            self.start.append(start_position)

            self.end.append(end_position)
            self.end.append(end_position)
            self.end.append(end_position)

            self.gene_type.append('gene')
            self.gene_type.append('transcript')
            self.gene_type.append('exon')

        final_genome = SeqRecord(Seq(str(self.genome)), id='1_CircularRNA', description='DoubleSeqWith50N')

        self.newGenome.append(final_genome)

        # Dump junction into tem file
        pickle.dump(self.junction, open('{}/junction'.format(tmp_file_location), 'wb'))
        pickle.dump(self.junction_name_dic, open('junction_name_dic', 'wb'))

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
        self.gff.to_csv(name, sep='\t', header=False, index=False, escapechar='"', doublequote=0)
        subprocess.call('sed -i \'s/""/"/g\' {}'.format(name), shell=True)


def main():
    parse = argparse.ArgumentParser(description='This script helps to clean reads and map to genome')
    parse.add_argument('-y', dest="yaml", required=True)
    parse.add_argument('-map-only', dest='map_only', required=False)
    args = parse.parse_args()

    yamlfile = args.yaml
    file = open(yamlfile)
    fileload = yaml.load(file)
    tmp_file_location = fileload['tmp_file_location']
    circrnas = fileload['circrnas']
    name = fileload['genome_name']
    info = Genome(circrnas)
    info.make_genome(tmp_file_location)
    info.make_gff_file()
    info.to_gff('{}/{}.gff'.format(tmp_file_location, name))
    info.to_fasta('{}/{}.fa'.format(tmp_file_location, name))


if __name__ == '__main__':
    main()