#/usr/bin/python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import subprocess
import pickle
import argparse
import yaml


class Genome(object):

    def __init__(self, circ_rnas, tmp_file_location):
        self.circ_rnas = circ_rnas
        self.genome = ''
        self.genome_1 = ''
        self.genome_2 = ''
        self.genome_3 = ''
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

            # The same behavior was repeated three times because they
            # correspond to genes, transcripts, and exons, respectively.
            # Writing in this way can look more clear, as well as the following.

            # -----gene-----
            self.circ_name.append(
                'gene_id "{}";'
                ' transcript_id "{}";'
                ' exon_number "1";'
                ' gene_name "{}";'
                ' gene_source "araport11";'
                ' gene_biotype "protein_coding";'
                ' transcript_source "araport11";'
                ' protein_id "{}";'
                ' protein_version "1";'
                .format(circrna.id, circrna.id, circrna.id, circrna.id))

            # -----transcript-----
            self.circ_name.append(
                'gene_id "{}";'
                ' transcript_id "{}";'
                ' exon_number "1";'
                ' gene_name "{}";'
                ' gene_source "araport11";'
                ' gene_biotype "protein_coding";'
                ' transcript_source "araport11";'
                ' protein_id "{}";'
                ' protein_version "1";'
                .format(circrna.id, circrna.id, circrna.id, circrna.id))

            # -----exon-----
            self.circ_name.append(
                'gene_id "{}";'
                ' transcript_id "{}";'
                ' exon_number "1";'
                ' gene_name "{}";'
                ' gene_source "araport11";'
                ' gene_biotype "protein_coding";'
                ' transcript_source "araport11";'
                ' protein_id "{}";'
                ' protein_version "1";'
                .format(circrna.id, circrna.id, circrna.id, circrna.id))

            # This step aims to cut candidate sequences which are longer than 500bp
            print('No.',n,'length:',len(circrna))
            
            # time.sleep(2)
            # time was used to debug here for developer.
            
            # The length of each list is set the limitation of 30000
            # If list is too large, it will cost much time to append() sequence.
            
            # The length of each circRNA was no more than 500bp
            # Because we only focus on the junction site of each circRNA
            # If the length is longer than 500bp, the longer part of circRNA will be cut off.
            if n <= 30000:
                if len(circrna) <= 500:
                    self.genome += (circrna.seq * 2 + polyN)
                    self.increase_length += len(circrna) * 2 + N_number
                else:
                    self.genome += (circrna.seq[-500:] + circrna.seq[:500] + polyN)
                    self.increase_length += (1000+N_number)
            elif 30000 < n <= 60000:
                if len(circrna) <= 500:
                    self.genome_1 += (circrna.seq * 2 + polyN)
                    self.increase_length += len(circrna) * 2 + N_number
                else:
                    self.genome_1 += (circrna.seq[-500:] + circrna.seq[:500] + polyN)
                    self.increase_length += (1000+N_number)
            elif 60000 < n <= 90000:
                if len(circrna) <= 500:
                    self.genome_2 += (circrna.seq * 2 + polyN)
                    self.increase_length += len(circrna) * 2 + N_number
                else:
                    self.genome_2 += (circrna.seq[-500:] + circrna.seq[:500] + polyN)
                    self.increase_length += (1000+N_number)
            else:
                if len(circrna) <= 500:
                    self.genome_3 += (circrna.seq * 2 + polyN)
                    self.increase_length += len(circrna) * 2 + N_number
                else:
                    self.genome_3 += (circrna.seq[-500:] + circrna.seq[:500] + polyN)
                    self.increase_length += (1000+N_number)



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
                self.junction.append(start_position + len(circrna))
                self.junction_name_dic[jun] = str(circrna.id)
            else:
                jun = int(start_position + 500)
                self.junction.append(start_position + 500)
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
        self.genome = self.genome+self.genome_1+self.genome_2+self.genome_3
        final_genome = SeqRecord(Seq(str(self.genome)), id='1_CircularRNA', description='DoubleSeqWith50N')

        self.newGenome.append(final_genome)

        # Dump junction into tem file
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
        newlist = []
        for i in self.junction:
            newlist.append(i)
            newlist.append(i)
            newlist.append(i)
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
    fileload = yaml.load(file)
    tmp_file_location = fileload['tmp_file_location']
    circ_rnas_file = fileload['circrnas']
    name = fileload['genome_name']
    try:
        subprocess.call('mkdir -p {}'.format(tmp_file_location), shell=True)
    except:
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
