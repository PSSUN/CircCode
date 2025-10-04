#/usr/bin/python3

import time
import subprocess
import argparse
import yaml
import pickle
import pandas as pd
import os
from pathlib import Path
from multiprocessing import Pool


# make index of genome and ribo-seq reads
def make_index(thread, genome, ribosome, tmp_file_location, genome_fasta, name):
    """Build genome and rRNA indices required for downstream mapping."""
    tmp_path = Path(tmp_file_location)
    genome = Path(genome)
    ribosome = Path(ribosome)
    genome_fasta = Path(genome_fasta)
    ribo_name = ribosome.stem

    print(
        'STAR --runThreadN {} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {}'.format(
            thread, tmp_path / '', genome
        )
    )

    subprocess.call(
        './requiredSoft/STAR --runThreadN {} --runMode genomeGenerate --genomeDir {} --genomeFastaFiles {}'.format(
            thread, tmp_path / '', genome
        ),
        shell=True,
    )

    subprocess.call(
        'bowtie-build --threads {} {} {}/{}'.format(thread, ribosome, tmp_path, ribo_name),
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT,
    )

    subprocess.call(
        'bowtie-build --threads {} {} {}/{}'.format(thread, genome_fasta, tmp_path, genome_fasta.name),
        shell=True,
    )

    print('bowtie-build --threads {} {} {}/{}'.format(thread, genome_fasta, tmp_path, genome_fasta.name))
    print('-' * 100)
    print(get_time(), 'Make index successfully!')
    print('-' * 100)

def deal_raw_data(
    genome,
    raw_read,
    ribosome,
    thread,
    trimmomatic,
    riboseq_adapters,
    tmp_file_location,
    genome_fasta,
    ribotype,
):
    print(get_time(), 'Start cleaning rawreads...')
    tmp_path = Path(tmp_file_location)
    ribo_name = Path(ribosome).stem

    if ribotype == 'sra':
        read_name = Path(raw_read).stem
        subprocess.call('fastq-dump {} -O {}'.format(raw_read, tmp_file_location), shell=True)
        input_fastq = tmp_path / f'{read_name}.fastq'
    elif ribotype == 'fastq.gz':
        read_name = Path(raw_read).name.replace('.fastq.gz', '')
        input_fastq = raw_read
    else:
        read_name = Path(raw_read).stem
        input_fastq = raw_read

    without_rrna_reads = f'{read_name}.clean.without.rRNA.fastq'
    unmaped_reads = f'{read_name}.clean.without.rRNA.unmapped.fastq'
    clean_fastq = tmp_path / f'{read_name}.clean.fastq'

    print(get_time(), 'Loading reads form:', raw_read)
    print('-' * 100)

    subprocess.call(
        'java -jar {} SE -threads {} -phred33 {} {} '
        'ILLUMINACLIP:{}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16'.format(
            trimmomatic, thread, input_fastq, clean_fastq, riboseq_adapters
        ),
        shell=True,
    )
    		            
    
    # Map clean reads to ribosome sequence by bowtie
    print('remove rRNA')
    subprocess.call(
        'bowtie -p {} -norc --un {} {} {} > {}.map_to_rRNA.sam'.format(
            thread,
            tmp_path / without_rrna_reads,
            tmp_path / ribo_name,
            clean_fastq,
            tmp_path / ribo_name,
        ),
        shell=True,
    )
    print('remove liner sequence')
    print(
        'bowtie -p {} -norc --un {} {} {} > {}.map_to_genome.sam'.format(
            thread,
            tmp_path / unmaped_reads,
            tmp_path / Path(genome_fasta).name,
            tmp_path / without_rrna_reads,
            tmp_path / ribo_name,
        )
    )

    subprocess.call(
        'bowtie -p {} -norc --un {} {} {} > {}.map_to_genome.sam'.format(
            thread,
            tmp_path / unmaped_reads,
            tmp_path / Path(genome_fasta).name,
            tmp_path / without_rrna_reads,
            tmp_path / ribo_name,
        ),
        shell=True,
    )

    print('-' * 100)
    print(get_time(), 'Finished clean process.')
    print('-' * 100)
    global cleanreads
    cleanreads = str(tmp_path / unmaped_reads)

    genome_name = str(genome).split('/')[-1].split('.')[0]
    print(get_time(), 'Start mapping...')
    print('command:')
    print(
        './requiredSoft/STAR --runThreadN {} --outSAMtype BAM SortedByCoordinate --alignIntronMax 10 --genomeDir {} --readFilesIn {} --outFileNamePrefix {}'.format(
            thread,
            tmp_path / '',
            tmp_path / unmaped_reads,
            tmp_path / 'all_bam' / read_name,
        )
    )
    # Path to tophat2 result:
    tophat_result = tmp_file_location+'/'+read_name+'_tophat_result'

    subprocess.call(
        './requiredSoft/STAR --runThreadN {} --outSAMtype BAM SortedByCoordinate --alignIntronMax 10 --genomeDir {} --readFilesIn {} --outFileNamePrefix {}'.format(
            thread,
            tmp_path / '',
            tmp_path / unmaped_reads,
            tmp_path / 'all_bam' / read_name,
        ),
        shell=True,
    )

    print(get_time(), 'Finished mapping')
    print('-'*100)    
    print(get_time(), 'Start analysing...')


def find_reads_on_junction(tmp_file_location, merge_result_name):
    tmp_path = Path(tmp_file_location)
    with open(tmp_path / 'junction_name_dic', 'rb') as f:
        jun_name_dic = pickle.load(f)

    merge_result_file = tmp_path / merge_result_name
    junction_file = tmp_path / 'junction'

    merge_result = pd.read_csv(merge_result_file, sep='\t', header=None, names=['a', 'b', 'c', 'd'])
    with open(junction_file, 'rb') as f:
        junction = pickle.load(f)

    reads_jun = []
    circ_id = []
    frames = []

    for j in junction:
        matches = merge_result[(merge_result.b < j) & (j < merge_result.c)]
        if not matches.empty:
            frames.append(matches)
            reads_jun.append(j)
            circ_id.append(jun_name_dic[j])

    result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['a', 'b', 'c', 'd'])

    print('dump data...')
    try:
        pickle.dump(circ_id, open(tmp_path / f'{merge_result_name}.trans_circ_id', 'wb'))
        pickle.dump(result, open(tmp_path / f'{merge_result_name}.RCRJ_result', 'wb'))
        pickle.dump(reads_jun, open(tmp_path / f'{merge_result_name}.reads_jun', 'wb'))
    except Exception as exc:
        print('Error while dumping RCRJ_result', exc)

    result.to_csv(tmp_path / f'{merge_result_name}.RCRJ_result.csv', sep='\t', header=0, index=False)
    print('find_reads_on_junction finashed. ')

def remove_tmp_file():
    subprocess.call('mkdir -p reads', shell=True)
    subprocess.call('mv *.clean.without.rRNA.fastq ./reads', shell=True)
    subprocess.call('mkdir -p result', shell=True)
    subprocess.call('rm *.fastq', shell=True)


def get_time():
    now = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time.time()))
    strnow = '[{tim}]'.format(tim=now)
    return strnow


def filter_and_map_reads():
    start = time.time()
    make_index()
    deal_raw_data()
    remove_tmp_file()
    stop = time.time()
    print('-' * 100)
    print(get_time(), 'Finished all pipline')
    hours = int(int((stop - start) / 60) / 60)
    print(get_time(), 'Totally run', hours, 'hours')
    print('-' * 100)

    
def bamtobed(bamfile,tmp_file_location):
    subprocess.call('bedtools bamtobed -bed12 -i {} > {}.bam2bedresult.bed'.
		    format(tmp_file_location+'/all_bam/'+bamfile, 
			   tmp_file_location+'/'+bamfile),
		    shell=True)
    subprocess.call('bedtools merge -i {} -c 1 -o count > {}'.
		    format(tmp_file_location+'/'+bamfile+'.bam2bedresult.bed',
			   tmp_file_location+'/'+bamfile+'.merge_result'),
		    shell=True)


def main():
    parse = argparse.ArgumentParser(description='This script helps to clean reads and map to genome')
    parse.add_argument('-y', dest="yaml", required=True)
    args = parse.parse_args()

    yamlfile  = args.yaml
    file      = open(yamlfile)
    fileload  = yaml.load(file, Loader=yaml.FullLoader)
    raw_reads = fileload['raw_reads']
    thread    = fileload['thread']
    ribosome  = fileload['ribosome_fasta']
    trimmomatic = fileload['trimmomatic_jar']
    riboseq_adapters = fileload['riboseq_adapters']
    tmp_file_location = fileload['tmp_file_location']
    genome_fasta = fileload['genome_fasta']
    name = fileload['genome_name']
    merge = fileload['merge']
    ribotype = fileload['ribotype']
    
    
    genome = '{}/{}.fa'.format(tmp_file_location, name)
    
    subprocess.call('mkdir -p {}'.format(tmp_file_location+'/all_bam'),
                    shell=True)

    make_index(thread,
               genome,
               ribosome,
               tmp_file_location,
               genome_fasta,
               name)

    # use multiprocess to deal raw reads
    for raw_read in raw_reads:
        deal_raw_data(genome,raw_read,ribosome,thread,trimmomatic,riboseq_adapters,tmp_file_location,genome_fasta,ribotype)

    if merge == 'T':
        print('-'*20)
        print('merge result...')
        print('-'*20)     
        subprocess.call('samtools merge -f {} {}'.format(tmp_file_location+'/all_bam/total.bam',
                                                         tmp_file_location+'/all_bam/*Aligned.sortedByCoord.out.bam'),
                        shell=True)
        subprocess.call('bedtools bamtobed -bed12 -i {} > {}/bamtobed_result.bed'.format(tmp_file_location+'/all_bam/total.bam',
                                                                                         tmp_file_location),
                        shell=True)
        subprocess.call('bedtools merge -i {}/bamtobed_result.bed -c 1 -o count > {}/merge_result'.format(tmp_file_location,
                                                                                                          tmp_file_location),
                        shell=True)
    else:
        print('-'*20)
        print('not merge...')
        print('-'*20)
        result_bam_list = list(filter(lambda x:x[-4:]=='.bam',
                                      os.listdir(tmp_file_location+'/all_bam')))
        
                                       
        print('bam_list:', result_bam_list)                              
        p2 = Pool(len(result_bam_list))
        for bamfile in result_bam_list:
            p2.apply_async(bamtobed,args=(bamfile,tmp_file_location))
        p2.close()
        p2.join() 
        
    
    if merge == 'T':
        print('analysis junction reads...')
        merge_result_name = 'merge_result'        
        find_reads_on_junction(tmp_file_location,merge_result_name)
    else:
        print('analysis junction reads...')
        

        merge_files = list(filter(lambda x:x[-12:] == 'merge_result', os.listdir(tmp_file_location)))
        p3 = Pool(thread)
        for merge_result_name in merge_files:
            p3.apply_async(find_reads_on_junction,args=(tmp_file_location, merge_result_name))
        p3.close()
        p3.join()         
        
                   
if __name__ == '__main__':
    main()

