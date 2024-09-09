#!/usr/bin/env python3

# PURPOSE: 
# COPRO-Seq analysis, generate read count and relative abundance given sequencing reads and reference genomes
#   1. Trim low quality and unpaired reads using trimmomatic
#   2. Generate a single index for reference genomes, quantify reads using kallisto
#   3. Generate reads count and relative abundance table

# EXAMPLE:
# python run_coproseq.py -d /Volumes/ReyLab_QZ10T/COPROseq/data_20240830 -r genome -s metadata_sample.csv -g metadata_genome.csv -o coproseq_out -p 8


import os
import datetime
import argparse
import subprocess
import pandas as pd

from Bio import SeqIO


def run_command(command_string, stdout_path = None):
    
    # checks if a command ran properly. If it didn't, then print an error message then quit
    print('run_coproseq.py, run_command: ' + command_string)
    if stdout_path:
        f = open(stdout_path, 'w')
        exit_code = subprocess.call(command_string, stdout=f, shell=True)
        f.close()
    else:
        exit_code = subprocess.call(command_string, shell=True)

    if exit_code != 0:
        print('run_coproseq.py: Error, the command:')
        print(command_string)
        print('failed, with exit code ' + str(exit_code))
        exit(1)



def run_trimmomatic(data_path, meta_sample, num_processors=4):
    
    # Given raw sequencing read fastq files
    # DO: remove low quality and unpaired reads, trimmed reads in new fastq files

    sample_list = meta_sample['sample_id'].to_list()

    for sample_id in sample_list:
        
        r1_file = meta_sample[meta_sample['sample_id'] == sample_id]['r1_file'].to_list()[0]
        r2_file = meta_sample[meta_sample['sample_id'] == sample_id]['r2_file'].to_list()[0]
        
        # check f1 and f2 file
        if os.path.isfile(r1_file) and os.path.isfile(r2_file):
            
            # run Trimmomatic
            run_command('java -jar /Users/rootqz/bin/trimmomatic-0.39.jar PE -threads {} \
                                    {} \
                                    {} \
                                    {}/{}_R1_trimmed_paired.fastq \
                                    {}/{}_R1_trimmed_unpaired.fastq \
                                    {}/{}_R2_trimmed_paired.fastq \
                                    {}/{}_R2_trimmed_unpaired.fastq \
                                    SLIDINGWINDOW:4:20 MINLEN:50'.format(num_processors, r1_file, r2_file, data_path, sample_id, data_path, sample_id, data_path, sample_id, data_path, sample_id))
            

def rename_genome(genome_path, meta_genome):
    
    # Given genome rederence fna files
    # DO: rename contig seq id by adding genome id
            
    genome_rename = []

    for ncbi_gcf_id in meta_genome['ncbi_gcf_id'].to_list():
        
        genome_id = meta_genome[meta_genome['ncbi_gcf_id'] == ncbi_gcf_id]['genome_id'].to_list()[0]
        genome_fna = genome_path+'/{}_genomic.fna'.format(ncbi_gcf_id)
        
        if os.path.isfile(genome_fna):
            genome_seqs = SeqIO.parse(genome_fna, 'fasta')
        
            for genome_seq in genome_seqs:
                genome_seq.id = genome_id + '|' + genome_seq.id 
                genome_rename = genome_rename + [genome_seq]
                
    SeqIO.write(genome_rename, genome_path+'/genome_ref_renamed.fasta', 'fasta')



def kallisto_index(genome_path):

    # generate kallisto index

    if os.path.isfile(genome_path+'/genome_ref_renamed.fasta'):
        run_command('kallisto index -i {}/genome_ref_renamed.idx {}/genome_ref_renamed.fasta'.format(genome_path,genome_path))



def kallisto_quant(genome_path, data_path, meta_sample, num_processors=4):
    
    # Mapping reads to reference genomes using 'kallisto quant'
            
    sample_list = meta_sample['sample_id'].to_list()
    genome_index = '{}/genome_ref_renamed.idx'.format(genome_path)

    for sample_id in sample_list:
        
        r1_trimmed = data_path+'/{}_R1_trimmed_paired.fastq'.format(sample_id)
        r2_trimmed = data_path+'/{}_R2_trimmed_paired.fastq'.format(sample_id)
        
        # check f1 and f2 file
        if os.path.isfile(r1_trimmed) and os.path.isfile(r2_trimmed):
            
            if not os.path.exists('tmp/kallisto_out_{}'.format(sample_id)):
                run_command('mkdir tmp/kallisto_out_{}'.format(sample_id))
            
            # run kallisto
            run_command('kallisto quant -i {} -o tmp/kallisto_out_{} {} {} -t {}'.format(genome_index, sample_id, r1_trimmed, r2_trimmed, num_processors))



def merge_result(out_path, meta_sample):
    
    # merge kallisto map result

    sample_list = meta_sample['sample_id'].to_list()

    for sample_id in sample_list:
        
        sample_tsv = pd.read_csv('tmp/kallisto_out_{}/abundance.tsv'.format(sample_id), sep = '\t', na_values=None)
        sample_tsv[['genome_id', 'contig_name']] = sample_tsv['target_id'].str.split(pat="|", expand=True)

        sample_abundance = sample_tsv[['genome_id', 'est_counts']].groupby('genome_id').sum()
        sample_abundance.columns = [sample_id]
        
        if sample_id == sample_list[0]:
            abundance_out = sample_abundance
        else:
            abundance_out =  pd.concat([abundance_out, sample_abundance], axis=1)

            
    abundance_out = abundance_out.transpose()
    abundance_out.index.name = 'sample_id'    
    abundance_out.to_csv(out_path+'/mapped_read_count.csv', sep=',', index=True)


    # relative abundance, sum of row to 1
    abundance_out.div(abundance_out.sum(axis=1), axis=0).to_csv(out_path+'/relative_abundance.csv', sep=',', index=True)


def read_stats(out_path, data_path, meta_sample):

    sample_list = meta_sample['sample_id'].to_list()

    read_stats = []
    mapped_read_count = pd.read_csv(out_path+'/mapped_read_count.csv', sep=',')
    mapped_read_count = mapped_read_count.set_index('sample_id')

    for sample_id in sample_list:
        
        sample_seqs = SeqIO.parse(data_path+'/{}_R1_trimmed_paired.fastq'.format(sample_id), 'fastq')
        total_read = sum(1 for _ in sample_seqs)
        mapped_read = sum(mapped_read_count.loc[[sample_id]].transpose()[sample_id].to_list())
        mapped_rate = mapped_read/total_read
        
        read_stats.append([sample_id, total_read, mapped_read, mapped_rate])
        

    read_stats_out = pd.DataFrame(read_stats, columns=['sample_id', 'total_read', 'mapped_read', 'mapped_rate'])
    read_stats_out.to_csv(out_path+'/mapped_read_stats.csv', sep=',',index = False)



def main(args):

    data_path = args.data_directory
    genome_path = args.genome_reference_directory
    meta_sample_file = args.sample_metadata
    meta_genome_file = args.genome_metadata
    out_path = args.output_directory
    num_processors = args.processors

    meta_sample = pd.read_csv(meta_sample_file, sep = ',')
    meta_genome = pd.read_csv(meta_genome_file, sep = ',')

    # check if the output dir exists
    if not os.path.isdir(out_path):
        os.makedirs(out_path)

    # check if the tmp/ dir exists
    if not os.path.isdir('tmp'):
        os.makedirs('tmp')


    with open('coproseq.log', 'w') as log_file:

        # trim raw reads
        print('[{}] Start to trim raw reads using trimmomatic'.format(datetime.datetime.now()), file=log_file)
        run_trimmomatic(data_path, meta_sample, num_processors)
        print('[{}] Reads trimming done'.format(datetime.datetime.now()), file=log_file)

        # rename and merge reference genome fna files
        print('[{}] Rename contigs of reference genomes and concatenate into single fna file'.format(datetime.datetime.now()), file=log_file)
        rename_genome(genome_path, meta_genome)
        print('[{}] Rename done'.format(datetime.datetime.now()), file=log_file)

        # build kallisto index
        print('[{}] Start to build kallisto index'.format(datetime.datetime.now()), file=log_file)
        kallisto_index(genome_path)
        print('[{}] kallisto index done'.format(datetime.datetime.now()), file=log_file)

        # map reads to genome reference
        print('[{}] Start to quatify reads to reference genomes using kallisto quant'.format(datetime.datetime.now()), file=log_file)
        kallisto_quant(genome_path, data_path, meta_sample, num_processors)
        print('[{}] kallisto quant done'.format(datetime.datetime.now()), file=log_file)

        # merge kallisto map result
        print('[{}] Merge reads counts and generate relative abundance table'.format(datetime.datetime.now()), file=log_file)
        merge_result(out_path, meta_sample)

        # get reads mapping stats
        print('[{}] Start to quatify total and mapped reads in each sample'.format(datetime.datetime.now()), file=log_file)
        read_stats(out_path, data_path, meta_sample)
        print('[{}] Reads stats done'.format(datetime.datetime.now()), file=log_file)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('-d', '--data_directory',
                        help='directory path for raw sequencing reads fastq files',
                        type=str,
                        default='')
    parser.add_argument('-r', '--genome_reference_directory',
                        help='directory path for reference genomes fna files',
                        type=str,
                        default='')
    parser.add_argument('-s', '--sample_metadata',
                        help='sample metadata .csv file',
                        type=str,
                        default='')
    parser.add_argument('-g', '--genome_metadata',
                        help='reference genome .csv file',
                        type=str,
                        default='')
    parser.add_argument('-o', '--output_directory',
                        help='output directory',
                        type=str,
                        default='')
    parser.add_argument('-p', '--processors',
                        help='Number of CPU',
                        type=int,
                        default=1)

    args = parser.parse_args()

    main(args)