## COPROSeq

This pipeline is for COPROSeq (community profiling by sequencing) analysis, given sequencing reads and reference genomes to generate read count and relative abundance matrix.


## Requirements

The pipeline was written using Python but requires additional dependencies. 

### Dependencies
1. Download the source from GitHub repo: [COPROseq-master.zip](https://github.com/qijunz/COPROseq/archive/refs/heads/master.zip).
2. Decompress the source and make sure the main scrip `run_coproseq.py` is in your working directory.
3. Download and install [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic). If you download binary, make sure file `trimmomatic-0.39.jar` is in `$PATH`, in some cases this doesn't work (e.g. Mac). Alternatively you can specify `trimmomatic-0.39.jar` file path in `-u`/`--utility_directory` argument.
4. Download and install [Java](https://www.java.com/en/download/help/download_options.html) which is required to run trimmomatic.
5. Download and install [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#obtaining-bowtie-2). Make sure the binary files `bowtie2` is in `$PATH`. You can check using bash command `which bowtie2`.
6. Download and install [samtools](https://www.htslib.org/download/). Make sure the binary files `samtools` is in `$PATH`. You can check using bash command `which samtools`.
7. Download and install [kallisto](https://pachterlab.github.io/kallisto/download).

### Prepare utility files

1. Make sure `trimmomatic-0.39.jar` is in your utility folder.
2. Download [mouse genome bowtie2 index](https://uwmadison.box.com/s/zxf5t0f6hhi3or5ghy4qd5cu1fxkdvra) file and put it in utility folder. Decompress in your utility folder using `tar -zxf Mus_musculus_GRCm38_Rel98_bowtie2_index.tar.gz`.

### Prepare metadata files
This pipeline requires two metadata files.

1. Sample metadata file.
- This `.csv` file should include your sample metadata information. 
- First column is your sample ID, 2nd and 3rd is the absolute path of raw sequencing reads fastq file. One example is shown in `demo/metadata_sample.csv`

2. Reference genome metadata file.
- This `.csv` file should include your reference genome metadata information. The COPROseq is the analysis to map seuqnecing reads to bacteria genomes from a defined community. The whole genome sequencing should be available. 
- First column is genome ID, which will be used for downstream analysis, you can customize genome_id using your preference. 2nd is the RefSeq ID of assembly from NCBI. Usually refseq genome ID start with `GCF`. Check these IDs and make sure they are the same to `.fna` files in genome reference directory. One example is shown in `demo/metadata_genome.csv`


## How to run

### Input files
1. All sequencing fastq files should be in same directory, which will be specified in `-d`/`--data_directory` argument. One example of sequencing files can be downloaded from [demo/data/](https://uwmadison.box.com/s/aiqp2utkqu9p9ydqbm8z3knjqhhh7dq8).
2. All genome fna files (downloaded from NCBI refseq) should be in same directory, which will be specified in `-r`/`--genome_reference_directory` argument. One example of genome files can be downloaded from [demo/genome/](https://uwmadison.box.com/s/ie7qhxt4eibgbuglyu3omjibjose1h13).

### Run script
One example to run:
```bash
$ python run_coproseq.py -r genome -s metadata_sample.csv -g metadata_genome.csv -u /Volumes/ReyLab_QZ10T/COPROseq/util -o coproseq_out -p 8
```

You can check details of arguments by
```bash
$ python run_coproseq.py -h
```

to get

```bash
usage: run_coproseq.py [-h] -r GENOME_REFERENCE_DIRECTORY -s SAMPLE_METADATA
                       -g GENOME_METADATA -o OUTPUT_DIRECTORY -u
                       UTILITY_DIRECTORY [-p PROCESSORS]

optional arguments:
  -h, --help            show this help message and exit
  -r GENOME_REFERENCE_DIRECTORY, --genome_reference_directory GENOME_REFERENCE_DIRECTORY
                        the directory path of reference genomes fna files
  -s SAMPLE_METADATA, --sample_metadata SAMPLE_METADATA
                        sample metadata .csv file
  -g GENOME_METADATA, --genome_metadata GENOME_METADATA
                        reference genome .csv file
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        the directory path to store output files
  -u UTILITY_DIRECTORY, --utility_directory UTILITY_DIRECTORY
                        the directory path of utility files
  -p PROCESSORS, --processors PROCESSORS
                        Number of CPU
```

### Output files

- `mapped_read_count.csv` is the estimated mapped reads count to each genome.
- `mapped_read_stats.csv` is the statistic of reads including the number of total reads, the number of mapped reads and the mapping rate.
- `relative_abundance.csv` is the matrix of relative abundance of each genome in each sample. This is calculated by normalizing sum of row from file `mapped_read_count.csv`. 
- `coproseq.log` is log file recording the run time of each step.
- `tmp/` is the folder to store intermediate files including clean reads and kallisto output files.