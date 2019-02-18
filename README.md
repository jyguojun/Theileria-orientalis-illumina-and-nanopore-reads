# Evaluation of *Theileria orientalis* genome assembly methods using nanopore sequencing and analysis of variation between genomes
This repository contains:
1) all codes and commands used in the methodology
2) supplementary results

## Abstract
*Theileria orientalis* (Apicomplexa: Piroplasmida) is a non-lymphocyte transforming tick-borne haemoparasite of cattle that causes ill-thrift and anaemia. In recent years, clinical outbreaks of *T. orientalis* caused by the pathogenic genotype of this parasite (Ikeda) have been increasingly observed throughout the Asia Pacific.  Currently, there are no  available vaccines for this disease, although a live vaccine based on a benign genotype (Buffeli) has been proposed to provide cross-protection against Ikeda. However, our recent genomic studies using illumina short reads of three *T. orientalis* genotypes, Ikeda, Buffeli and a low pathogenic genotype (Chitose) have revealed substantial genetic divergence, perhaps at the species level. As short read technology is unable to effectively resolve the structure of the genomes, we continued to investigate the isolates using previously generated Illumina short reads combined with nanopore long reads. In this study, we sequenced the three isolates with a R9.4.1 MinION flow cell and tested four different hybrid assembly methods, Flye, Canu, Unicycler and Masurca. Flye and Canu assemblies were further processed with Nanopolish and five iterations of Pilon using Illumina reads. Different combinations of the assemblers were trialed and evaluated in order to determine the best pipeline for *T. orientalis* genome assembly. Evaluations with Quast and MUMmer revealed Unicycler to be the best assembler for *T. orientalis* Ikeda, and Flye for genotypes, Chitose and Buffeli. Alignments to the *T. orientalis* (Shintoku) reference sequence revealed potential structural variation in the apathogenic Buffeli genotype. The detailed methodology and results from this study will be presented and discussed, including the genome annotation and findings of the variation between the pathogenic and apathogenic *T. orientalis* types.

## Commands used in this chapter

### Basecalling nanopore reads
**[Albacore](https://community.nanoporetech.com/protocols/albacore-offline-basecalli/v/abec_2003_v1_revan_29nov2016/linux)** (Requires Nanopore account)

`read_fast5_basecaller.py -i fast5_input_directory -s output_directory -f FLO-MIN106 -k SQK-LSK108 -o fastq -t 24`

### Trimming nanopore reads
**[Porechop](https://github.com/rrwick/Porechop)**

`porechop -i input_albacore_reads.fastq -o output_porechopped_reads.fastq`

### Long read assembly (Nanopore reads only)
**[Canu](https://github.com/marbl/canu)**

`canu -p output_prefix -d canu_output_directory genomeSize=9.0m -nanopore-raw albacore_reads.fastq useGrid=false`

**[Flye](https://github.com/fenderglass/Flye)**

Normal coverage:

`flye --nano-raw albacore_reads.fastq -g 9m -o flye_output_directory -i 5 -t 24`

40X reduced coverage for initial contig assembly:

`flye --nano-raw albacore_reads.fastq -g 9m -o flye_output_directory_40x -i 5 -t 24 --asm-coverage 40`

### Hybrid assembly (Nanopore and illumina reads)
**[Unicycler](https://github.com/rrwick/Unicycler)**

`unicycler -1 illumina_1.fastq -2 illumina_2.fastq -l porechopped_reads.fastq -o output_directory --linear_seqs 6 -t 24`

**[Masurca](https://github.com/alekseyzimin/masurca)**

Masurca requires a configuration file to run (use # to ignore parameters in configuration file). To run Masurca `bash assemble.sh`, use the shell script `assemble.sh` generated from `masurca configuration_file`

Example of the configuration file used in this assembly:
```
# DATA is specified as type {PE,JUMP,OTHER,PACBIO} and 5 fields:
# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads
# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be
# innies, i.e. --->.<---, and JUMP are assumed to be outties
# <---.--->. If there are any jump libraries that are innies, such as
# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads
# are optional for PE libraries and mandatory for JUMP libraries. Any
# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first
# converted into Celera Assembler compatible .frg files (see
# http://wgs-assembler.sourceforge.com)
DATA
#Illumina paired end reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#if single-end, do not specify <reverse_reads>
#MUST HAVE Illumina paired end reads to use MaSuRCA
PE= i2 362 145 illumina_1.fastq illumina_2.fastq
#Illumina mate pair reads supplied as <two-character prefix> <fragment mean> <fragment stdev> <forward_reads> <reverse_reads>
#JUMP= sh 3600 200  /FULL_PATH/short_1.fastq  /FULL_PATH/short_2.fastq
#pacbio OR nanopore reads must be in a single fasta or fastq file with absolute path, can be gzipped
#if you have both types of reads supply them both as NANOPORE type
#PACBIO=/FULL_PATH/pacbio.fa
NANOPORE=porechopped_reads.fastq
#Other reads (Sanger, 454, etc) one frg file, concatenate your frg files into one if you have many
#OTHER=/FULL_PATH/file.frg
END

PARAMETERS
#set this to 1 if your Illumina jumping library reads are shorter than 100bp
EXTEND_JUMP_READS=0
#this is k-mer size for deBruijn graph values between 25 and 127 are supported, auto will compute the optimal size based on the read data and GC content
GRAPH_KMER_SIZE = auto
#set this to 1 for all Illumina-only assemblies
#set this to 0 if you have more than 15x coverage by long reads (Pacbio or Nanopore) or any other long reads/mate pairs (Illumina MP, Sanger, 454, etc)
USE_LINKING_MATES = 0
#specifies whether to run mega-reads correction on the grid
#USE_GRID=0
#specifies queue to use when running on the grid MANDATORY
#GRID_QUEUE=all.q
#batch size in the amount of long read sequence for each batch on the grid
#GRID_BATCH_SIZE=300000000
#use at most this much coverage by the longest Pacbio or Nanopore reads, discard the rest of the reads
LHE_COVERAGE=25
#set to 1 to only do one pass of mega-reads, for faster but worse quality assembly
MEGA_READS_ONE_PASS=0
#this parameter is useful if you have too many Illumina jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms 
LIMIT_JUMP_COVERAGE = 300
#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically. 
#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms.
CA_PARAMETERS =  cgwErrorRate=0.15
#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if Illumina coverage >100
KMER_COUNT_THRESHOLD = 1
#whether to attempt to close gaps in scaffolds with Illumina data
CLOSE_GAPS=1
#auto-detected number of cpus to use
NUM_THREADS = 28
#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage
JF_SIZE = 72000000
#set this to 1 to use SOAPdenovo contigging/scaffolding module.  Assembly will be worse but will run faster. Useful for very large (>5Gbp) genomes from Illumina-only data
SOAP_ASSEMBLY=0
END
```

### Polishing reads(Nanopore and illumina reads)
Only Flye assemblies were subjected to polishing. To polish the reads, a combination of programs (**[BWA](https://github.com/lh3/bwa)**, **[Nanopolish](https://github.com/jts/nanopolish)** and **[Pilon](https://github.com/broadinstitute/pilon)** were used. 

`nanopolish index -d fast5_files_directory albacore_reads.fasta`

`bwa index flye_scaffolds.fasta`

`bwa mem -x ont2d -t 20 flye_scaffolds.fasta albacore_reads.fasta | samtools sort -o reads.sorted.bam -T reads.tmp`

`samtools index reads.sorted.bam`

`python nanopolish_makerange.py flye_scaffolds.fasta | parallel --results nanopolish.results -P 8 \
	nanopolish variants --consensus -o flye_assembly.polished.{1}.vcf -w {1} -r albacore_reads.fasta -g flye_scaffolds.fasta -b reads.sorted.bam -t 6 --min-candidate-frequency 0.1 -q dcm,dam`

`nanopolish vcf2fasta -g flye_scaffolds.fasta flye_assembly.polished.*.vcf > ikeda_flye_polished_genome.fa`
