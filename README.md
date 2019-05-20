# *Physochlaina orientalis* assembly and analysis

## Scripts used during the assembly and repeat analysis of the *P. orientalis* mitochondrial genome ##

  The scripts here were used to analyze the mitochondrial genome of the Solanaceae *Phyoschlaina orientalis*. Scripts are more a combination of different programs than a program itself. So you need to check out if you have installed all pre-requisites in your path before running them. All scripts are made in bash and tested in mac OS, so there is no warranty they worked correctly on Linux. Nonetheless, they should work with minimal changes. 
  
   Please be noticed that I am not a professional bioinformatic, and therefore these scripts are extremally "home-made". However, they do the job. If you used them please cite: *"The Complete Organelle Genomes of Physochlaina orientalis: Insights into Short Sequence Repeats across Seed Plant Mitochondrial Genomes"*  

For any doubts write to: gandini.carolin@gmail.com

## Installation

Just download the complete folder. Give permissions to all scripts within and add files to the bash profile. Then you can just run calling the script.

## 1- extend_contigs_subtracting_reads.sh: get a subset of reads subtracting unwanted reads and extend using SSAKE

  This script allows you to extend contigs individually by first subtract unwanted reads. 
  
**PRE-REQUISITES (used versions are within parenthesis):**
  
  - BWA (0.7.15-r1140): https://sourceforge.net/projects/bio-bwa/files/
  - samtools (1.4): https://sourceforge.net/projects/samtools/files/
  - seqtk (1.2-r101-dirty): https://github.com/lh3/seqtk
  - SSAKE (3.8.5): http://www.bcgsc.ca/platform/bioinfo/software/ssake **do not add to PATH**
  
**INPUTS:**
  
  You should first create the following files (please see Materials and Methods of *"The Complete Organelle Genomes of Physochlaina orientalis: Insights into Short Sequence Repeats across Seed Plant Mitochondrial Genomes"* for more information):
  
  - fasta file for subset: reads mapping this fasta (or multifasta) file will be used to make the subset
  - fasta file for extension: contigs in this dasta file will be extended
  - fasta file for substract reads: reads mapping this fasta (or multifasta) file will be subtracted from the subset
  - reads 1: paired-end reads file 1 in fastq format
  - reads 2: paired-end reads file 2 in fastq format
  - threads: number of threads to use
  - read format: indicate the format of reads names within the fastq file: **A**, if paired reads are denoted as read_name/1 and read_name/2 or **B**, if not. You can check it by doing: 

```  
head reads1.fastq
```

  - insert size: indicate the insert size of paired-end reads
  - PATH to SSAKE folder: **complete path to SSAKE folder inside the script**

**Run in the terminal**

```  
mkdir [extension folder]
cd PATH/[extension folder]

PATH/extend_contigs_subtracting_reads.sh -h ## print help

"Usage: extend_contigs_subtracting_reads.sh [fasta file for subset] [fasta file for extension] [fasta file for subtract reads] [reads 1] [reads 2] [# of threads to use] [read format; A, if pair end reads names end with /1 and /2, B, if not] [insert size]

```

**OUTPUTS:**

You will 2 folders:
  
  - subset
    - mt_reads_mapped_titles_sorted.txt: names of reads mapping to the mitochondrial genome
    - cp_reads_mapped_titles_sorted.txt: names of reads mapping to the chloroplast genome
    - reads_mapped_titles.txt: names of reads after substraction of reads mapping the chloroplast genome
    - subse4ssake_reads1.fq and subse4ssake_reads2.fq: subset of reads used for ssake extension
    
  - extension
    - extended_contigs.fa: fasta file with contigs after the extension process
    - extension.txt: information about the extension length of each contig

## 2- extend_contigs_withsubset.sh

This script is similar to the previous one, however, no substraction is done. 

**Run in the terminal**

```  
mkdir [extension folder]
cd PATH/[extension folder]

PATH/extend_contigs_withsubset.sh -h

"Usage: extend_contigs_withsubset.sh [fasta file for subset] [fasta file for extension] [reads 1] [reads 2] [# of threads to use] [read format; A, if pair end reads names end with /1 and /2, B, if not] [insert size]

```
    
## 3- subset_reads_bowtie.sh: get a subset of reads for any fasta or multifasta file

This script allows you to get a subset of reads for any fasta or multifasta file. This is useful to reduce memory consumption of other programs (as for example Consed) and therefore allows you to work easily in a conventional PC. 

  **PRE-REQUISITES (used versions are within parenthesis):**
  
  - bowtie2 (2.2.2): http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  - samtools (1.4): https://sourceforge.net/projects/samtools/files/
  - seqtk (1.2-r101-dirty): https://github.com/lh3/seqtk

**INPUTS:**

  - fasta: fasta or multifasta file for which you need the subset
  - reads 1: paired-end reads file 1 in fastq format
  - reads 2: paired-end reads file 2 in fastq format
  - threads: number of threads to use
  - read format: indicate the format of reads names within the fastq file: **A**, if paired reads are denoted as read_name/1 and read_name/2 or **B**, if not. You can check it by doing: 

```  
head reads1.fastq
```

**Run in the terminal**

```  
PATH/subset_reads_bowtie.sh -h

"Usage: subset_reads_bowtie.sh [fasta file] [reads 1] [reads 2] [# of threads to use] [read format; A, if pair end reads names end with /1 and /2, B, if not]"
```

**OUTPUTS:**

  - fasta_name_subset_1.fq and fasta_name_subset_2.fq: subset of reads in fastq files 
  - fasta_name.bam: bam file of aligned reads

## 4- get_repeats.sh: analyze repeats and short repeats (< 100 bp) families using VSEARCH

  This script reports genome content of all repeat types, defined as tandem repeats (TR), short repeats (SR: < 100 bp), intermediate repeats (IntR/IR: 100-1000 bp) and large repeats (LR: > 1000bp). In addition, a clustering/family analysis of SR is done using the algorithm included in the VSEARCH tool. 
  This script allows you to analyze multiple species. However, each species should be in a different fasta file (.fa or .fasta). If the species have 2 or more chromosomes (or scaffolds), these should be placed in a multifasta file.  
  
  **PRE-REQUISITES (used versions are within parenthesis):**
  
  - BLAST+ (2.7.1+): https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  - trf (4.09): https://tandem.bu.edu/trf/trf.download.html
  - seqkit (0.10.1): https://github.com/shenwei356/seqkit
  - vsearch (2.8.1): https://github.com/torognes/vsearch
  - mafft (7.305b): https://mafft.cbrc.jp/alignment/software/
  - bedtools (2.26.0): https://bedtools.readthedocs.io/en/latest/content/installation.html
  - bedops (2.4.35): https://bedops.readthedocs.io/en/latest/
  - R
  
  **INPUTS:** you will need to create a folder with all fasta files to analyze. Then, indicate:
  
  - threads: number of threads to use
  - mincluster: minumun number of sequences in a cluster to re-cluster [integer]
  - getal: get aligments using mafft? [Y or N]
  - clusterall: cluster all analyzed species? [Y or N], usefull if you have more than one species
  
**Run in the terminal**
  
```  
cd PATH/[folder with fasta files]
PATH/get_repeats.sh -h

"cd to folder with genomes in fasta files. Each genome should be in individual files. If a genome is composed by 2 or more chromosomes/scaffolds should be placed in a multifasta file

Usage: get_repeats.sh [# of threads to use] [minimum number of sequences in a cluster to re-cluster] [get alignments? Y or N] [cluster all? Y or N]"
```

  **OUTPUTS:** for each fasta a folder will be created. Within you will find:

  * chr folder: with BLAST and Tandem Repeat Finder results
  * clusters_SR folder: with short repeat (<100 bp) clustering analysis
  * repeats_stats.txt: resume of all types of repeats and clustering analysis (including genome coverages)
  * BED files for (remember that in BED files start coords start at 0). Names for each repeat represent the start from repeat pair 1 and the start from repeat pair 2:
      * all_coords: all repeats
      * LR_coords: large repeats
      * IR_coords: intermediate repeats
      * SR_coords: short repeats
      * SR_with_cluster: short repeats with cluster name in which was assigned
      * TR_coords: tandem repeats
      
      
 * repeat_table.txt: a resume table of all analyzed species
    
