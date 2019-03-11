# *Physochlaina orientalis*

## Scripts used during the assembly and analysis of the *P. orientalis* mitochondrial genome ##

  The scripts posted here were used to analyze the mitochondrial genome of the Solanaceae *Phyoschlaina orientalis*. Scripts are more a combination of different programs than a program itself. So you need to check out if you have installed all pre-requisites in your path before running them. All scripts are made in bash and tested in mac OS, so there is no warranty they worked correctly on Linux. Nonetheless, they should work with minimal changes. 
  
   Please be noticed that I am not a professional bioinformatic, and therefore these scripts are extremally "home-made". However, they do the job. If you used them please cite: "XXX"  

For any doubts write to: gandini.carolin@gmail.com

### 1- extend_contigs.sh: Get a subset of reads and extend using SSAKE (deleting plastidial reads)

  This script allows you to extend contigs individually using a mitochondrial subset of reads. 
  
**PRE-REQUISITES (used versions are within parenthesis):**
  
  - BWA (0.7.15-r1140): https://sourceforge.net/projects/bio-bwa/files/
  - samtools (1.4): https://sourceforge.net/projects/samtools/files/
  - seqtk (1.2-r101-dirty): https://github.com/lh3/seqtk
  - SSAKE (3.8.5): http://www.bcgsc.ca/platform/bioinfo/software/ssake **do not add to PATH**
  
**INPUTS:**
  
  You should first create the following files (please see Materials and Methods of XXX for more information):
  
  - mt4subse: mitochondrial contigs for the subset in fasta format, this should be all contigs with BLAST mitochondrial hits plus all contigs within the mitochondrial coverage.
  - mt4extension: mitochondrial contigs for the extension in fasta format, all contigs with mitochondrial BLAST mitochondrial hits
  - cp: chloroplast genome or chloroplast contigs in fasta format
  - reads1: pair-end reads file 1 in fastq format
  - reads2: pair-end reads file 2 in fastq format
  - threads: number of threads to use
  - reads_format: indicate the format of reads names within the fastq file: **A**, if pair reads are denoted as read_name/1 and read_name/2 or **B**, if pair reads are denoted as read_name 1:N:0 and read_name 2:N:0. You can check it by doing: 

```  
head reads1.fastq
```
**Run in the terminal**

```  
mkdir [extension folder]
cd PATH/[extension folder]
PATH/extend_contigs.sh PATH/mt4subset PATH/mt4extension PATH/cp PATH/reads1 PATH/reads2 threads reads_format
```

**OUTPUTS:**

The file extended_contigs.fa contain contigs after the extension process. 

### 2- subset_reads.sh: Get a subset of reads

This script allows you to get a subset of reads for any fasta or multifasta file. This is useful to reduce memory consumption of other programs and therefore allows you to work easily in a conventional PC. 

  **PRE-REQUISITES (used versions are within parenthesis):**
  
  - bowtie2 (2.2.2): http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  - samtools (1.4): https://sourceforge.net/projects/samtools/files/
  - seqtk (1.2-r101-dirty): https://github.com/lh3/seqtk

**INPUTS:**

  - fasta: fasta or multifasta file for which you need the subset
  - reads1: pair-end reads file 1 in fastq format
  - reads2: pair-end reads file 2 in fastq format
  - threads: number of threads to use

**Run in the terminal**

```  
mkdir [subset folder]
cd PATH/[subset folder]
PATH/subset_reads.sh PATH/fasta PATH/reads1 PATH/reads2 threads
```

**OUTPUTS:**

  - fasta_name_subset_1.fq and fasta_name_subset_2.fq: subset of reads in fastq files 
  - fasta_name.bam: bam file of aligned reads

### 3- get_repeats.sh: Analyze repeats and short repeats (< 100 bp) families using VSEARCH

  This script allows you to analyze multiple fasta files (.fa). Each species should be in a different file. If the species have 2 or more chromosomes (or scaffolds), these should be placed in a multifasta file.  
  
  **PRE-REQUISITES (used versions are within parenthesis):**
  
  - BLAST+ (2.7.1+): https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
  - trf (4.09): https://tandem.bu.edu/trf/trf.download.html
  - seqkit (0.10.1): https://github.com/shenwei356/seqkit
  - vsearch (2.8.1): https://github.com/torognes/vsearch
  - mafft (7.305b): https://mafft.cbrc.jp/alignment/software/
  - bedtools (2.26.0): https://bedtools.readthedocs.io/en/latest/content/installation.html
  - bedops (2.4.35): https://bedops.readthedocs.io/en/latest/
  
  **INPUTS:** you will need to create a folder with all fasta files to analyze (.fa extension). Then, 
  
**Run in the terminal**
  
```  
cd PATH/[folder with files]
PATH/get_repeats.sh
```

  **OUTPUTS:** for each fasta a folder would be created. Within you will find:

  * chr folder: with BLAST and Tandem Repeat Finder results
  * clusters_SR folder: with short repeat (<100 bp) clustering analysis
  * repeats_stats.txt: resume of all types of repeats and clustering analysis (including genome coverages)
  * BED files for (remember that in BED files start coords start at 0). Except for SR_with_cluster, names for each repeat represent the start from pair 1 and the start from pair 2:
    * all_coords: all repeats
    * LR_coords: large repeats
    * IR_coords: intermediate repeats
    * SR_coords: short repeats
    * SR_with_cluster: short repeats with cluster name in which was assigned
    * TR_coords: tandem repeats
 * repeat_table.txt: a resume table of all analyzed species
    
