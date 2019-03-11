# *Physochlaina orientalis*

## Scripts used during the assembly and analysis of the *P. orientalis* mitochondrial genome ##

  The scripts posted here were used to analyze the mitochondrial genome of the Solanaceae *Phyoschlaina orientalis*. Scripts are more a combination of different programs than a program itself. So you need to check out if you have installed all pre-requisites in your path before running them. All scripts are made in bash and tested in mac OS, so there is no warranty they worked correctly on Linux. Nonetheless, they should work with minimal changes. 
  
   Please be noticed that I am not a professional bioinformatic, and therefore these scripts are extremally "home-made". However, they do the job. If you used them please cite: "XXX"  

For any doubts write me to: gandini.carolin@gmail.com

### 1- extend_contigs.sh: Get a subset of reads and extend using SSAKE (deleting plastidial reads)

**PRE-REQUISITES (used versions are within parenthesis):**
  
  - BWA (0.7.15-r1140): https://sourceforge.net/projects/bio-bwa/files/
  - samtools (1.4): https://sourceforge.net/projects/samtools/files/
  - seqtk (1.2-r101-dirty): https://github.com/lh3/seqtk
  
### 2- subset_reads.sh: Get a subset of reads

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
  
  **INPUTS:** you will need to create a folder with all fasta files to analyze (.fa extension). Then, you just run in the terminal
  
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
    
