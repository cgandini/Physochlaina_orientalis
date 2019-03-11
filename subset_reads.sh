#!/bin/bash

f=$1 #contig fasta file
R1=$2 #reads file 1
R2=$3 #reads file 2
threads=$4 #number of threads to use

if ! [ -x "$(command -v bowtie2)" ]; then
  echo 'Error: bowtie2 is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v samtools)" ]; then
  echo 'Error: samtools is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v seqtk)" ]; then
  echo 'Error: seqtk is not installed.' >&2
  exit 1
fi

bowtie2-build-s $f "${f%%.*}" 
bowtie2-align-s -a -p $threads --local --very-sensitive-local -q -x "${f%%.*}" -1 $R1 -2 $R2 --no-unal -S "${f%%.*}".sam 
samtools view --threads $threads -h -b -S "${f%%.*}".sam | samtools sort > "${f%%.*}".bam
samtools view "${f%%.*}".bam | awk -F'\t' '{print $1}' | perl -ne 'print if ! $x{$_}++' > reads_mapped_titles_uniq.txt 
seqtk subseq $R1 reads_mapped_titles_uniq.txt > "${f%%.*}"_subset_1.fq
seqtk subseq $R2 reads_mapped_titles_uniq.txt > "${f%%.*}"_subset_2.fq

rm *.sam
