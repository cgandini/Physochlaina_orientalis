#!/bin/bash

f=$1 #contig fasta file
r1=$2 #reads file 1
r2=$3 #reads file 2
threads=$4 #number of threads to use
readformat=$5 #A, if pair end reads are with /1 and /2, B, if not

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
bowtie2-align-s -a -p $threads --local --very-sensitive-local -q -x "${f%%.*}" -1 $r1 -2 $r2 --no-unal -S "${f%%.*}".sam 
samtools view --threads $threads -h -b -S "${f%%.*}".sam | samtools sort > "${f%%.*}".bam


if [ $readformat == "A" ]; then
	samtools view "${f%%.*}".bam | awk -F'\t' '{print $1}' | perl -pe 's/\/\d$//g' | perl -ne 'print if ! $x{$_}++' > reads_mapped_titles.txt 
	perl -pe 's/\n/\/1\n/g' reads_mapped_titles.txt > reads_mapped_titles_1.txt
	perl -pe 's/\n/\/2\n/g' reads_mapped_titles.txt > reads_mapped_titles_2.txt
	seqtk subseq $r1 reads_mapped_titles_1.txt > "${f%%.*}"_subset_1.fq
	seqtk subseq $r2 reads_mapped_titles_2.txt > "${f%%.*}"_subset_2.fq
	rm reads_mapped_titles_*.txt
else 
	samtools view "${f%%.*}".bam | awk -F'\t' '{print $1}' | perl -ne 'print if ! $x{$_}++' > reads_mapped_titles_uniq.txt 
	seqtk subseq $r1 reads_mapped_titles.txt > "${f%%.*}"_subset_1.fq
	seqtk subseq $r2 reads_mapped_titles.txt > "${f%%.*}"_subset_2.fq
fi

rm *.sam
