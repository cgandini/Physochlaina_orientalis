#!/bin/bash

if [ "$1" == "-h" ] ; then
    echo -e "Usage: "$(basename $0)" [fasta file] [reads 1] [reads 2] [# of threads to use] [read format; A, if pair end reads names end with /1 and /2, B, if not]"
    exit 1
fi

f="$(cd $(dirname "$1") && pwd -P)/$(basename "$1")" #contig fasta file
r1="$(cd $(dirname "$2") && pwd -P)/$(basename "$2")" #reads file 1
r2="$(cd $(dirname "$3") && pwd -P)/$(basename "$3")" #reads file 2
threads=$4 #number of threads to use
readformat=$5 #A, if pair end reads names end with /1 and /2, B, if not
fasta_name="$(basename "$1")"

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

mkdir subset_bowtie

cd subset_bowtie

cp $f ./fasta.fa

bowtie2-build-s fasta.fa fasta 
bowtie2-align-s -a -p $threads --local --very-sensitive-local -q -x fasta -1 $r1 -2 $r2 --no-unal -S fasta.sam 
samtools view --threads $threads -h -b -S fasta.sam | samtools sort > fasta.bam


if [ $readformat == "A" ]; then
	samtools view fasta.bam | awk -F'\t' '{print $1}' | perl -pe 's/\/\d$//g' | perl -ne 'print if ! $x{$_}++' > reads_mapped_titles.txt 
	perl -pe 's/\n/\/1\n/g' reads_mapped_titles.txt > reads_mapped_titles_1.txt
	perl -pe 's/\n/\/2\n/g' reads_mapped_titles.txt > reads_mapped_titles_2.txt
	seqtk subseq $r1 reads_mapped_titles_1.txt > fasta_subset_1.fq
	seqtk subseq $r2 reads_mapped_titles_2.txt > fasta_subset_2.fq
	rm reads_mapped_titles_*.txt
else 
	samtools view fasta.bam | awk -F'\t' '{print $1}' | perl -ne 'print if ! $x{$_}++' > reads_mapped_titles.txt 
	seqtk subseq $r1 reads_mapped_titles.txt | perl -pe 's/(@.*) .*/$1\/1/g' > "${fasta_name%%.*}"_subset_1.fq
	seqtk subseq $r2 reads_mapped_titles.txt | perl -pe 's/(@.*) .*/$1\/2/g' > "${fasta_name%%.*}"_subset_2.fq
fi

rm *.sam *.bam *.bt2
