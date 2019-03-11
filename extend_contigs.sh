#!/bin/bash

mt=$1 #mitochondrial contigs for subset
mt2=$2 #mitochondrial contigs for extension
cp=$3 #chloroplast or chloroplast contigs
r1=$4 #reads file 1
r2=$5 #reads file 2
threads=$6 #number of threads to use
PATH2SSAKE=$7 #PATH to SSAKE folder
readformat=$8 #A if pair end reads are with /1 and /2, B if pair end reads are with 1:N:0 and 2:N:0

if ! [ -x "$(command -v bwa)" ]; then
  echo 'Error: bwa is not installed.' >&2
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

mkdir subset

cd subset

	bwa index $mt
	bwa mem -t $threads -B 5 -O 10 -E 15 $mt $r1 $r2 > mt_reads_aligned.sam 
	samtools view --threads $threads -h -b -S mt_reads_aligned.sam > mt_reads_aligned.bam 
	samtools view --threads $threads -b -F 4 mt_reads_aligned.bam > mt_reads_mapped.bam 
	samtools view mt_reads_mapped.bam > mt_reads_mapped.txt 
	awk -F'\t' '{print $1}' mt_reads_mapped.txt > mt_reads_mapped_titles.txt 
	perl -ne 'print if ! $x{$_}++' mt_reads_mapped_titles.txt > mt_reads_mapped_titles_uniq.txt 
	sort mt_reads_mapped_titles_uniq.txt > mt_reads_mapped_titles_sorted.txt
	
	rm mt_reads_aligned.sam mt_reads_aligned.bam mt_reads_mapped.bam mt_reads_mapped.txt
	
	bwa index $cp
	bwa mem -t $threads -B 5 -O 10 -E 15 -T 120 $cp $r1 $r2 > cp_reads_aligned.sam 
	samtools view --threads $threads -h -b -S cp_reads_aligned.sam > cp_reads_aligned.bam 
	samtools view --threads $threads -b -F 4 cp_reads_aligned.bam > cp_reads_mapped.bam 
	samtools view cp_reads_mapped.bam > cp_reads_mapped.txt 
	awk -F'\t' '{print $1}' cp_reads_mapped.txt > cp_reads_mapped_titles.txt 
	perl -ne 'print if ! $x{$_}++' cp_reads_mapped_titles.txt > cp_reads_mapped_titles_uniq.txt 
	sort cp_reads_mapped_titles_uniq.txt > cp_reads_mapped_titles_sorted.txt
	
	rm cp_reads_aligned.sam cp_reads_aligned.bam cp_reads_mapped.bam cp_reads_mapped.txt
	
	awk 'NR == FNR { list[tolower($0)]=1; next } { if (! list[tolower($0)]) print }' cp_reads_mapped_titles_sorted.txt mt_reads_mapped_titles_sorted.txt > reads_mapped_titles.txt 
	
if [ $8 == "A" ]; then
	seqtk subseq $r1 reads_mapped_titles.txt > subset4ssake_reads1.fq 
	seqtk subseq $r2 reads_mapped_titles.txt > subset4ssake_reads2.fq
else 
	seqtk subseq $r1 reads_mapped_titles.txt | perl -pe 's/ 1:N:0/\/1/g' > subset4ssake_reads1.fq 
	seqtk subseq $r2 reads_mapped_titles.txt | perl -pe 's/ 1:N:0/\/1/g' > subset4ssake_reads2.fq
fi  

cd ..

mkdir SSAKE_extension

cd SSAKE_extension

	cp ../subset/subset4ssake_reads1.fq  ./r1.fq
	cp ../subset/subset4ssake_reads2.fq  ./r2.fq

	$PATH2SSAKE/tools/TQSfastq.py -f r1.fq -c 30 -t 20
	$PATH2SSAKE/tools/TQSfastq.py -f r2.fq -c 30 -t 20

	$PATH2SSAKE/tools/makePairedOutput2UNEQUALfiles.pl r1.fq_T20C30E33.trim.fa r2.fq_T20C30E33.trim.fa 755

	mv paired.fa paired.fasta
	mv unpaired.fa unpaired.fasta

	rm r1.fq r2.fq *.trim.fa

	cp $mt2 ./contigs2extend.fa

	seqkit split -i -f -O split contigs2extend.fa

	mv ./split/*.fa ./

	rm -r contigs2extend.fa split

	for f in *.fa; do
		$PATH2SSAKE/SSAKE -f paired.fasta -s $f -i 0 -j 15 -u 0 -h 0 -w 5 -m 50 -o 5 -r 0.7 -t 0 -q 0 -y 0 -c 1 -z 100 -p 1 -e 0.75 -k 2 -a 0.5 -x 20 -g unpaired.fasta
		rm *.log & rm *.merged* & rm *.pairing* & rm *.readpo* & rm *.short & rm *.singlets & rm *.csv & rm *.scaff*
		rm $f
	done

	cat *.contigs > extended_contigs.fa

	rm *.contigs