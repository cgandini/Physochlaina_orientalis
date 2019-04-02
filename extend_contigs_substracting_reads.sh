!/bin/bash

c1=$1 #reads mapping this fasta will be used to make the subset
c2=$2 #contigs in this fasta will be extended
c3=$3 #reads mapping this fasta will be subtracted from the subset
r1=$4 #reads file 1
r2=$5 #reads file 2
threads=$6 #number of threads to use
readformat=$7 #A, if pair end reads are with /1 and /2, B, if not
isize=$8 #insert size of paired-end reads
PATH2SSAKE=$9 #path to SSAKE folder 

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

rm -r subset

mkdir subset

cd subset

	bwa index $c1
	bwa mem -t $threads -B 5 -O 10 -E 15 $c1 $r1 $r2 > c1_reads_aligned.sam 
	samtools view --threads $threads -h -b -S c1_reads_aligned.sam > c1_reads_aligned.bam 
	samtools view --threads $threads -b -F 4 c1_reads_aligned.bam > c1_reads_mapped.bam 
	samtools view c1_reads_mapped.bam > c1_reads_mapped.txt 
	awk -F'\t' '{print $1}' c1_reads_mapped.txt > c1_reads_mapped_titles.txt 
	perl -ne 'print if ! $x{$_}++' c1_reads_mapped_titles.txt > c1_reads_mapped_titles_uniq.txt 
	
	rm c1_reads_aligned.sam c1_reads_aligned.bam c1_reads_mapped.bam c1_reads_mapped.txt
	
	bwa index $c2
	bwa mem -t $threads -B 5 -O 10 -E 15 -T 120 $c2 $r1 $r2 > c3_reads_aligned.sam 
	samtools view --threads $threads -h -b -S c3_reads_aligned.sam > c3_reads_aligned.bam 
	samtools view --threads $threads -b -F 4 c3_reads_aligned.bam > c3_reads_mapped.bam 
	samtools view c3_reads_mapped.bam > c3_reads_mapped.txt 
	awk -F'\t' '{print $1}' c3_reads_mapped.txt > c3_reads_mapped_titles.txt 
	perl -ne 'print if ! $x{$_}++' c3_reads_mapped_titles.txt > c3_reads_mapped_titles_uniq.txt 
	
	rm c3_reads_aligned.sam c3_reads_aligned.bam c3_reads_mapped.bam c3_reads_mapped.txt
		
if [ $readformat == "A" ]; then
	perl -pe 's/\/\d$//g' c1_reads_mapped_titles_uniq.txt | sort -u > c1_reads_mapped_titles_sorted.txt
	perl -pe 's/\/\d$//g' c3_reads_mapped_titles_uniq.txt | sort -u > c3_reads_mapped_titles_sorted.txt
	awk 'NR == FNR { list[tolower($0)]=1; next } { if (! list[tolower($0)]) print }' c3_reads_mapped_titles_sorted.txt c1_reads_mapped_titles_sorted.txt > reads_mapped_titles.txt 
	perl -pe 's/\n/\/1\n/g' reads_mapped_titles.txt > reads_mapped_titles_1.txt
	perl -pe 's/\n/\/2\n/g' reads_mapped_titles.txt > reads_mapped_titles_2.txt
	seqtk subseq $r1 reads_mapped_titles_1.txt > subset4ssake_readc1.fq 
	seqtk subseq $r2 reads_mapped_titles_2.txt > subset4ssake_reads2.fq
	rm reads_mapped_titles_*.txt
else 
	sort c1_reads_mapped_titles_uniq.txt > c1_reads_mapped_titles_sorted.txt
	sort c3_reads_mapped_titles_uniq.txt > c3_reads_mapped_titles_sorted.txt
	awk 'NR == FNR { list[tolower($0)]=1; next } { if (! list[tolower($0)]) print }' c3_reads_mapped_titles_sorted.txt c1_reads_mapped_titles_sorted.txt > reads_mapped_titles.txt 
	seqtk subseq $r1 reads_mapped_titles.txt > subset4ssake_reads1.fq 
	seqtk subseq $r2 reads_mapped_titles.txt > subset4ssake_reads2.fq
fi

rm c1_reads_mapped_titles_uniq.txt c3_reads_mapped_titles_uniq.txt

cd ..

rm -r extension

mkdir extension

cd extension

	cp ../subset/subset4ssake_reads1.fq  ./r1.fq
	cp ../subset/subset4ssake_reads2.fq  ./r2.fq

	$PATH2SSAKE/tools/TQSfastq.py -f r1.fq -c 30 -t 20
	$PATH2SSAKE/tools/TQSfastq.py -f r2.fq -c 30 -t 20

	$PATH2SSAKE/tools/makePairedOutput2UNEQUALfiles.pl r1.fq_T20C30E33.trim.fa r2.fq_T20C30E33.trim.fa $isize

	mv paired.fa paired.fasta
	mv unpaired.fa unpaired.fasta

	rm r1.fq r2.fq *.trim.fa

	cp $c2 ./contigs2extend.fa

	seqkit split -i -f -O split contigs2extend.fa

	mv ./split/*.fa ./

	rm -r contigs2extend.fa split

	for f in *.fa; do
		$PATH2SSAKE/SSAKE -f paired.fasta -s $f -i 0 -j 15 -u 0 -h 0 -w 5 -m 50 -o 5 -r 0.7 -t 0 -q 0 -y 0 -c 1 -z 100 -p 1 -e 0.75 -k 2 -a 0.5 -x 20 -g unpaired.fasta
		rm *.log *.merged* *.pairing* *.readpo* *.short *.singlets *.csv *.scaff*
		rm $f
	done

	cat *.contigs > extended_contigs.fa

	rm *.contigs