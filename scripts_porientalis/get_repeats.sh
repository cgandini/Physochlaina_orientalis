#!/bin/bash

if [ "$1" == "-h" ] ; then
    echo -e "cd to folder with genomes in fasta files. Each genome should be in individual files. If a genome is composed by 2 or more chromosomes/scaffolds should be placed in a multifasta file\nUsage: "$(basename $0)" [# of threads to use] [minimum number of sequences in a cluster to re-cluster] [get alignments? Y or N] [cluster all? Y or N]"
    exit 1
fi

threads=$1
mincluster=$2 ## minimum number of sequences in a cluster to re-cluster
getal=$3 ## Y or N
clusterall=$4 ## Y or N

SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if ! [ -x "$(command -v makeblastdb)" ]; then
  echo 'Error: makeblastdb is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v blastn)" ]; then
  echo 'Error: blastn is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v trf)" ]; then
  echo 'Error: trf is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v seqkit)" ]; then
  echo 'Error: seqkit is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v vsearch)" ]; then
  echo 'Error: vsearch is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v mafft)" ]; then
  echo 'Error: mafft is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v bedtools)" ]; then
  echo 'Error: bedtools is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v bedops)" ]; then
  echo 'Error: bedops is not installed.' >&2
  exit 1
fi

if ! [ -x "$(command -v R)" ]; then
  echo 'Error: R is not installed.' >&2
  exit 1
fi

for f in *.fa; do ## each genome should be in individual files. If a genome is composed by 2 or more chromosomes or scaffolds should be placed in a multifasta file. 

name=$(head -1 $f | perl -pe 's/>//g')
rm -rf ./"${name}" repeat_table.txt
mkdir ./"${name}"
mkdir ./"${name}"/chr
cp $f ./"${name}"/chr/"${name}".fa

done
 
for D in ./*/; do

cd $D

cat ./chr/*.fa | seqkit seq -u -w 0 > genome.fa

## Calculate repeats

makeblastdb -in genome.fa -dbtype nucl -parse_seqids

cd chr

for f in *.fa; do

blastn -word_size 7 -ungapped -perc_identity 80 -evalue 0.001 -query $f -db ../genome.fa -out "${f%%.*}"_blast.txt -outfmt "6 qstart qend sstart send sstrand length pident evalue bitscore mismatch qcovhsp qseqid sseqid"

name=$(head -1 $f | perl -pe 's/>//g')

trf $f 2 7 7 80 10 50 2000 -h -m -d

cat *.dat | while read line; do
	if [[ "$line" =~ ^Sequence.* ]]; then
		chr=$(echo $line | perl -pe 's/Sequence: //g') 
  elif [[ "$line"  =~ ^[0-9] ]]; then
  	 echo $chr $line >> TR_coords.tmp
  fi
 done

done

cd ..

name=$(head -1 genome.fa | perl -pe 's/>//g')
length=$(seqkit fx2tab --length genome.fa --name --only-id | awk '{sum += $2} END {print sum}')

cat ./chr/*_blast.txt > repeats_blast.txt
cat ./chr/TR_coords.tmp | cut -f1-3 -d ' ' | perl -pe 's/ /\t/g' > TR_coords.bed
rm ./chr/TR_coords.tmp

awk -v OFS='\t' '{if($3<$4) print $0; else print $1,$2,$4,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' repeats_blast.txt | awk -v OFS='\t' '{if($1<$3) print $0; else print $3,$4,$1,$2,$5,$6,$7,$8,$9,$10,$11,$13,$12}' | sort -k1,1n | awk '!seen[$1,$2,$3,$4,$5,$12,$13]++' | awk '{if($11<100) print $0}' | awk '{if(!($1==$3 && $2==$4)) print $0}' > repeats_blast.sorted.txt

awk -v OFS='\t' '{if($6<100) print $0}' repeats_blast.sorted.txt | awk -v OFS='\t' '{if($5=="minus") print $12,$1-1,$2,"1","1","-\n"$13,$3-1,$4,"2","2","-"; else print $12,$1-1,$2,"1","1","+\n"$13,$3-1,$4,"2","2","+"}' | sort -k1,1 -k2,2n > SR_coords.sorted.bed

seqkit subseq --bed SR_coords.sorted.bed genome.fa | perl -pe 's/\:.*//g' | perl -pe 's/-/_/g' | seqkit rmdup -n > SR.fa

rm SR_coords.sorted.bed

## Calculate clusters for all sequences

mkdir clusters_SR

cd clusters_SR

vsearch --cluster_fast ../SR.fa --threads $threads --centroids SR_unique.fa --id 1 --sizeout --strand both --uc clusters_unique.uc --sizeorder --minseqlength 20 --log vsearch_UN_CL.txt --clusters UN_CL_ 

vsearch --cluster_size SR_unique.fa --threads $threads --centroids centroids.fa --id 0.80 --sizeout --sizein --strand both --uc clusters.uc --sizeorder --minseqlength 20 --log vsearch_CL.txt --clusters CL_

for f in CL*; do

	seqs=$(grep "size" $f | perl -pe 's/.*size=(\d*)/$1/g' | awk '{sum += $1} END {print sum}') 
	echo -e ""${f}"\t"${seqs}"" >> CL_abundance.txt

done

sort -k2,2nr CL_abundance.txt > cluster_abundance.txt

awk -v mincluster=$mincluster '{if($2>=mincluster) print $0}' cluster_abundance.txt > clusters_min"${mincluster}"seq.txt

while IFS= read -r line; do
	a=( $line )
	CL=$(echo ${a[0]})
	seqs=$(echo ${a[1]})
	cp "${CL}" cluster_"${seqs}".fa
done < clusters_min"${mincluster}"seq.txt

clusters_minseq=$(cat clusters_min"${mincluster}"seq.txt | wc -l)

if [ "$getal" == "Y" ]; then

if [ "$clusters_minseq" -gt 0 ]; then

	mkdir alignments

	for f in cluster_*.fa; do

		name=$(head -1 ../genome.fa | perl -pe 's/>//g')
		CL=$(echo $f | perl -pe 's/CL_(\d*)_(\d*)\.fa.*/$1/g')
		seqs=$(echo $f | perl -pe 's/CL_(\d*)_(\d*)\.fa.*/$2/g')
		perl -pe "s/>.*(\;size=\d*)/>${name}_${seqs}_${CL}\1/g" $f > "${f}".rename.fa
		mafft --adjustdirectionaccurately "${f}".rename.fa > ./alignments/"${f}".align.fa
		rm $f

	done

fi

fi

# Get repeat bed_file with clusters names

	for f in CL_*; do

	grep ">" $f | perl -pe 's/>(.*);size=(\d*)/$1\t$2/g' > seq_"${f}".txt
	CL1=$(echo $f | perl -pe 's/CL_//g')
	size1=$(awk -v CL1=$CL1 '{if($1=="C" && $2==CL1) print $3}' clusters.uc)

			while IFS= read -r line; do
				a=($line)
				seq=$(echo ${a[0]})
				size2=$(echo ${a[1]})
				CL2=$(awk -v size2=$size2 -v seq=$seq '{if($1=="C" && $3==size2 && $9==seq) print $0}' clusters_unique.uc | cut -f2)
				cat UN_CL_"${CL2}" | perl -pe "s/>(.*)_(\d*)_(\d*);size=\d*/>\1_\2_\3_${CL1}_${size1}/g" >> SR_repeats.fa
				cat UN_CL_"${CL2}" | grep ">" | perl -pe 's/>(.*)_(\d*)_(\d*);size=\d*/$1\t$2\t$3/g' | awk -v size1=$size1 -v CL1=$CL1 '{print $0"\t"size1"_"CL1}' | uniq >> SR_repeats.bed
			done < seq_"${f}".txt

	done

if [ "$clusters_minseq" -gt 0 ]; then

seqkit fx2tab SR_repeats.fa | perl -pe 's/(.*_\d*_\d*_\d*)_(\d*)/$1\t$2/g' | awk -F '\t' -v mincluster=$mincluster -v OFS='\t' '{if($2>=mincluster) print $1"_"$2,$3}' | seqkit tab2fx > clusters_min"${mincluster}"seq.fa

fi

mv clusters_unique.uc vsearch_1round.uc
mv clusters.uc vsearch_2round.uc

rm cluster_*.fa

for f in seq_*.txt; do

	rm $f

done 

for f in CL*; do 

	rm $f

done 

for f in UN_CL*; do 

	rm $f

done 

cd ..

## Get repeat stats

awk -v OFS='\t' '{if($6>=1000) print $0}' repeats_blast.sorted.txt  | awk -v OFS='\t' '{if($5=="minus") print $12,$1-1,$2,$1"_"$3,"LR","+\n"$13,$3-1,$4,$1"_"$3,"LR","-"; else print $12,$1-1,$2,$1"_"$3,"LR","+\n"$13,$3-1,$4,$1"_"$3,"LR","+"}' | sort -k1,1 -k2,2n | awk '{if($2>=0) print $0;}' > LR_coords.bed 

awk -v OFS='\t' '{if($6<1000 && $6>=100) print $0}' repeats_blast.sorted.txt | awk -v OFS='\t' '{if($5=="minus") print $12,$1-1,$2,$1"_"$3,"IR","+\n"$13,$3-1,$4,$1"_"$3,"IR","-"; else print $12,$1-1,$2,$1"_"$3,"IR","+\n"$13,$3-1,$4,$1"_"$3,"IR","+"}' | sort -k1,1 -k2,2n | awk '{if($2>=0) print $0;}' > IR_coords.bed 

awk -v OFS='\t' '{if($6<100) print $0}' repeats_blast.sorted.txt | awk -v OFS='\t' '{if($5=="minus") print $12,$1-1,$2,$1"_"$3,"SR","+\n"$13,$3-1,$4,$1"_"$3,"SR","-"; else print $12,$1-1,$2,$1"_"$3,"SR","+\n"$13,$3-1,$4,$1"_"$3,"SR","+"}' | sort -k1,1 -k2,2n | awk '{if($2>=0) print $0;}' > SR_coords.bed 

cat *_coords.bed | sort -k1,1 -k2,2n | cut -f1-3 | awk '{if($3>0) print $0}' > all_coords.bed

cat IR_coords.bed SR_coords.bed TR_coords.bed | sort -k1,1 -k2,2n | cut -f1-3 > IR_SR_TR_coords.bed

echo -e "\n######## genome length ########\n" > repeats_stats.txt

echo -e "${PWD##*/}" >> repeats_stats.txt
echo -e "chromosomes/scaffolds: $(grep ">" genome.fa | wc -l)" >> repeats_stats.txt
echo -e "total (bp):\t$(seqkit fx2tab --length genome.fa --name --only-id | awk '{sum += $2} END {print sum}')"  >> repeats_stats.txt

rm genome.*

echo -e "\n######## total number of repeat pairs ########\n" >> repeats_stats.txt
echo -e "total:\t$(cat repeats_blast.sorted.txt | wc -l)" >> repeats_stats.txt
echo -e "unique:\t$(awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' repeats_blast.sorted.txt| sort | uniq | wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i all_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i all_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt

echo -e "\n######## repeat pairs >= 1000 ########\n" >> repeats_stats.txt
echo -e "total:\t$(awk '{if($11<100 && $6>=1000) print $0}' repeats_blast.sorted.txt | wc -l)" >> repeats_stats.txt
echo -e "unique:\t$(awk '{if($11<100 && $6>=1000) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | sort | uniq| wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i LR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i LR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt

echo -e "\n######## repeat pairs 100-1000 ########\n" >> repeats_stats.txt
echo -e "total:\t$(awk '{if($11<100 && $6<1000 && $6>=100) print $0}' repeats_blast.sorted.txt | wc -l )" >> repeats_stats.txt
echo -e "unique:\t$(awk '{if($11<100 && $6<1000 && $6>=100) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | sort | uniq | wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i IR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i IR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')">> repeats_stats.txt

echo -e "\n######## repeat pairs < 100 ########\n" >> repeats_stats.txt
echo -e "total:\t$(awk '{if($11<100 && $6<100) print $0}' repeats_blast.sorted.txt | wc -l )" >> repeats_stats.txt
echo -e "unique:\t$(awk '{if($11<100 && $6<100) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | sort | uniq | wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt

echo -e "\n######## tandem repeats ########\n" >> repeats_stats.txt
echo -e "total:\t$(cat TR_coords.bed | wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i TR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i TR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt

echo -e "\n######## repeat pairs < 1000 ########\n" >> repeats_stats.txt
echo -e "total:\t$(awk '{if($11<1000 && $6<1000) print $0}' repeats_blast.sorted.txt | wc -l )" >> repeats_stats.txt
echo -e "unique:\t$(awk '{if($11<1000 && $6<1000) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | sort | uniq | wc -l)" >> repeats_stats.txt
echo -e "coverage (bp):\t$(bedtools merge -i IR_SR_TR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "coverage (%):\t$(bedtools merge -i IR_SR_TR_coords.bed  | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')"  >> repeats_stats.txt

echo -e "\n######## coverage intersections ########\n" >> repeats_stats.txt
echo -e "TR n SR (bp):\t$(bedops --intersect TR_coords.bed SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "TR n SR (%):\t$(bedops --intersect TR_coords.bed SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt
echo -e "SR n IR (bp):\t$(bedops --intersect IR_coords.bed SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "SR n IR (%):\t$(bedops --intersect IR_coords.bed SR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt
echo -e "(TR,SR,IR) n LR (bp):\t$(bedops --intersect IR_SR_TR_coords.bed LR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2}' | awk '{sum += $4} END {print sum}')" >> repeats_stats.txt
echo -e "(TR,SR,IR) n LR (%):\t$(bedops --intersect IR_SR_TR_coords.bed LR_coords.bed | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk -v genome_l=$length '{sum += $4} END {print sum/genome_l*100}')" >> repeats_stats.txt


echo -e "\n######## number/coverage of repeats per size ########\n" >> repeats_stats.txt

tot_R=$(awk '{if($11<100) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | wc -l | perl -pe 's/ //g')

for f in 20-39 40-59 60-79 80-99 100-199 200-299 300-499 500-999 1000-10000000; do 

	awk -v OFS='\t' -v size="${f}" '{split(size, subf, "-"); if($6>=subf[1] && $6<=subf[2]) print $12,$1,$2"\n"$13,$3,$4}' repeats_blast.sorted.txt | sort -k1,1 -k2,2n > "${f}".tmp
	echo -e "total "${f}":\t$(cat "${f}".tmp | wc -l)" >> repeats_stats.txt
	echo -e ""${f}" (%):\t$(cat "${f}".tmp | wc -l | perl -pe 's/ //g' | awk -v tot_R=$tot_R '{print $0/tot_R*100}')" >> repeats_stats.txt
	echo -e "coverage "${f}" (bp):\t$(bedtools merge -i "${f}".tmp | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}' )" >> repeats_stats.txt
	rm *.tmp

done

echo -e "\n######## number/coverage of repeats <100 per size ########\n" >> repeats_stats.txt

SR_pairs=$(awk '{if($11<100 && $6<100) print $0}' repeats_blast.sorted.txt | awk -v OFS="\t" '{print $1"_"$2"\n"$3"_"$4}' | wc -l | perl -pe 's/ //g')

for f in 20-39 40-59 60-79 80-99; do 

	awk -v OFS='\t' -v size="${f}" '{split(size, subf, "-"); if($6>=subf[1] && $6<=subf[2]) print $12,$1,$2"\n"$13,$3,$4}' repeats_blast.sorted.txt | sort -k1,1 -k2,2n > "${f}".tmp
	echo -e "total "${f}":\t$(cat "${f}".tmp | wc -l)" >> repeats_stats.txt
	echo -e ""${f}" (SR %):\t$(cat "${f}".tmp | wc -l | perl -pe 's/ //g' | awk -v SR=$SR_pairs '{print $0/SR*100}')" >> repeats_stats.txt
	echo -e "coverage "${f}" (bp):\t$(bedtools merge -i "${f}".tmp | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}' )" >> repeats_stats.txt

rm *.tmp

done

echo -e "\n######## number/coverage of repeats <100 per identity ########\n" >> repeats_stats.txt

for f in 80-85 85-90 90-95 95-101; do 

	awk -v OFS='\t' -v pident="${f}" '{split(pident, subf, "-"); if($11<100 && $6<100 && $7>=subf[1] && $7<subf[2]) print $12,$1,$2"\n"$13,$3,$4}' repeats_blast.sorted.txt | sort -k1,1 -k2,2n > "${f}".tmp
	echo -e "total "${f}":\t$(cat "${f}".tmp | wc -l)" >> repeats_stats.txt
	echo -e ""${f}" (SR %):\t$(cat "${f}".tmp | wc -l | perl -pe 's/ //g' | awk -v SR=$SR_pairs '{print $0/SR*100}')" >> repeats_stats.txt
	echo -e "coverage "${f}" (bp):\t$(bedtools merge -i "${f}".tmp | awk -v OFS="\t" '{print $0,$3-$2+1}' | awk '{sum += $4} END {print sum}' )" >> repeats_stats.txt

rm *.tmp

done

## Get cluster stats

SR=$(grep ">" SR.fa | wc -l | perl -pe 's/ //g')

echo -e "\n######## repeat clustering ########\n" >> repeats_stats.txt
echo -e "SR to cluster:\t$(grep ">" SR.fa | wc -l | awk '{if ($0=="") print "0"; else print $0}')" >> repeats_stats.txt
echo -e "SR after cluster 100% identicals:\t$(grep ">" ./clusters_SR/SR_unique.fa | wc -l | awk '{if ($0=="") print "0"; else print $0}')\t$(grep ">" ./clusters_SR/SR_unique.fa | wc -l | perl -pe 's/ //g' | awk '{if ($0=="") print "0"; else print $0}' | awk -v SR=$SR '{print $0/SR*100}')" >> repeats_stats.txt
echo -e "SR in largest cluster:\t$(cut -f2 ./clusters_SR/cluster_abundance.txt | head -1 | awk '{if ($0=="") print "0"; else print $0}')\t$(cut -f2 ./clusters_SR/cluster_abundance.txt | head -1 | awk '{if ($0=="") print "0"; else print $0}' | awk -v SR=$SR '{print $0/SR*100}')" >> repeats_stats.txt

for f in 5000 1000 500 200 100 50 20 10; do

	echo -e "SR in cluster (>"${f}"):\t$(awk -v size="${f}" '{if($2>size) sum += $2} END {print sum}' ./clusters_SR/cluster_abundance.txt | awk '{if ($0=="") print "0"; else print $0}')\t$(awk -v size="${f}" '{if($2>size) sum += $2} END {print sum}' ./clusters_SR/cluster_abundance.txt | awk '{if ($0=="") print "0"; else print $0}' | awk -v SR=$SR '{print $0/SR*100}')" >> repeats_stats.txt

done

echo -e "10 largest clusters:\n$(head -10 ./clusters_SR/cluster_abundance.txt | awk '{if ($0=="") print "0"; else print $0}')" >> repeats_stats.txt

rm repeats_blast.* IR_SR_TR_coords.bed

cd ..

done

## Get repeat table 

echo -e "species name\t#chromosomos\tlength (bp)\tTotR pairs\t# of TotR\tTotR cov (%)\tLR pairs\t# of LR\tLR cov (%)\tIR pairs\t# of IR\tIR cov (%)\tSR pairs\t# of SR\tSR cov (%)\tSR in largest cluster (%)\tSR in cluster >1000 sequences (%)\tSR in cluster >50 sequences (%)\tSR in cluster >10 sequences\tTR\tTR cov (%)\tTR n SR (%)\tSR n IR (%)\t(TR,SR,IR) n LR\t20-39 bp (TotR %)\t40-59 bp (TotR %)\t60-79 bp (TotR %)\t80-99 bp (TotR %)\t100-199 bp (TotR %)\t200-299 bp (TotR %)\t300-499 bp (TotR %)\t500-999 bp (TotR %)\t>1000 bp (TotR %)\t20-39 bp (SR %)\t40-59 bp (SR %)\t60-79 bp (SR %)\t80-99 bp (SR %)\t80-85 pident (SR %)\t85-90 pident (SR %)\t90-95 pident (SR %)\t95-100 pident (SR %)" > repeat_table.txt

for f in ./*/repeats_stats.txt; do

echo -e "$(awk 'NR==4' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==5' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==6' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==10' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==11' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==13' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==17' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==18' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==20' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==24' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==25' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==27' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==31' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==32' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==34' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==122' $f | cut -f3 | perl -pe 's/.*:\s+//g')\t$(awk 'NR==124' $f | cut -f3 | perl -pe 's/.*:\s+//g')\t$(awk 'NR==128' $f | cut -f3 | perl -pe 's/.*:\s+//g')\t$(awk 'NR==130' $f | cut -f3 | perl -pe 's/.*:\s+//g')\t$(awk 'NR==38' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==40' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==52' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==54' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==56' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==61' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==64' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==67' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==70' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==73' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==76' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==79' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==82' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==85' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==91' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==94' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==97' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==100' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==106' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==109' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==112' $f | perl -pe 's/.*:\s+//g')\t$(awk 'NR==115' $f | perl -pe 's/.*:\s+//g')" >> repeat_table.txt

done

if [ "$clusterall" == "Y" ]; then

## Calculate clusters interspecies

rm -r clusterall

mkdir clusterall

cd clusterall

cat ../*/clusters_SR/clusters_min"${mincluster}"seq.fa > all_repeats_min"${mincluster}"seqs.fa

vsearch --cluster_fast all_repeats_min"${mincluster}"seqs.fa --threads $threads --centroids SR_unique.fa --id 0.80 --sizein --sizeout --strand both --uc clusters_unique.uc --sizeorder --minseqlength 20 --log vsearch.txt --clusters CL_ 

mkdir combinations

for f in CL*; do

	seqs=$(grep "size" $f | perl -pe 's/.*size=(\d*)/$1/g' | awk '{sum += $1} END {print sum}') 
	echo -e ""${f}"\t"${seqs}"" | sort -k2,2nr >> CL_abudance.txt
	echo "${f%%.*}" > resume.txt
	grep ">" $f | perl -pe 's/>(.*)_\d*_\d*_(\d*)_(\d*)\;.*/$1_$3_$2/g' | sort | uniq -c | perl -pe 's/^ *([^ ]+) +/\1\t/'  >> resume.txt

	CL=$(echo "${f%%.*}")

	grep ">" $f |perl -pe 's/>(.*)_\d*_\d*_(\d*)_(\d*)\;.*/$1_$3_$2/g' | sort | uniq -c | perl -pe 's/^ *([^ ]+) +/\1\t/' | awk -v CL=$CL -v OFS="\t" '{print $0,CL}' >> seqs_in_cluster.txt

	grep ">" $f | perl -pe 's/>(.*)_\d*_\d*_(\d*)_(\d*)\;.*/$1_$3_$2/g' | sort | uniq -c |

	awk -v CL=$CL -v OFS="\t" '
	        {
	                A[++c] = $2
	        }
	        END {
	                for ( i = 1; i <= c; i++ )
	                {
	                        for ( j = 1; j <= c; j++ )
	                        {
	                                print A[j], A[i], CL
	                        }
	                }
	        }
	' > ./combinations/"${f%%.*}".txt

	rm $f

done


for f in ./combinations/*.txt; do 

         cat $f | perl -pe 's/>//g' >> combinations_all.tmp
         cat $f | perl -pe 's/>//g' | perl -pe 's/_\d*_\d*\t/\t/g' | cut -f1-2 | awk '{if ($1!=$2) print $0}' >> combinations_interspecies.tmp
         cat $f | perl -pe 's/>//g' | perl -pe 's/_(\d*_\d*)\t/\t$1\t/g' | cut -f1-4 | awk '{if ($1==$3) print $0}' >> combinations_intraspecies.tmp

done

for f in *.tmp; do

sort -k1,1 $f |  awk '{ if (!seen[$0]++) print $0}' > "${f%%.*}".txt
rm $f

done

perl -pe 's/_(\d*_\d*)\t/\t$1\t/g' combinations_all.txt | cut -f1-4 | awk -v OFS='\t' '{if($1==$3) print $1"_"$2,$3"_"$4,"same"; else print $1"_"$2,$3"_"$4,"other"}' |  awk '{ if (!seen[$0]++) print $0}' > combinations_all_R.txt

cut -f1 combinations_all_R.txt | uniq| perl -pe 's/_(\d*)_(\d*)\n/\t$1\t$2\n/g' | sort -k1,1 -k2,2nr | perl -pe 's/\t/_/g' > levels_R.txt

Rscript $SCRIPTPATH/tools/get_cluster_graph.R combinations_all_R.txt levels_R.txt

rm combinations_all_R.txt levels_R.txt resume.txt seqs_in_cluster.txt 

cd ..

fi



