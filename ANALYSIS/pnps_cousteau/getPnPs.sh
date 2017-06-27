#! /bin/bash

# Compute per species and per gene PnPs
echo 'Load Python'
ml Python
echo 'Load BEDTools'
ml BEDTools

directory=$1
db_ann=$2 # /directory/freeze9 or /directory/freeze11.marine.regenomes
db_fasta=$3

directory=/science/paolil/TARA_metaSNV/tara_1_metaSNV
db_ann=/nfs/home/paolil/metaSNV/db/freeze11.marine.repgenomes/freeze11.marine.repgenomes.annotations.txt
db_fasta=/nfs/home/paolil/metaSNV/db/freeze11.marine.repgenomes/freeze11.marine.repgenomes
#directory=/science/paolil/GUT_metaSNV/GUT_3_metaSNV
#db_ann=/nfs/home/paolil/metaSNV/db/freeze_DB_metaSNP/freeze4.annotations.txt
#db_fasta=/nfs/home/paolil/metaSNV/db/freeze_DB_metaSNP/freeze4.genomes.RepGenomesv4.fna

mkdir $directory/pnps
mkdir $directory/pnps/fasta
mkdir $directory/pnps/genes

#files=$(ls $directory/filtered/pop/ | cut -f1 -d'.')
#files=1471461
files=1096769

#loop over species
for species in $files;
do
echo $species
#cut -f1 $directory/filtered/pop/$species.filtered.freq | cut -f2 -d':' | grep -v - | grep -v ^$ | sort -u > $directory/pnps/genes/$species.genes
#grep -f $directory/pnps/genes/$species.genes $db_ann | awk 'BEGIN { FS = OFS = "\t" } { print $3,$7,$8+1,$3 ":" $2,0, $9}' > $directory/pnps/genes/$species.bed
#bedtools getfasta -fi $db_fasta -bed $directory/pnps/genes/$species.bed -s -name > $directory/pnps/fasta/$species.fa
python computePnPs.py $directory/filtered/pop/$species.filtered.freq
done

#python computePnPs_full.py $directory/filtered/pop/ gut_samples


