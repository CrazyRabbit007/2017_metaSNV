#! /bin/bash

# Script to extract universal genes from the metaSNP output

dir='~/TARA_POPULATION_GENETICS/DATA/TARA_metaSNV/tara_G_metaSNV/'
dirT='~/TARA_POPULATION_GENETICS/DATA/TARA_metaSNV/tara_T_metaSNV/'

mkdir $dir/filtered_match
mkdir $dir/filtered_match/pop/
cd $dir/filtered_match/pop/

LIST=$(ls $dir/filtered/pop/ | grep .freq | cut -f1 -d.)

echo -e "\nList of the $(echo "$LIST"| grep -c .) taxonomic ID found :\n$LIST\n\nFor each we are going to extract the genes expressed."

for TAX_ID in $LIST;
do
echo -e "\n Organism $TAX_ID \n"

head -n 1 $dir/filtered/pop/"$TAX_ID".filtered.freq > temporaire.1
# Get list of genes from metaT
cut -f2 -d: $dirT/filtered/pop/"$TAX_ID".filtered.freq | grep -v TARA_ | grep -v - | sort -u > temporaire.genes

echo -e "We found $(grep -c . temporaire.genes) expressed genes. \n"

grep -f temporaire.genes $dir/filtered/pop/"$TAX_ID".filtered.freq > temporaire.2
cat temporaire.1 temporaire.2 > "$TAX_ID".filtered.freq
rm temporaire*

done
