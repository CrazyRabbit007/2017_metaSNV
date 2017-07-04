#! /bin/bash

LIST=$(ls /science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered/pop/ | grep .freq | cut -f1 -d.)

echo -e "\nList of the $(echo "$LIST"| grep -c .) taxonomic ID found :\n$LIST\n\nFor each we are going to extract the common genes."

for TAX_ID in $LIST;
do
echo -e "\n Organism $TAX_ID \n"
perl -lane 'print unless $_ =~ /\t-1/' /science/paolil/TARA_metaSNP/tara_subset_genes_metaSNP/filtered_universal/pop/"$TAX_ID".universal.filtered.freq > "$TAX_ID".universal.common.filtered.freq

COMMON_POS=$(perl -lane 'print unless $_ =~ /\t-1/' /science/paolil/TARA_metaSNP/tara_subset_genes_metaSNP/filtered_universal/pop/"$TAX_ID".universal.filtered.freq | cut -f1 | cut -f2 -d: | grep -v "^$" | grep -c .)
POS=$(cut -f1 /science/paolil/TARA_metaSNP/tara_subset_genes_metaSNP/filtered_universal/pop/"$TAX_ID".universal.filtered.freq | cut -f2 -d: | grep -v "^$" | grep -c .)


echo -e "We found $COMMON_POS common positions out of $POS positions. \n \n  **********";

done
