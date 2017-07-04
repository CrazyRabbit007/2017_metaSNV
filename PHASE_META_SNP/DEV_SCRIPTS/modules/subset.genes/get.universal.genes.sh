#! /bin/bash

# Script to extract universal genes from the metaSNP output

dir = $1  # ex : ~/TARA_metaSNP_light/tara_4_metaSNP/

mkdir $dir/filtered_universal
mkdir $dir/filtered_universal/pop/

cd $dir/filtered_universal/pop/

LIST=$(ls $dir/filtered/pop/ | grep .freq | cut -f1 -d.)

echo -e "\nList of the $(echo "$LIST"| grep -c .) taxonomic ID found :\n$LIST\n\nFor each we are going to extract the universal genes."

for TAX_ID in $LIST;
do
echo -e "\n Organism $TAX_ID \n"
head -n 1 $dir/filtered/pop/"$TAX_ID".filtered.freq > temporaire.1
cut -f2 -d: $dir/filtered/pop/"$TAX_ID".filtered.freq | grep -v TARA_ | grep -v - | cut -f3 -d. | cut -f1 -d_ | sort -u | \
grep -f - /nfs/cds/freeze11/freeze11_mgs/All.ids.txt | cut -f3 -d. | sed 's/$/:/' | \
grep -f -  $dir/filtered/pop/"$TAX_ID".filtered.freq > temporaire.2

cat temporaire.1 temporaire.2 > "$TAX_ID".universal.filtered.freq
rm temporaire*

N_GENES=$(cut -f2 -d: $dir/filtered/pop/"$TAX_ID".filtered.freq | grep -v TARA_ | grep -v - | cut -f3 -d. | cut -f1 -d_ | sort -u | \
grep -f - /nfs/cds/freeze11/freeze11_mgs/All.ids.txt | cut -f3 -d. | sed 's/$/:/' | \
grep -f -  $dir/filtered/pop/"$TAX_ID".filtered.freq | cut -f1 | cut -f2 -d: | sort -u | grep -c .)


COMMON_POS=$(cut -f2 -d: $dir/filtered/pop/"$TAX_ID".filtered.freq | grep -v TARA_ | grep -v - | cut -f3 -d. | cut -f1 -d_ | sort -u | \
grep -f - /nfs/cds/freeze11/freeze11_mgs/All.ids.txt | cut -f3 -d. | sed 's/$/:/' | \
grep -f -  $dir/filtered/pop/"$TAX_ID".filtered.freq | cut -f1 | cut -f2 -d: | grep -c .)



POS=$(cut -f1 $dir/filtered/pop/"$TAX_ID".filtered.freq | cut -f2 -d: | grep -v "^$" | grep -c .)


echo -e "We found $N_GENES genes accounting for $COMMON_POS out of $POS positions.\n \n  **********";

done

