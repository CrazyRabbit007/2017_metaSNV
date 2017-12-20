#!/bin/bash
filtered_freeze=/nfs/home/paolil/mOTUS_Paper/DATA/metaSNV_res/hmp.75.freeze11.metasnv/filtered-m20-d10-b40-p0.9/pop/
filtered_mOTUs=/nfs/home/paolil/mOTUS_Paper/DATA/metaSNV_res/hmp.75.motu.metasnv/filtered-m20-d10-b80-p0.9/pop/
map_file=hmp.75.m20.map.taxID.txt
master_map=/nfs/cds/mOTUv2.20170607/mOTUv2_db/padded/aug.centroids/mOTU.v2b.master.map.tsv
ID_map=/nfs/cds/mOTUv2.20170607/mOTUv2_db/padded/mOTULG_specI.centroids/Gene2mOTULG_specIv2.ID.map
linkage_file=/nfs/cds/mOTUv2.20170607/analysis/data/mOTU_linkage_groups+extended-specI__PPV-0.80.tsv
taxonomy=/science/paolil/DATABASES/metaSNV_freeze/freeze11/freeze11.genomes.taxonomy 

echo -e "\n***************\n* From freeze *\n***************"
freeze_sp=$(ls $filtered_freeze | cut -f1 -d.)

specI_freeze=$(for i in $(echo $freeze_sp); do grep $'\t'$i'\.' $master_map; done | sort -k4 -u | cut -f5 | sort -u)
ids_freeze=$(for i in $(echo $specI_freeze); do grep $i$'\t' $ID_map; done | sort -u -k2 | cut -f3)

num=$(for i in $(echo $ids_freeze); do ls $filtered_mOTUs | grep $i.filt; done | cut -f1 -d.)
num_s=$(for i in $(echo "$num"); do grep $'\t'$i$ $ID_map; done | cut -f2 | sort -u)
num_f=$(for i in $(echo "$num_s"); do grep $i$ $master_map; done | cut -f4 | cut -f1 -d. | sort -u)
num_freeze=$(for i in $(echo "$num_f"); do ls $filtered_freeze | grep $i.filt; done)

missing=$(comm -3 <(echo "$ids_freeze") <(echo "$num"))
missing_s=$(for i in $missing; do grep $'\t'$i$ $ID_map;done | sort -k2 -u | cut -f2)
missing_f=$(grep $missing_s$ $master_map | cut -f4 | sort -u | cut -f1 -d.)

echo 'Found '$(echo "$freeze_sp" | wc -l)' freeze genomes in metaSNV of which '$(echo "$specI_freeze" | wc -l)' have a SpecI in the master map.'
echo $(echo "$ids_freeze" | wc -l)' of the specI have an ID from the ID map (if some are missing, this is likely due to having less than 5 genes among the 10 selected for mOTUs).'
echo $(echo "$num" | wc -l)' were detected in mOTUs by metaSNV.'
#echo -e 'The missing id :\n'"$missing"'\nbeing the representative genomes :\n'"$missing_f"

echo -e "\n**************\n* From mOTUs *\n**************"
ids=$(ls $filtered_mOTUs | cut -f1 -d. | sort -n)

specI=$(for i in $(echo $ids); do grep $'\t'$i\$ $ID_map; done | sort -k2 -u | cut -f2)
specI_found=$(for i in $(echo $specI); do grep $i\$ $master_map; done | sort -k5 -u | wc -l)
taxa=$(for i in $(echo $specI); do grep $i\$ $master_map; done | sort -k4 -u | cut -f4 | cut -f1 -d. | sort -u)

freeze=$(for i in $(echo $taxa); do ls $filtered_freeze | grep $i.filt; done)

echo 'Found '$(echo "$specI" | wc -l)' mOTUs in metaSNV of which '$(echo "$specI" | grep spec | wc -l)' are specI and '"$specI_found"' have a taxonomic assignments ('$(echo "$taxa" | wc -l)' taxa in total) in the master map.'
echo $(echo "$freeze" | wc -l)' were found in freeze by metaSNV.'

echo -e "\nTesting..."
test "$freeze" = "$num_freeze" && echo 'The '$(echo "$freeze" | wc -l)' reference genomes found are indeed the same.'

for i in $(echo $ids); do echo $i$'\t'$(echo $(grep $'\t'$i\$ $ID_map | sort -k2 -u | cut -f2)$ | grep -f - $master_map | sort -k4 -u | cut -f4 | cut -f1 -d'.' | grep -f - <(ls $filtered_freeze) | cut -f1 -d.) ; done > $map_file

#sleep 3

#echo -e "\n**********************\n* Taxonomy of the matching  *\n**********************"
#for i in $(echo "$freeze" | cut -f1 -d.); do grep ^$i$'\t' $taxonomy; done | sort -k2

#echo -e "\n*****************\n* Num of seq ? *\n*****************"
#for i in $(echo "$freeze" | cut -f1 -d.); do echo;echo $i : $(grep $'\t'$i $master_map | wc -l);grep $'\t'$i $master_map; done

#sleep 3

#echo -e "\n************************************\n* Taxonomy of the freeze missing *\n************************************"
#for i in $(comm -3 <(echo "$freeze_sp") <(echo "$freeze" | cut -f1 -d.)); do grep ^$i$'\t' $taxonomy; done | sort -k2

#echo -e "\n****************\n* Num of seq ? *\n****************"
#for i in $(comm -3 <(echo "$freeze_sp") <(echo "$freeze" | cut -f1 -d.)); do echo;echo $i : $(grep $'\t'$i $master_map | wc -l);grep $'\t'$i $master_map; done

#echo -e "\n******************\n* linkage file ? *\n******************"
#for i in $specI_freeze; do grep $i$ $linkage_file; done
#for i in $specI_freeze; do grep $i$ $linkage_file; done | sort -k2 -u | wc -l

#sleep 3
# Not working :
#echo -e "\n********************************\n* Taxonomy of the new 3 taxIDs *\n********************************"
#for i in $(comm -1 <(echo "$freeze_sp") <(echo "$freeze" | cut -f1 -d.)); do grep ^$i$'\t' $taxonomy; done | sort -k2









