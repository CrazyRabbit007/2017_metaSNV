#! /bin/bash

#  Script.sh
#
#  Created by Lucas Paoli on 10/02/2017. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#
#  Extra coverage values
# 

# Arg : directory (ex : ~/TARA_metaSNP_light/tara_4_metaSNP/
dir = $1

mkdir $dir/statistics
mkdir $dir/statistics/cov

cd $dir/statistics/cov

names=$(ls $dir/filtered/pop | cut -f1 -d.) 
for species in $names; do grep $species $dir/cov/*.summary | sort -k 2n > $species.cov.tab; done
for species in $names; do grep $species $dir/cov/*.detail | sort -k 2n > $species.cov.detail.tab; done
