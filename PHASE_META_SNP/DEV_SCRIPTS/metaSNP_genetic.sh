#! /bin/bash
############################################
#  metaSNP Step III:   `Post-Processing    #
############################################

# Custom script for the metagenomic SNP calling pipeline (metaSNP)
# Copyright (c) 2017 Lucas Paoli
# Licenced under the GNU General Public License (see LICENSE) 
#
# Helper script for metaSNP genetic distance computation

###########################################
# Working & metaSNP directory directory
metaSNP_dir=$1
OUT=$2
# /!\ Add some fail safe / tests
###########################################
# Script used to compute the diversity
script=$metaSNP_dir/src/computeDiversity.R

# global output directory for the genetic distances
mkdir $OUT/distances_genetic

# store species names
species=$(ls $OUT/filtered/pop/ | cut -f1 -d.)
###########################################


#############################
# Nucleotidic diversity All, Non synonymous, Synonymous
#############################
# output directory
mkdir $OUT/distances_genetic/pi/
out_pi="$OUT/distances_genetic/pi/"

mkdir $OUT/distances_genetic/pi.N/
out_pi_N="$OUT/distances_genetic/pi.N/"


# Create a file to store all the comman lines :
#loop over species
for i in $species;
do
# Call Rscript as : script -outdir -species -[dominant FALSE] -[syn NA]
echo Rscript $script $OUT $out_pi $i 
echo Rscript $script $OUT $out_pi_N $i "N";
done




