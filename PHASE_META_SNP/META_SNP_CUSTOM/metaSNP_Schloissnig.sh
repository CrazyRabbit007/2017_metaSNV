#! /bin/bash
############################################
#  metaSNP Step III:   `Post-Processing    #
############################################

# Custom script for the metagenomic SNP calling pipeline (metaSNP)
# Copyright (c) 2017 Lucas Paoli
# Licenced under the GNU General Public License (see LICENSE) 
#
# Helper script for metaSNP Schloissnig et al. (2013) analysis

###########################################
# Working & metaSNP directory directory
metaSNP_dir=$1
OUT=$2
DEF=$3
VCOV=$4
HCOV=$5
# /!\ Add some fail safe / tests
###########################################
# Script used to compute the Analysis
script=$metaSNP_dir/src/computeSchloissnig.R

# global output directory for the genetic distances
mkdir $OUT/statistics

# store species names
names=$(ls $OUT/filtered/pop/ | cut -f1 -d.)
###########################################


#############################
# Get coverage values 
#############################
mkdir $OUT/statistics/cov 

for species in $names; do grep $species $OUT/cov/*.summary | sort -k 2n > $OUT/statistics/cov/$species.cov.tab; done

#############################
# Get coverage & SNP statistics 
#############################

echo -e "\nGet coverage & SNP statistics\n"

Rscript $metaSNP_dir/src/computeSchloissnig.R $OUT $DEF $VCOV $HCOV

#############################
# Plot the Schloissnig figure 
#############################

echo -e "\nPlot the Schloissnig figure\n"

Rscript $metaSNP_dir/src/plotSchloissnig.R "$OUT/statistics"
