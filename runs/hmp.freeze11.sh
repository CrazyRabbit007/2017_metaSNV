#! /bin/bash

#  Script.sh
#
#  Created by Lucas Paoli on 15/04/2017. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#
# Purpose : Calling metaSNP/metaSNV, to be used on the Tara Oceans data.
#	metaSNP can be found here : "git clone git@git.embl.de:rmuench/metaSNP.git" or http://metasnp.embl.de/index.html
#
# Requires Boost-1.53.0 or above, samtools-1.19 or above and Python-2.7 or above. 
#	Make sure to edit the metaSNP "SETUP" file with the path towards those dependencies 
# Optionally R (ape, ggplot2, gridExtra, clue) can be needed for the downstream analysis
#
# This requires metaSNP, python and samtools to be environmental variables. Assuming python and samtools already are :


###########################################
echo -e "\n\n*************************\n\n"
echo "0. LOADING MODULES"
echo -e "\n\n*************************\n\n"
###########################################

ml SAMtools
ml HTSlib
ml Boost
ml Python

metaSNV_dir=~/DEV_METASNV/metaSNV

export PATH=$metaSNV_dir:$PATH

######################
# DEFINING VARIABLES #
######################

# Input Files
SAMPLES=../../DATA/hmp.new.freeze.samples

# Output Directory
OUT=../../DATA/metaSNV_res/hmp.new.freeze.metasnv # use "output" not "output/"

# DATABASE
# Fasta file
FASTA=/nfs/home/paolil/DEV_METASNV/metaSNV/db/freeze11/freeze11.genomes.representatives.fa
# Genes annotation
GENE_CLEAN=/nfs/home/paolil/DEV_METASNV/metaSNV/db/freeze11/freeze11.annotations.txt

# THREADS
threads=16

###########################################
echo -e "\n\n*************************\n\n"
echo "1. COVERAGE COMPUTATION"
echo -e "\n\n*************************\n\n"
###########################################

#metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads $threads --n_splits $threads --db_ann "${GENE_CLEAN}" --print-commands > cov.jobs 

# JOB PARRALELLISATION
#jnum=$(grep -c "." cov.jobs) # Store the number of jobs
#/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 cov.jobs # Create a file per job
#qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array 

###########################################
echo -e "\n\n*************************\n\n"
echo "2. SNV CALLING"
echo -e "\n\n*************************\n\n"
###########################################

# Repeat command :

#metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads $threads --n_splits $threads --db_ann "${GENE_CLEAN}" --print-commands | grep 'samtools mpileup' > snp.jobs
#metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads $threads --n_splits $threads --print-commands | grep 'samtools mpileup' > snp.jobs

# JOB PARRALELLISATION
#jnum=$(grep -c "." snp.jobs) # Store the number of jobs
#/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 snp.jobs # Create a file per job
#qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array

###########################################
echo -e "\n\n*************************\n\n"
echo "3. POST PROCESSING"
echo -e "\n\n*************************\n\n"
###########################################

# Filtering :
python ~/DEV_METASNV/metaSNV_Filtering_2.0.py "${OUT}" -m 20 -d 10 -b 40 -p 0.9 --n_threads $threads

# Compute distances :
python ~/DEV_METASNV/metaSNV_DistDiv.py --filt $OUT/filtered-m20-d10-b40-p0.9 --dist --n_threads $threads



