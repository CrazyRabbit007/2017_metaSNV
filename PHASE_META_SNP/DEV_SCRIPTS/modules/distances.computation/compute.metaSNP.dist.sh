#! /bin/bash

#  Script.sh
#
#  Created by Lucas Paoli on 10/02/2017. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#
# Purpose : Calling metaSNP, to be used on the Tara Oceans data.
#	metaSNP can be found here : "git clone git@git.embl.de:rmuench/metaSNP.git" or http://metasnp.embl.de/index.html
#
# Requires Boost-1.53.0 or above, samtools-1.19 or above and Python-2.7 or above. 
#	Make sure to edit the metaSNP "SETUP" file with the path towards those dependencies 
# Optionally R (ape, ggplot2, gridExtra, clue) can be needed for the downstream analysis
#
# This requires metaSNP, python and samtools to be environmental variables. Assuming python and samtools already are :

export PATH=/science/paolil/metaSNP/:$PATH

######################
# DEFINING VARIABLES #
######################

# /!\ TO DO : load variables from $1, $2 etc...
# /!\ TO DO : parallelize the steps that need to

# Input Files
#BAMFILES="/science/paolil/TARA_metaSNP/tara_bamfiles_test/"
#BAMFILES="/science/paolil/metaSNP/EXAMPLE/samples"
BAMFILES=$(for x in $(cat /science/paolil/TARA/tara.all.prok.180); do ls /science/paolil/TARA/$x/reads.filtered.freeze11.marine.repgenomes.solexaqa/$x*.unique.sorted.bam; done)

# Output Directory
#OUT="tara_test_metaSNP" # use "output" not "output/"
OUT="tara_4_metaSNP" # use "output" not "output/"

# DATABASE
# Fasta file
FASTA="/science/paolil/freeze11.marine.repgenomes/freeze11.marine.repgenomes"
#FASTA="/science/paolil/metaSNP/db/RepGenomesv9.fna"
#FASTA="/science/paolil/metaSNP/raw_db/freeze11.genomes.representatives.fa"
# Definition file
DEF="/science/paolil/freeze11.marine.repgenomes/freeze11.marine.repgenomes.len.def.bed"
#DEF="/science/paolil/metaSNP/db/Genomev9_definitions"
#DEF="/science/paolil/metaSNP/raw_db/freeze11.len.def.bed"
	# Note: The genome definition file is in the database you downloaded. If you setup your 
	# own database you need a file with the contig ranges in BED format for each contig. 
	# (Fields: Contig_id, contigStart, contigEnd).
# Genes annotation
#GENE_CLEAN="/science/paolil/metaSNP/db/RefOrganismDB_v9_gene.clean"
#GENE_CLEAN="/science/paolil/freeze_DB_metaSNP/freeze11.annotations.txt"
GENE_CLEAN="/science/paolil/freeze11.marine.repgenomes/freeze11.marine.repgenomes.annotations.txt"



###########################################
echo -e "\n\n*************************\n\n"
echo "3. POST-PROCESSING"
echo -e "\n\n*************************\n\n"
###########################################

echo -e  "\nDISTANCE COMPUTATION\n"

#Use metaSNP_distances.sh (R) or metaSNP_distances.py (python) to compute the pair-wise 
#distances between samples for each species in the filtered discovery set. metaSNP_distances.sh 
#also generates a PCoA projection of the SNP space if the data allows for it.

#python metaSNP_distances.py "${OUT}"/filtered/pop/657317.filtered.freq "${OUT}"/distance.py/

#or

ml R # Load R, does not exist locally

# BUG : Install locally the R libraries : "ape", "clue" and "gridExtra"
# BUG : Change R directory in the file metaSNP_distances.sh

echo "$(metaSNP_distances.sh "${OUT}"/pan_universal_common/ "${OUT}"/distances_pangenome/)"
echo -e "\n\n***\n\n"
eval "$(metaSNP_distances.sh "${OUT}"/pan_universal_common/ "${OUT}"/distances_pangenome/)" # Run a single Job at once
#Note: This script generates command lines for parallel computing.

# JOB PARRALELLISATION
#metaSNP_distances.sh "${OUT}"/filtered/pop/ "${OUT}"/distances/ > dist.jobs # Store the jobs in a file
#jnum=$(grep -c "." dist.jobs) # Store the number of jobs
#/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 dist.jobs # Create a file per job
#qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array



