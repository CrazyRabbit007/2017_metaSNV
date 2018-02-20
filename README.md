# In-house customized version *under development*

**This repo is based on metaSNV :**
- [MetaSNV Paper](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0182392)
- [MetaSNV Compagnon website](http://metasnv.embl.de/index.html)
- [MetaSNV GitLab](https://git.embl.de/costea/metaSNV) 

## Structure of the repository

1. **Original metaSNV structure** : 
    - README_metaSNV.md
    - metaSNV.py
    - metaSNV_post.py
    - src/
    - db/
    - LICENSE
    - Makefile
    - SETUPFILE
    - EXAMPLE/
    - getRefDB.sh
    
2. **Code edits** : 
    - metaSNV.py :
![alt text](https://github.com/LucasPaoli/2017_metaSNV/blob/master/figures/Code_changes.2.png "metaSNV.py edit")
    - collapse_coverages.py :
![alt text](https://github.com/LucasPaoli/2017_metaSNV/blob/master/figures/Code_changes.1.png "collapse_coverages.py edit")

    
3. **Additional files** : 
    - metaSNV_Filtering_2.0.py : Parallelized filtering step
    - metaSNV_DistDiv.py : Computes pairwise distances, nucleotide diversity and FST
    - motus.remove.padded.sh : Removes padded regions from filtered files
    - match.motu.freeze.ids.sh : Produces a map between mOTUs ids and freeze ids

4. **Running the pipeline** : 
    - Determining parameters :
        - m, b and d : ![alt text](https://github.com/LucasPaoli/2017_metaSNV/blob/master/figures/Tara.motu.coverages.pdf "Example for Tara/mOTUs parameters")
        - c and p : default for c and 90% for p to keep regions mostly shared between all samples.
        - Summary : ![alt text](https://github.com/LucasPaoli/2017_metaSNV/blob/master/figures/Parameters_summary.pdf)


## Example script :

    - [Scripts](https://github.com/LucasPaoli/2017_metaSNV/blob/master/runs) : [Example for HMP/mOTUs](https://github.com/LucasPaoli/2017_metaSNV/blob/master/runs/hmp.motus.sh)

```
###########################################
echo -e "\n\n*************************\n\n"
echo "0. LOADING MODULES"
echo -e "\n\n*************************\n\n"
###########################################

ml HTSlib
ml Boost
ml Python

metaSNV_dir=~/DEV_METASNV/metaSNV

export PATH=$metaSNV_dir:$PATH

######################
# DEFINING VARIABLES #
######################

# Input Files
SAMPLES=../../DATA/tara.75.motu.samples

# Output Directory
OUT=../../DATA/metaSNV_res/tara.75.motu.metasnv # use "output" not "output/"

# DATABASE
# Fasta file
FASTA=/nfs/home/paolil/DEV_METASNV/metaSNV/db/mOTUs/mOTU.v2b.centroids.reformatted.padded
# Genes annotation
GENE_CLEAN=/nfs/home/paolil/DEV_METASNV/metaSNV/db/mOTUs/mOTU.v2b.centroids.reformatted.padded.annotations

# THREADS
threads=16

###########################################
echo -e "\n\n*************************\n\n"
echo "1. COVERAGE COMPUTATION"
echo -e "\n\n*************************\n\n"
###########################################

metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads $threads --n_splits $threads --db_ann "${GENE_CLEAN}" --print-commands > cov.jobs 

# JOB PARRALELLISATION
jnum=$(grep -c "." cov.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 cov.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array

###########################################
echo -e "\n\n*************************\n\n"
echo "2. SNV CALLING"
echo -e "\n\n*************************\n\n"
###########################################

# Repeat command :

metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads $threads --n_splits $threads --db_ann "${GENE_CLEAN}" --print-commands | grep 'samtools mpileup' > snp.jobs

# JOB PARRALELLISATION
jnum=$(grep -c "." snp.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 snp.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array

###########################################
echo -e "\n\n*************************\n\n"
echo "3. POST PROCESSING"
echo -e "\n\n*************************\n\n"
###########################################

# Filtering :
python ~/DEV_METASNV/metaSNV_Filtering_2.0.py "${OUT}" -m 20 -d 10 -b 80 -p 0.9 --n_threads $threads

# Remove Padding :
/nfs/home/paolil/mOTUS_Paper/DATA/motus.remove.padded.sh $OUT/filtered-m20-d10-b80-p0.9/pop

# Compute distances :
python ~/DEV_METASNV/metaSNV_DistDiv.py --filt $OUT/filtered-m20-d10-b80-p0.9 --dist --n_threads $threads
