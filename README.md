# metaSNV-LP : In-house custom version under development

This project relies on metaSNV :
- Paper : Under revision
- [Compagnon website](http://metasnv.embl.de/index.html)
- [GitLab](https://git.embl.de/costea/metaSNV) : Check for latest update of the core code. Last checked on 2017-07-05. 

## Structure of the project

1. **PHASE_META_SNP** : Archive of the first software draft (early 2017)
2. **metaSNV** : metaSNV code
3. **ANALYSIS**, **DATA**, **RESULTS** : Scripts under work, data to test the and their output
4. **metaSNV_*.py** : Polished and working scripts, to use after the two first steps of the original metaSNV

## Tutorial ([Original metaSNV tutorial](http://metasnv.embl.de/tutorial.html))

### Run the two first steps of metaSNV :

- **Variables :**
````bash
OUT=/project/directory/
SAMPLES=/path/bamfiles_names_list
FASTA=/path/db/database.fasta
GENE_CLEAN=/path/db/annotations
````

- **Coverage estimation :**
````bash
metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads 8 --n_splits 40 --db_ann "${GENE_CLEAN}" --print-commands > cov.jobs
````
Submit the command lines :
````sh
jnum=$(grep -c "." cov.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 cov.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array 
````

- **Variant calling :**
````
metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads 8 --n_splits 40 --db_ann "${GENE_CLEAN}" --print-commands | grep 'samtools mpileup' > snp.jobs
````
Submit the command lines :
````
jnum=$(grep -c "." snp.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 snp.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array
````
