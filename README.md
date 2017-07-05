# metaSNV-LP : In-house custom version under development

**This project relies on metaSNV :**
- Paper : Under revision
- [Compagnon website](http://metasnv.embl.de/index.html)
- [GitLab](https://git.embl.de/costea/metaSNV) : Check for latest update of the core code. Last checked on 2017-07-05. 

**Purpose :**

Identify variable genomic positions from aligned metagenomic data (bamfiles), describe the populations (samples) and compute genetic distances between said populations.

**Input :**
- Bamfiles
- Database
    - Reference genomes (fasta)
    - Gene annotations (bed)

## Structure of the repository

1. **PHASE_META_SNP** : Archive of the first software draft (early 2017)
2. **metaSNV** : metaSNV code
3. **ANALYSIS**, **DATA**, **RESULTS** : Scripts under work, data to test them and their output
4. **metaSNV_*.py** : Polished and working scripts, to use after the two first steps of the original metaSNV
    - metaSNV_filtering.py : Filtering step
    - metaSNV_join.py : Join datasets (needed when working both on MetaG and MetaT)
    - metaSNV_universal.py : Extract universal genes
    - metaSNV_stats.py : Computes descriptive statistics
    - metaSNV_DistDiv.py : Computes pairwise distances, diversity and FST
    - metaSNV_pnps.py : Computes pnps per genome and per gene

## Tutorial ([Original metaSNV tutorial](http://metasnv.embl.de/tutorial.html))

- **Variables :**

`````bash
OUT=/project/directory/
SAMPLES=/path/bamfiles_names_list
FASTA=/path/db/database.fasta
GENE_CLEAN=/path/db/annotations
`````

### 1. Run the two first steps of metaSNV :

- **Coverage estimation :**

````bash
metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads 8 --n_splits 40 --db_ann "${GENE_CLEAN}" --print-commands > cov.jobs
````
Submit the command lines :
````bash
jnum=$(grep -c "." cov.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 cov.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array 
````

- **Variant calling :**

````bash
metaSNV.py "${OUT}" "${SAMPLES}" "${FASTA}" --threads 8 --n_splits 40 --db_ann "${GENE_CLEAN}" --print-commands | grep 'samtools mpileup' > snp.jobs
````
Submit the command lines :
````bash
jnum=$(grep -c "." snp.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 snp.jobs # Create a file per job
qsub -sync y -V -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array
````

### 2. Filtering :
````bash
python metaSNV_filtering.py --help
````

### 3. [Optional] Joining datasets :
When working on MetaG and MetaT
````bash
python metaSNV_join.py --help
````

### 4. [Optional] Extract universal genes only :
````bash
python metaSNV_universal.py --help
````

### 5. Compute descriptive statistics :
````bash
python metaSNV_stats.py --help
````

### 6. Distance computation :
````bash
python metaSNV_DistDiv.py --help
````

### 7. Compute pnps values :
````bash
python metaSNV_pnps.py --help
````

