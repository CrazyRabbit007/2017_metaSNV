Custom version of metaSNP, with some changes : metaSNV
============================================== 

The original metaSNP is available [here](http://metasnp.embl.de/index.html)

I. Edits :
----------

1. Correct the allele frequency in the computeDistances.R script. 
```
0.05 istead of 0.5 in line 131
```
2. In metaSNP_distances.sh --> give R as an argument
3. Comment the filtering steps in computeDistances.R / reproducibility with genetic distances
4. Implementation of a genetic distance computation
5. Implementation of additional statistics (relative to Schloissnig 2013)

II. List of additional files :
------------------------------

1. metaSNP_genetic.sh # Calls the following scripts
2. src/computeDiversity.R
3. src/computeFST.R
4. metaSNP_Schloissing # Calls the following scripts
5. src/computeSchloissnig.R
6. src/plotSchloissnig

III. Extra R dependencies :
---------------------------

1. tidyverse
2. data.table

MetaSNP, a metagenomic SNP calling pipeline
============================================== 


The metaSNP pipeline performs variant calling on aligned metagenomic samples and enables species delineation with sub-species resolution.
For an elaborate tutorial visit the [website](http://metasnp.embl.de/index.html).

Download
======== 

Via Git:

    git clone git@git.embl.de:rmuench/metaSNP.git
    
or [download](https://git.embl.de/rmuench/metaSNP/repository/archive.zip?ref=master) a zip file of the repository.


Dependencies
============

* Boost-1.53.0 or above

* samtools-1.19 or above
 
* Python-2.7 or above
    * Pandas
    * Numpy

* R (optional)
    * ape
    * ggplot2
    * gridExtra


Setup & Compilation
===================
MetaSNP is mainly implemented in C/C++ and needs to be compiled.

I. Setup
--------

a)  Change the SAMTOOLS variable in the **SETUPFILE** to your path to samtools.
    If you don't have samtools, you can find it [here](http://samtools.sourceforge.net/):
    
            SAMTOOLS="/path/2/samtools" 
            
a)  Change the HTSLIB variable in the **SETUPFILE** to your path to htslib.
    If you don't have htslib, check your samtools folder first or download it [here](http://www.htslib.org/):
    
            HTSLIB="/path/2/htslib" 

b)  Change the BOOST variable in the **SETUPFILE** to your path to BOOST library.
    If you don't have boost, you can find boost [here](http://www.boost.org/users/download/): 
        
            BOOST="/path/2/boost" 
            
c)  Setup your own **reference database** or acquire the database we use.
    In order to acquire our database run the provided script in the parent directory:
    
            usage: ./getRefDB.sh 


II. Compilation:
----------------

1)   run `make` in the parent directory of metaSNP to compile qaTools and the snpCaller.
        
            usage: make
            

III. Environmental Variables:
----------------
We assume the metaSNP parent directory, ``samtools`` and ``python`` in your system path variable. 
You can set this variable temporarily by runing the following commands and permanently by putting these line into your .bashrc: 
    
            export PATH=/path/2/metaSNP:$PATH
            export PATH=/path/2/python:$PATH
            export PATH=/path/2/samtools:$PATH

Note: Replace '/path/2' with the corresponding global path.

[Tutorial](http://metasnp.embl.de/tutorial.html)
================================================
