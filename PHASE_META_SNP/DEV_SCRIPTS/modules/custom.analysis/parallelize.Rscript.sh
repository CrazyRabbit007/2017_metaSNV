#! /bin/bash


# JOB PARRALELLISATION for Rscript

ml R

# Script you want to parralellize :
script=/nfs/home/paolil/TARA_metaSNP_light/custom.analysis/scripts/Schloissnig.implement.enhanced.R

outdir="/nfs/home/paolil/TARA_metaSNP_light/custom.analysis/pi_dominant"

# Create a file to store all the comman lines :
rm R.jobs
touch R.jobs

species=$(ls /nfs/home/paolil/TARA_metaSNP_light/tara_4_metaSNP/filtered/pop/ | cut -f1 -d.)

for i in $species;
do
# Call Rscript as : Rscript -outdir -species -[dominant FALSE] -[syn NA] -[dist ""] -[dist.dir ""]
echo Rscript $script $outdir $i TRUE >> R.jobs; # Store the jobs in a file
done

jnum=$(grep -c "." R.jobs) # Store the number of jobs
/nfs/home/ssunagaw/bork.bin/job.creator.pl 1 R.jobs # Create a file per job
qsub -sync y -V -l h_vmem=8G -t 1-$jnum -pe smp 1 /nfs/home/ssunagaw/bork.bin/run.array.sh # Submit the array

