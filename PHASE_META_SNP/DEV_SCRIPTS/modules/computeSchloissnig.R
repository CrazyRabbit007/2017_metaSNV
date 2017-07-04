#############################################################
# Script generating a summary table for the coverage values #
#############################################################

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)


################
# Arguments
args = commandArgs(trailingOnly=TRUE)
wdirectory = args[1]
bedfile = args[2]
vcov.filter = args[3]
hcov.filter = args[4]
setwd(wdirectory)
################

################
# Libraries
library(tidyverse)
library(data.table)
################

# List species
files = list.files("filtered/pop")
names = gsub(".filtered.freq","",files)

############
# COVERAGE #
############

# Initialize table
all.cov.table = data.frame(sample=character(),
                        species = character(),
                        avg.cov=numeric(),
                        hcov.1X=numeric(),
                        hcov.2X=numeric())

# Loop over species
for (species in names){

	# Read the custom coverage file
	cov.table = read.delim(paste0("statistics/cov/",species,".cov.tab"),header=F)
	# Edit the first column name to be human readable
	cov.table[,1]=gsub(paste0(wdirectory,"cov\\/"),"",cov.table[,1])
	
	# Match names of the species table with summary table
	names(cov.table)=c("sample","avg.cov","hcov.1X","hcov.2X")
	# Seperate first column in sample + species
	cov.table <- cov.table %>% separate(sample, c("sample","species"),sep =":")

	# Select for rows matching avg and horizontal coverage filters
	cov.table = cov.table[cov.table$avg.cov>=vcov.filter,]
	cov.table = cov.table[cov.table$hcov.1X>=hcov.filter,]
	
	# Storing speceis data in summary table
	all.cov.table = rbind(all.cov.table,cov.table)
}

# Writing table
write.table(all.cov.table,"statistics/all.cov.stat.tab",sep="\t",row.names=F)

########
# SNPs #
########

kb.SNP.table = data.frame(species = character(),
                        total.SNP=numeric(),
                        total.multiallelic=numeric(),
                        N.SNP=numeric(),
                        N.multiallelic=numeric(),
                        S.SNP=numeric(),
                        S.multiallelic=numeric(),
                        len.genome=numeric(),
                        SNP.per.Kb=numeric())


# Loading bedfile describing reference genomes
ref.genomes = fread(bedfile)
names(ref.genomes) = c("Taxo","Start","Length")

print(as.character(ref.genomes[1,"Taxo"]))
test.name = strsplit(as.character(ref.genomes[1,"Taxo"]),"\\.")[[1]]
print(test.name)

if (length(test.name)==3) {
	ref.genomes <- ref.genomes %>% separate(Taxo, c("Species","TaxID","Contig"),sep ="\\.")
} else if (length(test.name)==2) {
	ref.genomes <- ref.genomes %>% separate(Taxo, c("Species","Contig"),sep ="\\.")
} else {
	print("Mmmmh. The script can't deal with your database naming convention.")
}

ref.genomes[, ref.length:=sum(Length),by=Species]
ref.genomes = ref.genomes[!duplicated(ref.genomes$Species)]


print("loop over Species")
for (i in names){

	print(i)

	table.species = fread(paste0("filtered/pop/",i,".filtered.freq"))

	# Formatting the first column :
	# with annotation : taxonomy:gene:position:SNP:codon
	# without annotation : taxonomy:-:position:SNP:.
	table.species[,V1:=gsub(":\\.",":-[-]",table.species[,V1])] 
	# Splitting the first column into each variable
	table.species <- table.species %>% separate(V1, c("taxo.ann","gene.ann","position","SNP","significance"),sep =":")
	# Splitting the codon information
	table.species <- table.species %>% separate(significance, c("significance","codon"),sep ="\\[") # Can it work with NA ?
	table.species[,codon:=gsub("\\]","",table.species[,codon])]

	# When the reference genome is in several chromosomes, you need to add the gene name to have unique positions per nucleotide.
	table.species[,position:=paste0(table.species$gene.ann,":",table.species$position)]

	# Subsetting for significance
	table.species.S=table.species[significance=="S"]
	table.species.N=table.species[significance=="N"]

	temp=data.frame(species=i,
                    total.SNP=length(unique(table.species[,position])),
                    total.multiallelic=(length(table.species[,position])/length(unique(table.species[,position]))),
                    N.SNP=length(unique(table.species.N$position)),
                    N.multiallelic=(length(table.species.N$position)/length(unique(table.species.N$position))),
                    S.SNP=length(unique(table.species.S$position)),
                    S.multiallelic=(length(table.species.S$position)/length(unique(table.species.S$position))),
                    len.genome=ref.genomes[which(ref.genomes$Species==i),ref.length])

	kb.SNP.table = rbind(kb.SNP.table,temp)
}

kb.SNP.table$SNP.per.Kb=kb.SNP.table$total.SNP/(kb.SNP.table$len.genome/1000)

write.table(kb.SNP.table,"statistics/table.SNP.stats.tab",sep="\t")
