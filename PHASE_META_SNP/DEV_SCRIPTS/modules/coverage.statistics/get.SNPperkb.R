###########################################################
# Script generating a summary table for SNP per kb values #
###########################################################

args = commandArgs(trailingOnly=TRUE)

wdirectory = args[1]
wdirectory=ifelse(is.na(directory),"~/TARA_metaSNP_light/tara_4_metaSNP/statistics",wdirectory)
setwd(wdirectory) 

directory = args[2] 
directory=ifelse(is.na(directory),"~/TARA_metaSNP_light/tara_4_metaSNP/filtered/pop/",directory)  

bedfile = args[3]
bedfile = ifelse(is.na(bedfile),"/science/paolil/freeze11.marine.repgenomes/freeze11.marine.repgenomes.len.def.bed",bedfile)


library(tidyverse)
library(gridExtra)
library(data.table)

files = list.files(directory)

names = gsub(".filtered.freq","",files)


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
ref.genomes <- ref.genomes %>% separate(Taxo, c("Species","TaxID","Contig"),sep ="\\.")
ref.genomes[, ref.length:=sum(Length),by=Species]
ref.genomes = ref.genomes[!duplicated(ref.genomes$Species)]


print("loop over Species")
for (i in names){

print(i)

table.species = fread(paste0(directory,i,".filtered.freq"))

# Formatting the names, replacing '-1' by 'NA'
names(table.species)=gsub("TARA_|_0.22-.*_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(table.species))

# Formatting the first column :
# with annotation : taxonomy:gene:position:SNP:codon
# without annotation : taxonomu:-:position:SNP:.
table.species[,V1:=gsub(":\\.",":-[-]",table.species[,V1])] # Check if it works
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

write.table(kb.SNP.table,"table.SNP.stats.tab",sep="\t")
