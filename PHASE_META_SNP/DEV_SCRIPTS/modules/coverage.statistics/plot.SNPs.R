############################################################################
# This script generates an histogram of the SNPs frequency for one species #
############################################################################

args = commandArgs(trailingOnly=TRUE)

species = args[1]
dist = args[2]
dir = args [3]

setwd("/science/paolil/TARA_metaSNP/tara_coverage/")

library(ggplot2)
library(reshape2)
library(gridExtra)

#species = "1096769"
#dist = ".core"
#dir = "_universal"


files = list.files("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered/pop/")

names = gsub(".filtered.freq","",files)

for (i in names){

species = i

SNP.table = read.delim(paste0("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered_universal/pop/",species,".core.filtered.freq"))
names(SNP.table)=gsub("TARA_","",names(SNP.table))
names(SNP.table)=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table))
names(SNP.table)=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table))
SNP.table=SNP.table[,2:length(SNP.table)]

SNP.table.all = read.delim(paste0("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered/pop/",species,".filtered.freq"))
names(SNP.table.all)=gsub("TARA_","",names(SNP.table.all))
names(SNP.table.all)=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table.all))
names(SNP.table.all)=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table.all))
SNP.table.all=SNP.table.all[,2:length(SNP.table.all)]

freq.table = t(as.data.frame(colSums(SNP.table==-1)/nrow(SNP.table)))

write.table(freq.table, file = paste0("minus.1.freq/",species,".freq.tab"), sep = "\t",row.names=F )

freq.table.m = melt(freq.table)

freq.table.all = t(as.data.frame(colSums(SNP.table.all==-1)/nrow(SNP.table.all)))

write.table(freq.table.all, file = paste0("minus.1.freq/",species,".all.freq.tab"), sep = "\t",row.names=F)

freq.table.all.m = melt(freq.table.all)

pdf(paste0("plots/Freq.",species,".pdf"))
grid.arrange(
ggplot(freq.table.m) + geom_bar(aes(x=row.names(freq.table.m),y=value),stat="identity") + theme(axis.text.x=element_text(angle = 90,size=6)) +
        ylab("-1 Frequency") + xlab("stations") + ggtitle(paste0("Universal genes")), 
ggplot(freq.table.all.m) + geom_bar(aes(x=row.names(freq.table.all.m),y=value),stat="identity") + theme(axis.text.x=element_text(angle = 90,size=6)) +
	ylab("-1 Frequency") + xlab("stations") + ggtitle(paste0("Full genome")),
ncol=1)
dev.off()

SNP.table.m = melt(cbind(SNP.table,id=as.numeric(row.names(SNP.table))),id="id")

pdf(paste0("plots/Hist.",species,".universal.pdf"), width = 10, height = length(SNP.table))
grid.arrange(ggplot(SNP.table.m) + geom_bar(aes(x=id,y=value),stat="identity",width=1) + facet_wrap(~variable,ncol=1))
dev.off()
}

species = args[1]

if (species != FALSE){
SNP.table = read.delim(paste0("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered",dir,"/pop/",species,dist,".filtered.freq"))
names(SNP.table)=gsub("TARA_","",names(SNP.table))
names(SNP.table)=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table))
names(SNP.table)=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(SNP.table))
SNP.table=SNP.table[,2:length(SNP.table)]

SNP.table.m = melt(cbind(SNP.table,id=as.numeric(row.names(SNP.table))),id="id")

pdf(paste0("plots/Hist.",species,dist,".pdf"), width = 10, height = length(SNP.table))
grid.arrange(ggplot(SNP.table.m) + geom_bar(aes(x=id,y=value),stat="identity",width=1) + facet_wrap(~variable,ncol=1))
dev.off()}






