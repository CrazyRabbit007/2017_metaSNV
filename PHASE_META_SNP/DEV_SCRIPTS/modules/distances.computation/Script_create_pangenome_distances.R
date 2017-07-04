##############################################################
#  metaSNP : Creating pangenome table (universal.common)     #
##############################################################

##########
# Variables :
#args=c("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/4_Run_genes_ann/pangenome/pan.genome.universal.common.filtered.freq",
#""/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/4_Run_genes_ann/pangenome/pan.genome.universal.common.distances.freq)
args <- commandArgs(trailingOnly = TRUE)
input=as.character(args[1])
output=as.character(args[2])
##########

##########
# Libraries :
library(clue)
library(ape)
library(ggplot2)
library(gridExtra)
library(plyr)
library(ggrepel)
library(reshape2)
##########

####################
# I - Loading data #
####################

pangenome.data=read.delim(input,row.names=1)
names(pangenome.data)=gsub("TARA_","",names(pangenome.data))
names(pangenome.data)=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(pangenome.data))
names(pangenome.data)=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(pangenome.data))

#pangenome.data[pangenome.data==-1]=NA

pangenome.dist=as.data.frame(as.matrix(dist(t(pangenome.data),method="manhattan")))/nrow(pangenome.data)*100

################################
# II - Writing pangenome table #
################################

write.table(pangenome.dist,file=output,row.names=F,quote=F,sep="\t")

