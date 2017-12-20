#############################################################
# Script generating a summary table for the coverage values #
#############################################################

args = commandArgs(trailingOnly=TRUE)

wdirectory = args[1]
wdirectory=ifelse(is.na(directory),"~/TARA_metaSNP_light/tara_4_metaSNP/statistics",wdirectory)
setwd(wdirectory)

directory = args[2]
directory=ifelse(is.na(directory),"~/TARA_metaSNP_light/tara_4_metaSNP/filtered/pop/",directory)


library(tidyverse)
library(gridExtra)

files = list.files(directory)

names = gsub(".filtered.freq","",files)

all.cov.table = data.frame(station=character(),species = character(),avg.cov=numeric(),hcov.1X=numeric(),hcov.2X=numeric())

for (i in names){

species = i

cov.table = read.delim(paste0("cov/",species,".cov.tab"),header=F)
cov.table[,1]=gsub(".*TARA_|_0.22-.*.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam.cov.summary","",cov.table[,1])

names(cov.table)=c("station","avg.cov","hcov.1X","hcov.2X")
cov.table <- cov.table %>% separate(station, c("station","species"),sep =":")

cov.table = cov.table[cov.table$avg.cov>=10,]
cov.table = cov.table[cov.table$hcov.1X>=40,]

all.cov.table = rbind(all.cov.table,cov.table)
}

write.table(all.cov.table,"all.cov.stat.tab",sep="\t",row.names=F)
