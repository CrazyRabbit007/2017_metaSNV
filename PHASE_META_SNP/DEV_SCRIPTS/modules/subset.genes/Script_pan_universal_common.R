#############################################################
#  metaSNP : Step III   `Post-processing` - FREQUENCIES     #
#############################################################

# This code is part of the metagenomic SNP calling pipeline (metaSNP)

setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/")


library(clue)
library(ape)
library(ggplot2)
library(gridExtra)
library(plyr)
library(ggrepel)
library(reshape2)

############## PART 2 : universal_common ############

setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/4_distances/universal_common/")
getwd()

files = list.files("filtered/")

names = gsub(".universal.common.filtered.freq","",files)

stations = sort(unique(read.delim("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/3_plot_filter/proj/all_mann_pcoa_proj.tab")[,1]))
stations=gsub("TARA_","",stations)
stations=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",stations)
stations=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",stations)
stations = stations[grepl("_SRF", stations)]

pan.genome = data.frame(matrix(-1,0,length(stations)+1))
names(pan.genome) = c("X",as.character(stations))

for (i in names){
  temp = read.delim(paste0("filtered/",i,".universal.common.filtered.freq"))
  names(temp)=gsub("TARA_","",names(temp))
  names(temp)=gsub("_0.22.1.6_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(temp))
  names(temp)=gsub("_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(temp))
  temp = cbind(X=temp[,"X"],temp[,grepl("_SRF",names(temp))])
  
  if (nrow(temp)>0){

    temp.bind=as.data.frame(matrix(-1,nrow(temp),length(pan.genome)))
    names(temp.bind)=names(pan.genome)

    rows=nrow(pan.genome)
    pan.genome = rbind(pan.genome,temp.bind)
    for (j in names(temp)){
        pan.genome[(rows+1):(rows+nrow(temp)),j] = as.character(temp[,j])
    }
  }
}


