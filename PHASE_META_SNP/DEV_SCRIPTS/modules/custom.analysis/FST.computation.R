################################################
# Rscript computing FST of the MetaSNP output. #
################################################

# Requires the Pi tables (nucleotide diversity within and between samples)

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)

###########
# Arguments
# Rscript script "pi_directory" "output"
args = commandArgs(trailingOnly=TRUE)
pi_directory = args[1]
out_dir = args[2]
###########

###########
# Test parameters
species = "1002672"
###########

###########
# Working directory
#setwd(out_dir)
setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/8_genetic_distances/")
###########

###########
# Libraries
library(tidyverse)
library(data.table)
###########

###########
# Function
get.fst <- function(sample1,sample2,nucl.dist){
  sample1=as.character(sample1)
  sample2=as.character(sample2)
  pi.within = (nucl.dist[sample1,sample1]+nucl.dist[sample2,sample2])/2
  fst = 1-(pi.within/nucl.dist[sample1,sample2])
  return(fst)
}
###########

#files = list.files("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered/pop/")
files = list.files("pi/")

species = gsub(".genetic.distance.tab","",files)
for (i in species){
  
  #i = species
  
  # Loading the nucleotide distance
  nucl.dist = read.delim(paste0("pi/",i,".genetic.distance.tab"),row.names=1,check.names = F)
  
  # Creating the distance matrix
  fst.dist = nucl.dist
  fst.dist[lower.tri(fst.dist,diag=T)]=0
  
  # Melting
  fst.dist.m = as.data.table(melt(cbind(fst.dist,id=row.names(fst.dist)),id='id'))
  
  # Computing
  fst.dist.m[,value:=get.fst(id,variable,nucl.dist),by=row.names(fst.dist.m)]
  
  # Building back the matrix
  fst.dist=as.data.frame(spread(fst.dist.m,variable,value))
  row.names(fst.dist)=fst.dist$id
  fst.dist=fst.dist[,-1]
  fst.dist[upper.tri(fst.dist,diag=T)]=NA
  
  write.table(fst.dist,paste0('fst/',i,".FST.distance.tab"),sep="\t")
  
}
