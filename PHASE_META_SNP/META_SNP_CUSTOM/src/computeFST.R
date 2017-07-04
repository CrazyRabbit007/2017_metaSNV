################################################
# Rscript computing FST of the MetaSNP output. #
################################################

# Requires the Pi tables (nucleotide diversity within and between samples)

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)

###########
# Arguments
# Rscript script "directory"
args = commandArgs(trailingOnly=TRUE)
directory = args[1]
# directory = "/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/GUT_POPULATION_GENETICS/RESULTS/genetic_final.1/"
###########

###########
# Working directory
setwd(directory)
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

files = list.files("pi/")

species = gsub(".genetic.distance.tab","",files)

for (i in species){
  
  #i = species
  
  # Loading the nucleotide distance
  nucl.dist = read.delim(paste0("pi/",i,".genetic.distance.tab"),row.names=1,check.names = F)
  
  if (length(nucl.dist)>1){
    
    # Creating the distance matrix
    fst.dist = nucl.dist
    fst.dist[lower.tri(fst.dist,diag=T)]=0
    diag(fst.dist)<-NA
    
    # Melting
    fst.dist.m = as.data.table(melt(cbind(fst.dist,id=row.names(fst.dist)),id='id'))
    
    # Computing
    # Temporarly subsetting for values != NA
    res = fst.dist.m[!is.na(value)]
    res = data.table(rows=row.names(res), res)
    res[,value:=get.fst(id,variable,nucl.dist),by=rows]
    
    # Building the matrix back
    fst.dist.m[!is.na(value),value:=res[,value]]
    fst.dist=as.data.frame(spread(fst.dist.m,variable,value))
    row.names(fst.dist)=fst.dist$id
    fst.dist=fst.dist[,-1]
    
    write.table(fst.dist,paste0('FST/',i,".FST.distance.tab"),sep="\t")
  }
}
