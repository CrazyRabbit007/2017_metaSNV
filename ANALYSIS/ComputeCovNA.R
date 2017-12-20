##############################################################
# Rscript computing genetic diversity of the MetaSNP output. #
##############################################################

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)

###########
# Arguments
# Rscript script "WorkingDirect" "outdir"
args = commandArgs(trailingOnly=TRUE)
directory = args[1]
outdir = args[2]
###########

###########
# Test parameters
files = list.files(directory)
###########

###########
# Libraries
library(tidyverse)
library(data.table)
###########

###########
# Functions
fast.replace = function(DT,x=-1,y=NA) {
  # by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==x),j,y)
}
# Compute the pairwise sample distances
pairwise.comp <- function(id,variable,table.species){
  
  id.NA = which(is.na(table.species[,get(id)]))
  variable.NA = which(is.na(table.species[,get(variable)]))
  
  pairwise.NA=as.vector(na.omit(match(variable.NA,id.NA)))
  
  id.unique = id.NA[!id.NA%in%pairwise.NA]
  variable.unique=variable.NA[!variable.NA%in%pairwise.NA]
  
  return(min(length(id.unique),length(variable.unique))/min(length(na.omit(table.species[,get(id)])),length(na.omit(table.species[,get(variable)]))))
}
###########

for (i in files){
  
  print(paste0(directory,"/",i))
  
  # Loading the data
  table.species = fread(paste0(directory,"/",i))
  
  # Replacing '-1' by 'NA'
  fast.replace(table.species)
  
  # Formatting the first column :
  # with annotation : taxonomy:gene:position:SNP:codon
  # without annotation : taxonomy:-:position:SNP:.
  table.species[,V1:=gsub(":\\.",":[-]",table.species[,V1])]
  # Splitting the first column into each variable
  table.species <- table.species %>% separate(V1, c("taxo.ann","gene.ann","position","SNP","significance"),sep =":")
  # Splitting the codon information
  table.species <- table.species %>% separate(significance, c("significance","codon"),sep ="\\[")
  table.species[,codon:=gsub("\\]","",table.species[,codon])]
  
  # When the reference genome is in several chromosomes, you need to add the gene name to have unique positions per nucleotide.
  table.species[,position:=paste0(table.species$taxo.ann,":",table.species$gene.ann,":",table.species$position)]
  
  # Setting the position as a key (~ row.names for data.table)
  setkey(table.species,position)
  
  
  # Storing the number of samples (column 8 onward) and their name
  n = as.numeric(length(table.species[,7:length(table.species)]))
  samples = as.character(names(table.species)[7:length(table.species)])
  
  # creating the distance matrix
  gen.dist.species=as.data.frame(matrix(NA,ncol=n,nrow=n),row.names=samples)
  names(gen.dist.species)=row.names(gen.dist.species)
  gen.dist.species[lower.tri(gen.dist.species,diag=T)]=0
  
  # Melting the distance matrix
  gen.dist.melt = as.data.table(melt(cbind(gen.dist.species,id=row.names(gen.dist.species),stringsAsFactors=F),id="id",variable.factor=FALSE,value.factor=FALSE))
  gen.dist.melt[,variable:=as.character(variable)]
  
  if (nrow(gen.dist.melt)==1){
    gen.dist.melt[,value:=pairwise.comp(id,variable,table.species)]
  }else if(nrow(gen.dist.melt)==4){
    # Computing the nucleotide diversity when nrow <=3 (prevent a data.table bug)
    gen.dist.melt[,value:=pairwise.comp(id,variable,table.species),by=row.names(gen.dist.melt)]
  }else{
    # Temporarly subsetting for values != NA
    res = gen.dist.melt[!is.na(value)]
    res = data.table(rows=row.names(res), res)
    res[,value:=pairwise.comp(id,variable,table.species),by=rows]
    # Building the matrix back
    gen.dist.melt[!is.na(value),value:=res[,value]]
  }
  
  gen.dist.species=spread(gen.dist.melt,variable,value)
  gen.dist.species=as.data.frame(gen.dist.species[,2:length(gen.dist.species)])
  gen.dist.species[upper.tri(gen.dist.species)]=NA
  # Adding names
  row.names(gen.dist.species)=names(gen.dist.species)
  
  # Writing the file
  write.table(gen.dist.species,paste0(outdir,"/",i,".covest.tab"),sep="\t")
}