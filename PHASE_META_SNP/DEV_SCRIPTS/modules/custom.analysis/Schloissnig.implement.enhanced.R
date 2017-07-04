################################################################
# Rscript computing genetique diversity of the MetaSNP output. #
################################################################

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)

###########
# Arguments
# Rscript script "1002672" FALSE "dist" "dist.dir" 
args = commandArgs(trailingOnly=TRUE)
outdir = args[1]
species = args[2]
dominant = args[3]
syn = args[4]
dist = args[5]
dist.dir = args[6]
# Default values for non required parameters
dominant=ifelse(is.na(dominant),FALSE,dominant)
dist=ifelse(is.na(dist),"",dist)
dist.dir=ifelse(is.na(dist.dir),"",dist.dir)
###########

###########
# Test parameters
species = ifelse(species=="test","1002672",species)
#species = "1002672"
#dist = ".core"
#dist.dir = "_universal"
# Non parallelized parameters :
#files = list.files("/nfs/home/paolil/TARA_metaSNP_light/tara_4_metaSNP/filtered/pop/")
#species = gsub(".filtered.freq","",files)
###########

###########
# Working directory
setwd(outdir)
#setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/8_genetic_distances/")
###########

###########
# Libraries
library(tidyverse)
library(data.table)
###########

###########
# Functions
# Fast replacement for data.table :
fast.replace = function(DT,x=-1,y=NA) {
  # by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==x),j,y)
}
# Computing positional nucleotide diversity
nucleotide.dist <- function(sample1,sample2,dominant){
  temp = cbind(sample1=c(1 - sum(sample1),sample1),
               sample2=c(1 - sum(sample2),sample2))
  
  if (dominant){
    temp[temp<0.2]=0
    vect = 1 - colSums(temp)
    vect = vect / colSums(temp!=0)
    temp[,1] = temp[,1] + vect[1]*(temp!=0)[,1]
    temp[,2] = temp[,2] + vect[2]*(temp!=0)[,2]
  }
  
  sum.k = 0
  for (row in 1:nrow(temp)){
    sum.k=sum.k+sum(temp[row,"sample1"]*temp[-row,"sample2"])
  }
  return(sum.k)
}
# Compute the pairwise sample distances
pairwise.comp <- function(id,variable,table.species,dominant){
  sample1 = as.name(variable)
  sample2 = as.name(id)
  dist.species = table.species[,nucleotide.dist(eval(sample1),eval(sample2),dominant), by=position]$V1
  return(sum(dist.species,na.rm=T)/sum(!is.na(dist.species)))
  #return(dist.species)
}
###########

#for (i in species){ # If not parallelized
i = species
# Loading the data
table.species = fread(paste0("/nfs/home/paolil/TARA_metaSNP_light/tara_4_metaSNP/filtered",dist.dir,"/pop/",i,dist,".filtered.freq"))
#table.species = fread("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/4_Run_genes_ann/universal_common/filtered/1002672.universal.common.filtered.freq")

# Formatting the names, replacing '-1' by 'NA'
names(table.species)=gsub("TARA_|_0.22-.*_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam","",names(table.species))
fast.replace(table.species)

# Formatting the first column :
# with annotation : taxonomy:gene:position:SNP:codon
# without annotation : taxonomu:-:position:SNP:.
table.species[,V1:=gsub(":\\.",":[-]",table.species[,V1])] # Check if it works
# Splitting the first column into each variable
table.species <- table.species %>% separate(V1, c("taxo.ann","gene.ann","position","SNP","significance"),sep =":")
# Splitting the codon information
table.species <- table.species %>% separate(significance, c("significance","codon"),sep ="\\[") # Can it work with NA ?
table.species[,codon:=gsub("\\]","",table.species[,codon])]

# When the reference genome is in several chromosomes, you need to add the gene name to have unique positions per nucleotide.
table.species[,position:=paste0(table.species$gene.ann,":",table.species$position)]

# Setting the position as a key (~ row.names for data.table)
setkey(table.species,position)

# Subset for Synonymous or Non SNP
if (!is.na(syn)){
table.species=table.species[significance==syn]
}

# Storing the number of samples (column 7 onward) and their name
n = as.numeric(length(table.species[,7:length(table.species)]))
samples = as.character(names(table.species)[7:length(table.species)])

# creating the distance matrix
gen.dist.species=as.data.frame(matrix(NA,ncol=n,nrow=n),row.names=samples)
names(gen.dist.species)=row.names(gen.dist.species)
gen.dist.species[lower.tri(gen.dist.species,diag=T)]=0

# Melting the distance matrix
gen.dist.melt = as.data.table(melt(cbind(gen.dist.species,id=row.names(gen.dist.species),stringsAsFactors=F),id="id",variable.factor=FALSE,value.factor=FALSE))
gen.dist.melt[,variable:=as.character(variable)]

# Temporarly subsetting for vallu != NA
res = gen.dist.melt[!is.na(value)]

# Computing the nucleotide diversity
res[,value:=pairwise.comp(id,variable,table.species,dominant),by=row.names(res)]

# Building the matrix back
gen.dist.melt[!is.na(value),value:=res[,value]]
gen.dist.species=spread(gen.dist.melt,variable,value)
gen.dist.species=as.data.frame(gen.dist.species[,2:length(gen.dist.species)])
# Adding names
row.names(gen.dist.species)=names(gen.dist.species)

# Writing the file
dist=ifelse(dominant,paste0(dist,".dominant"),dist)
write.table(gen.dist.species,paste0(i,dist,syn,".genetic.distance.tab"),sep="\t")


# PROOF TEST
# gen.dist.fast = read.delim("1002672.gen.dist.matrix.fast.tab")
# gen.dist.species = read.delim("1002672.gen.dist.matrix.tab")
# all.equal(gen.dist.fast,gen.dist.species) # Ok

# ##### Alternative options :
# ## 1
# system.time({
# for (i in 1:nrow(res)){
#   sample1 = as.name(res[i,variable])
#   sample2 = as.name(res[i,id])
#   dist.species = table.species[,nucleotide.dist(eval(sample1),eval(sample2)), by=position]$V1
#   set(res,i,which(names(res)=="value"),(sum(dist.species,na.rm=T)/sum(!is.na(dist.species))))
# }
# })
# 
# ## 2
# system.time({
#   for (i in 1:nrow(res)){
#     sample1 = as.name(res[i,variable])
#     sample2 = as.name(res[i,id])
#     dist.species = table.species[,nucleotide.dist(eval(sample1),eval(sample2)), by=position]$V1
#     res[i,value:=(sum(dist.species,na.rm=T)/sum(!is.na(dist.species)))]
#   }
# })

