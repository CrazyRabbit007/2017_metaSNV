##############################################################
# Rscript computing genetic diversity of the MetaSNP output. #
##############################################################

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)

###########
# Arguments
# Rscript script "WorkingDirect" "outdir' "species" FALSE NA 
args = commandArgs(trailingOnly=TRUE)
directory = args[1]
outdir = args[2]
species = args[3]
syn = args[4]
###########

###########
# Test parameters
species = ifelse(species=="test","1002672",species)
# Non parallelized parameters :
#files = list.files(paste0(outdir,"filtered/pop/"))
#species = gsub(".filtered.freq","",files)
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
nucleotide.dist <- function(sample1, sample2){
  temp = outer(c(1 - sum(sample1), sample1),c(1 - sum(sample2), sample2))
  diag(temp)<-0
  return(sum(temp))
}
# Compute the pairwise sample distances
pairwise.comp <- function(id,variable,table.species){
  pairwise.NA=unique(c(which(is.na(table.species[,get(variable)])),
                       which(is.na(table.species[,get(id)]))))
  if (length(pairwise.NA)!=0){table.species = table.species[-pairwise.NA]}
  dist.species = table.species[,nucleotide.dist(get(variable),get(id)), by=position]$V1 
  return(sum(dist.species))
}
###########

#for (i in species){ # If not parallelized
i = species
# Loading the data
table.species = fread(paste0(directory,"/filtered/pop/",i,".filtered.freq"))

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

# Subset for Synonymous or Non SNP when asked for
if (!is.na(syn)){
  table.species=table.species[significance==syn]
}

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

system.time({
  if (nrow(gen.dist.melt)==1){
    gen.dist.melt[,value:=pairwise.comp(id,variable,table.species)]
  }else if(nrow(gen.dist.melt)==4){
    # Computing the nucleotide diversity when nrow <=3 (prevent a data.table bug)
    gen.dist.melt[,value:=pairwise.comp(id,variable,table.species),by=row.names(gen.dist.melt)]
  }else{
    # Temporarly subsetting for values != NA
    res = gen.dist.melt[!is.na(value)]
    res = data.table(rows=row.names(res), res)
    # Computing the nucleotide diversity
    res[,value:=pairwise.comp(id,variable,table.species),by=rows]
    # Building the matrix back
    gen.dist.melt[!is.na(value),value:=res[,value]]
  }
})

gen.dist.species=spread(gen.dist.melt,variable,value)
gen.dist.species=as.data.frame(gen.dist.species[,2:length(gen.dist.species)])
gen.dist.species[upper.tri(gen.dist.species)]=NA
# Adding names
row.names(gen.dist.species)=names(gen.dist.species)

# Dividing by minimal covered length of genome
genome.stat = read.delim(paste0(directory,"/statistics/table.SNP.stats.tab"), check.names = F)
coverage.table = read.delim(paste0(directory,"/statistics/all.cov.stat.tab"), check.names = F)
coverage.table$sample=gsub('.*/','',coverage.table$sample)
# Subset for the species of interest
coverage.table = coverage.table[which(coverage.table$species==i),]

# Creating a table to store the min coverage for each paiwise comparison
coverage.correction = gen.dist.species
coverage.correction[,] = NA
  
for (sample1 in names(coverage.correction)){
  pos.sample1 = which(names(coverage.correction)==sample1)
  for (sample2 in names(coverage.correction)[1:pos.sample1]){
    coverage.correction[sample1,sample2]=min(
      coverage.table[which(coverage.table$sample==paste0(sample1,'.cov.summary')),'hcov.1X'],
      coverage.table[which(coverage.table$sample==paste0(sample2,'.cov.summary')),'hcov.1X'])/100
  }
}

gen.dist.species = gen.dist.species / (coverage.correction * genome.stat[which(genome.stat$species==i), "len.genome"])

# Writing the file
syn=ifelse(is.na(syn),"",paste0(".",syn))
write.table(gen.dist.species,paste0(outdir,"/",i,syn,".genetic.distance.tab"),sep="\t")
