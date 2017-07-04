# Rscript implementing the analysis presented in Schloissnig et al 2013

args = commandArgs(trailingOnly=TRUE)

species = args[1]
dist = args[2]
dist.dir = args[3]
dist=ifelse(is.na(dist),"",dist)
dist.dir=ifelse(is.na(dist.dir),"",dist.dir)

library(ggplot2)
library(reshape2)
library(gridExtra)
library(splitstackshape)


#setwd("/science/paolil/TARA_metaSNP/custom.analysis/")
setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/8_genetic_distances/")

#files = list.files("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered/pop/")
#species = gsub(".filtered.freq","",files)

species = "1002672"
#dist = ".core"
#dist.dir = "_universal"

#for (i in species){
i = species

#table.species = read.delim(paste0("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered",dist.dir,"/pop/",i,dist,".filtered.freq"),row.names=1)
table.species = read.delim("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/DATA/4_Run_genes_ann/universal_common/filtered/1002672.universal.common.filtered.freq",row.names=1)
names(table.species)=gsub("TARA_|_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam|_0.22.3_G.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam)","",names(table.species))

table.species[table.species==-1]=NA

table.species = cbind(read.table(text=row.names(table.species),sep=":",col.names=c("taxo.ann","gene.ann","position","SNP","significance")),table.species)
table.species = data.frame(table.species[,1:5],codon=table.species$significance,table.species[,6:length(table.species)])
table.species$significance=gsub("\\[.*\\]","",table.species$significance)
table.species$codon=gsub("N\\[|S\\[|\\]","",table.species$codon)
names(table.species)[7:length(table.species)]=gsub("X","",names(table.species)[7:length(table.species)])

head(table.species,10)
ncol(table.species)

gen.dist.species=as.data.frame(matrix(NA,ncol=length(table.species[7:length(table.species)]),nrow=length(table.species[7:length(table.species)])),col.names = names(table.species)[7:length(table.species)], row.names = names(table.species)[7:length(table.species)])
names(gen.dist.species)=row.names(gen.dist.species)
head(gen.dist.species)


system.time({
for (sample1 in row.names(gen.dist.species)){
  n = which(row.names(gen.dist.species)==sample1)
  for (sample2 in row.names(gen.dist.species)[1:n]){
    sum = 0
    count = 0
    for (k in unique(table.species$position)){
      vect.pos.1 = table.species[table.species$position==k,sample1]
      vect.pos.2 = table.species[table.species$position==k,sample2]
      
      gen.pos.1 = 1 - sum(vect.pos.1)
      gen.pos.2 = 1 - sum(vect.pos.2)
      
      temp = cbind(sample1=c(gen.pos.1,vect.pos.1),
                   sample2=c(gen.pos.2,vect.pos.2))
      sum.k = 0
      for (row in 1:nrow(temp)){
        sum.k=sum.k+sum(temp[row,"sample1"]*temp[-row,"sample2"])
      }
      
      if (!is.na(sum.k)){
        sum = sum + sum.k
        count = count +1
      }	
    }
    gen.dist.species[sample1,sample2]=sum/count
  }
}

print(gen.dist.species)
})
write.table(gen.dist.species,paste0(i,dist,".gen.dist.matrix.tab"))
#} 


