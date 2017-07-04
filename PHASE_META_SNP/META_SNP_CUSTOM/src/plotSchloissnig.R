#################################################################
# Compare the nucleotide diversity distribution between species #
#################################################################

# Requires the Pi tables (nucleotide diversity within and between samples)

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)


############
# Arguments
args = commandArgs(trailingOnly=TRUE)
wdirectory = args[1]
setwd(wdirectory)
###########

###########
# Libraries
library(tidyverse)
library(data.table)
library(gridExtra)
###########

files = list.files("pi/")
species.all = gsub("..genetic.distance.tab","",files)

### Classical Pi

plot.dataset = data.frame(id=character(),variable=character(),value=numeric(),species=character())

for (species in species.all){
  nucl.dist = read.delim(paste0("pi/",species,"..genetic.distance.tab"),row.names=1,check.names = F)
  # Taking only the intra sample diversity :
  nucl.dist[lower.tri(nucl.dist)]<-NA
  nucl.dist.m=melt(cbind(nucl.dist,id=row.names(nucl.dist)),id='id',na.rm=T)
  plot.dataset = rbind(plot.dataset,cbind(nucl.dist.m,species=species))
}

### Pi N / Pi S ratio

plot.ratio = data.frame(id=character(),variable=character(),value=numeric(),species=character())

for (species in species.all){
  nucl.ratio.N = read.delim(paste0("pi.N/",species,"N.genetic.distance.tab"),row.names=1,check.names = F)
  nucl.ratio.S = read.delim(paste0("pi.S/",species,"S.genetic.distance.tab"),row.names=1,check.names = F)
  nucl.ratio.N[lower.tri(nucl.ratio.N)]<-NA
  nucl.ratio.S[lower.tri(nucl.ratio.S)]<-NA
  nucl.ratio = nucl.ratio.N / nucl.ratio.S
  nucl.ratio.m=melt(cbind(nucl.ratio,id=row.names(nucl.ratio)),id='id',na.rm=T)
  plot.ratio = rbind(plot.ratio,cbind(nucl.ratio.m,species=species))
}

### Coverage information

# Load coverage information
coverage.table = read.delim(paste0("all.cov.stat.tab"),check.names = F)

# Order by increasing coverage
coverage.table = as.data.table(coverage.table)
coverage.table[, cumulative := sum(avg.cov),by = species]
coverage.table[, mean := mean(avg.cov),by = species]
setorder(coverage.table,cumulative)
coverage.table=as.data.frame(coverage.table)

# Store order for plot purposes
order = as.factor(unique(coverage.table$species))
coverage.table$species=factor(coverage.table$species, levels = order)

## Number of SNP per kb
nSNP.table = read.delim(paste0("table.SNP.stats.tab"),check.names = F)

# Reordering tables
nSNP.table$species=factor(nSNP.table$species,levels = order)
plot.dataset$species=factor(plot.dataset$species,levels = order)
plot.ratio$species=factor(plot.ratio$species,levels = order)

########
# PLOT #
########

pdf("full.analysis.pdf",height=29.7/2,width=21/2)
grid.arrange(

  ## NUMBER SNP / KB
  ggplot(nSNP.table)+geom_point(aes(x=species,y=SNP.per.Kb))+
    theme_grey()+ylab("SNPs per Kb")+
    scale_y_log10(limits=c(1,NA))+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_smooth(aes(x=as.integer(species),y=SNP.per.Kb),size=0.5),
  
  ## PI
  ggplot(plot.dataset)+geom_boxplot(aes(x=species,y=value))+
    theme_grey()+ylab("Nucleotide diversity (pi)")+
    scale_y_log10()+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_smooth(aes(x=as.integer(species),y=value),size=0.5)+
    geom_hline(yintercept = 0.01,colour="#990000", linetype="dashed",size=0.5),
  
  ## PIN/PIS
  ggplot(plot.ratio)+geom_boxplot(aes(x=species,y=value))+
    theme_grey()+ylab("piN/piS")+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_smooth(aes(x=as.integer(species),y=value),size=0.5),
  
  ## nN / nS
  ggplot(nSNP.table)+geom_point(aes(x=species,y=N.SNP/S.SNP))+
    theme_grey()+ylab("nN/nS")+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank())+
    geom_smooth(aes(x=as.integer(species),y=N.SNP/S.SNP),size=0.5),
  
  ## Multiallelic 
  ggplot(nSNP.table)+
    geom_point(aes(x=species,y=total.multiallelic,shape="All"))+
    geom_point(aes(x=species,y=N.multiallelic,shape="N"))+
    geom_point(aes(x=species,y=S.multiallelic,shape="S"))+
    theme_grey()+ylab("Multiallelic SNP proportion")+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position = c(0, 1), 
          legend.justification = c(0, 1)),
  
  ## H COVERAGE
  ggplot(coverage.table)+
    geom_boxplot(aes(x=species,y=hcov.1X))+ylab("Horizontal coverage (%)")+
    geom_smooth(aes(x=as.integer(species),y=hcov.1X),size=0.5)+
    theme_grey()+
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank()),
  ## COVERAGE
  ggplot(coverage.table)+
    geom_boxplot(aes(x=species,y=avg.cov))+
    scale_y_log10()+ylab("Coverage\n(cumulative\nand distribution)")+
    geom_point(aes(x=species,y=cumulative),shape=15)+
    theme_grey()+
    theme(axis.text.x=element_text(angle = 90)) + xlab("Species"),
 
  ncol=1)

dev.off()


  
