#################################################################
# Compare the nucleotide diversity distribution between species #
#################################################################

# Requires the Pi tables (nucleotide diversity within and between samples)

############
# Working Directory :
args = commandArgs(trailingOnly=TRUE)
dir = args[1]
#dir = "/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/TARA_POPULATION_GENETICS/RESULTS/6_Resultats_metaSNV/tara_joined.G-m20-d10-b40.T-m20-d10-b17/Plot_allele_dist/"
setwd(dir)
# Data directory
outdir = paste0(gsub("/[^/]*$", "", gsub('/$','',dir)),'/Plot_allele_dist/')
dir.create(file.path(outdir))
###########

###########
# Libraries
library(tidyverse)
library(data.table)
library(gridExtra)
library(cowplot)
###########

# Fast replacement for data.table :
fast.replace = function(DT,x=-1,y=NA) {
  # by number (slightly faster than by name) :
  for (j in seq_len(ncol(DT)))
    set(DT,which(DT[[j]]==x),j,y)
}
diff = function(x) {
  return(x[1] - x[2])
}

dom = function(x) {
 if (x[1]>=0.7){
  return('dominant')}
 else if (x[1]<=0.3){
  return('rare')}
 else{
  return('inter')}
}

#############
# Load Data #
#############
dir = paste0(dir,'/pop/')
files = list.files(dir)
species = gsub('.filtered.freq','',files[grep('.filtered.freq',files)])

print('Species found :\n')
print(species)

for (sp in species){
  print(sp)
  # Loading the data
  freq.file = fread(paste0(dir,sp,'.filtered.freq'))
  # Replacing '-1' by 'NA'
  fast.replace(freq.file)
  
  freq.file.m = melt(freq.file,id='V1')
  
  freq.file.m[, dataset := gsub('.filtered.*|TARA_.*_0.22-.*_','',freq.file.m$variable)]
  freq.file.m[, station := gsub('.filtered.*|TARA_|_0.22-.*','',freq.file.m$variable)]
  freq.file.m[, grp := paste(freq.file.m$station,freq.file.m$V1,sep=':')]
  
  freq.file.m = freq.file.m[, col := diff(value), by=grp] 
  freq.file.m = freq.file.m[!is.na(col)]
  freq.file.m = freq.file.m[col!=0]
  freq.file.m = freq.file.m[, dom := dom(value), by=grp]

  p1 = ggplot(freq.file.m) + geom_density(aes(x = value, y = ..scaled.., fill = dataset), alpha = 0.3) + theme_grey() +
    xlab('Allele frequency') + ggtitle('Density of the allele frequency') + ylab('scaled density') +
    theme(legend.position = c(.95, .95), legend.justification = c("right", "top"))
  
  p2 = ggplot(freq.file.m) + geom_line(aes(x = dataset, y = value, group = grp, col = abs(col)), alpha= 0.3) +
    theme_grey() + ylab('Allele Frequency') + ggtitle('Variability of allele frequency from G to T') +
    theme(legend.position = 'right')
  
  p3 = ggplot(freq.file.m) + geom_density(aes(x = col, fill=dom), alpha = 0.3) + theme_grey() +
    xlab('Variability of allele frequency from G to T') + ggtitle('Density of the allele variability') +
    theme(legend.position = c(.95, .95), legend.justification = c("right", "top"))
  
  plot = plot_grid(p1,p2,p3, ncol=3, rel_widths = c(0.75,1,0.75))
  
  print('saving plot...')
  ggsave(paste0(outdir,sp,'.allele.freq.pdf'), plot = plot, width = 12, height = 6)
}

