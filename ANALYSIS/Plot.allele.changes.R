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
  diff = x[1] - x[2]
  if (!is.na(diff)){
    if (diff < 0){
      if (x[2] == 1){
        return('Fixed')
      }else{return('Higher')}
    }else if (diff > 0){
      if (x[2] == 0){
        return('Inactive')
      }else{return('Lower')}
    }else {
      if (x[2] == 0){
        return('No changes (Inact.)')
      }else if(x[2] == 1){
        return('No changes (Fixed)')
      }else{return('No changes')}}
}
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
  
  freq.file.m = freq.file.m[, Change := diff(value), by=grp] 
  freq.file.m = freq.file.m[!is.na(Change)]
  freq.file.m = freq.file.m[, dom := dom(value), by=grp]
  freq.file.m$Change = factor(freq.file.m$Change, levels = 
                           c('No changes (Fixed)','Fixed','Higher','No changes','Lower','Inactive','No changes (Inact.)'))

  plot = ggplot(freq.file.m) + geom_bar(aes(x = Change, fill=Change)) + theme_grey() +
    xlab('Variability of allele frequency from G to T') + ggtitle('Allele variability') +
    theme(legend.position = c(.95, .98), legend.justification = c("right", "top"),
    axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank()) + facet_wrap(~dom, nrow=3)
  
  
  print('saving plot...')
  ggsave(paste0(outdir,sp,'.allele.change.pdf'), plot = plot, width = 4, height = 10)
}

