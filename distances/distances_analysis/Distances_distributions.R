# =======================================================
# metaSNV part of the mOTUs paper
# Lucas Paoli, lucas.paoli@ens.fr ; lucas.paoli@gmail.com
# Plot distances distributions
# =======================================================

# ===== Libraries =====
library(tidyverse)
library(readxl)

# ===== Variables =====
rm(list=ls())
setwd("/Users/Lucas/Documents/Publications/2018_mOTUs/")
#source("Analysis/Load_data_for_privacy")
data_m = read_tsv('Data/summary_hmp_mOTUs.tsv')
data_m[which(data_m$id == "763901136-supragingivalplaque1"), 'SuperSite1'] = 'OR'
data_m[which(data_m$variable == "763901136-supragingivalplaque1"), 'SuperSite2'] = 'OR'
data_f = read_tsv('Data/summary_hmp_freeze11.tsv')
data_f[which(data_f$id == "763901136-supragingivalplaque1"), 'SuperSite1'] = 'OR'
data_f[which(data_f$variable == "763901136-supragingivalplaque1"), 'SuperSite2'] = 'OR'
hmp_snps = read_delim("Data/Distances-final/hmp.final.motu.n.SNPs.txt", col_names = F, delim = ' ')
hmp_snps$X2 = gsub('.filt.*','',hmp_snps$X2) ; hmp_snps$X1 = hmp_snps$X1 - 1
names(hmp_snps) = c('N.snps', 'Species')
super.sites = read_tsv("Data/mOTUv2.insert.scaled.floor.rm0.site.map.tsv")
id = read_tsv('Data/Distances-final//map_genomes_2_ref-motus.tsv', col_names = F) ; id$X2 = gsub('\\..*', '', id$X2)
subspecies = read_xlsx('Biblio/inline-supplementary-material-4.xlsx', sheet = 2)

# ===== Intersection between mOTUs and Freeze =====

data_compare = data_m[data_m$Species %in% id$X1, ]
names(data_compare)[4] = "Species.mOTUs"
data_compare_f = data_f[data_f$Species %in% id$X2, ]
data_compare_f$Species.mOTUs = NA

for (species in unique(data_compare_f$Species)){
  data_compare_f[which(data_compare_f$Species == species),]$Species.mOTUs = id[which(id$X2 == species)[1],]$X1
}

data_compare_f$Ref = paste0(data_compare_f$id,':',
                            data_compare_f$variable,':',
                            data_compare_f$Species.mOTUs)

data_compare$Ref = paste0(data_compare$id,':',
                          data_compare$variable,':',
                          data_compare$Species.mOTUs)

data_compare = left_join(data_compare,data_compare_f, by='Ref', suffix = c(".m", ".f"))
data_compare = data_compare[!is.na(data_compare$value.f),]

# ===== Figures =====

data_compare_plot = tibble(Distances=c(data_compare$value.m,data_compare$value.f),
                           Reference=rep(c('mOTUs','Freeze11'),each=nrow(data_compare)),
                           Individual=c(data_compare$Indiv.m,data_compare$Indiv.f),
                           Site=c(data_compare$SuperSite1.m,data_compare$SuperSite1.f))

p1a = ggplot(data_compare_plot) + 
  ggtitle('Distributions of distances', 
          subtitle = 'Matches between mOTUs and Freeze11 (includes all species)')+
  geom_density(aes(x=Distances,colour=Reference,fill=Reference), alpha=.5) +
  facet_grid(Individual~Site) + theme_bw() + theme(axis.title = element_blank())+
  theme(legend.position = 'bottom') +
  scale_fill_manual(values=c('#f6c47d','#5085b6')) +
  scale_color_manual(values=c('#f6c47d','#5085b6')) 

data_all_plot = tibble(Distances=c(data_m$value,data_f$value),
                       Reference=c(rep('mOTUs', nrow(data_m)),
                                   rep('Freeze11', nrow(data_f))),
                       Individual=c(data_m$Indiv,data_f$Indiv),
                       Species=c(data_m$Species, data_f$Species),
                       Site=c(data_m$SuperSite1,data_f$SuperSite1))

p1b = ggplot(data_all_plot) + 
  ggtitle('Distributions of distances',
          subtitle = 'All available data (includes all species)')+
  geom_density(aes(x=Distances,colour=Reference,fill=Reference), alpha=.5) +
  facet_grid(Individual~Site) + theme_bw() + theme(axis.title = element_blank())+
  scale_fill_manual(values=c('#f6c47d','#5085b6')) + theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'))

data_m$Type = ifelse(grepl('ref_mOTU', data_m$Species), 'Ref', 'Meta')
table(data_m[grepl('ref_mOTU', data_m$Species), 'SuperSite1'])
table(data_m[grepl('meta_mOTU', data_m$Species), 'SuperSite1'])

p1c = ggplot(data_m) + 
  ggtitle('Distributions of distances',
          subtitle = 'mOTUs : Ref vs Meta')+
  geom_density(aes(x=value,colour=Type,fill=Type), alpha=.5) +
  facet_grid(Indiv~SuperSite1) + theme_bw() + theme(axis.title = element_blank())+
  scale_fill_manual(values=c('#f6c47d','#5085b6')) + theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'))

ggsave('Results/distances.2.0-1A.pdf', p1a, height = 5, width = 5)
ggsave('Results/distances.2.0-1B.pdf', p1b, height = 5, width = 5)
ggsave('Results/distances.2.0-1C.pdf', p1c, height = 5, width = 5)
