# =======================================================
# metaSNV part of the mOTUs paper
# Lucas Paoli, lucas.paoli@ens.fr ; lucas.paoli@gmail.com
# Plot privacy distributions
# =======================================================

# ===== Libraries =====
library(tidyverse)
library(readxl)

# ===== Variables =====
rm(list=ls())
setwd("/Users/Lucas/Documents/Publications/2018_mOTUs/")
#source("Analysis/Load_data_for_privacy")
privacy_compare = read_tsv('Data/privacy_compare.tsv')
privacy_m = read_tsv('Data/privacy_mOTUs.tsv')
privacy_f = read_tsv('Data/privacy_freeze11.tsv')
hmp_snps = read_delim("Data/Distances-final/hmp.final.motu.n.SNPs.txt", col_names = F, delim = ' ')
hmp_snps$X2 = gsub('.filt.*','',hmp_snps$X2) ; hmp_snps$X1 = hmp_snps$X1 - 1
names(hmp_snps) = c('N.snps', 'Species')
super.sites = read_tsv("Data/mOTUv2.insert.scaled.floor.rm0.site.map.tsv")
id = read_tsv('Data/Distances-final//map_genomes_2_ref-motus.tsv', col_names = F) ; id$X2 = gsub('\\..*', '', id$X2)
subspecies = read_xlsx('Biblio/inline-supplementary-material-4.xlsx', sheet = 2)

# ===== Figures =====

privacy_compare_plot = tibble(Diff=c(privacy_compare$Value.m,privacy_compare$Value.f),
                              Reference=rep(c('mOTUs','Freeze11'),each=nrow(privacy_compare)),
                              Site=factor(c(privacy_compare$Super.Site.m,privacy_compare$Super.Site.f),labels=c('Oral','Stool','Skin','Vagina')))
p1a = ggplot(privacy_compare_plot) + 
  geom_vline(xintercept = 0, linetype=2) +
  ggtitle('Distributions of differences between self and closest other', 
          subtitle = 'Matches between mOTUs and Freeze11 (includes all species)')+
  geom_density(aes(x=Diff,colour=Reference,fill=Reference), alpha=.5) +
  facet_wrap(~Site) + theme_bw() + theme(axis.title = element_blank())+
  scale_fill_manual(values=c('#f6c47d','#5085b6')) + theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'))

privacy_all_plot = tibble(Diff=c(privacy_m$Value,privacy_f$Value),
                          Reference=c(rep('mOTUs', nrow(privacy_m)),
                                      rep('Freeze11', nrow(privacy_f))),
                          Species=c(privacy_m$Species.mOTUs, privacy_f$Species),
                          Site=factor(c(privacy_m$Super.Site,privacy_f$Super.Site),labels=c('Oral','Stool','Skin','Vagina')))
p1b = ggplot(privacy_all_plot) + 
  geom_vline(xintercept = 0, linetype=2) +
  ggtitle('Distributions of differences between self and closest other',
          subtitle = 'All available data (includes all species)')+
  geom_density(aes(x=Diff,colour=Reference,fill=Reference), alpha=.5) +
  facet_wrap(~Site) + theme_bw() + theme(axis.title = element_blank())+
  scale_fill_manual(values=c('#f6c47d','#5085b6')) + theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'))

privacy_m$Type = ifelse(grepl('ref_mOTU', privacy_m$Species.mOTUs), 'Ref', 'Meta')
table(privacy_m[grepl('ref_mOTU', privacy_m$Species.mOTUs), 'Super.Site'])
table(privacy_m[grepl('meta_mOTU', privacy_m$Species.mOTUs), 'Super.Site'])
p1c = ggplot(privacy_m) + 
  geom_vline(xintercept = 0, linetype=2) +
  ggtitle('Distributions of differences between self and closest other',
          subtitle = 'mOTUs : Ref vs Meta')+
  geom_density(aes(x=Value,colour=Type,fill=Type), alpha=.5) +
  facet_wrap(~Super.Site) + theme_bw() + theme(axis.title = element_blank())+
  scale_fill_manual(values=c('#f6c47d','#5085b6')) + theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'))
t.test(subset(privacy_m, Super.Site == 'Stool' & Type == 'Ref')$Value,
       subset(privacy_m, Super.Site == 'Stool' & Type == 'Meta')$Value)

p2a = ggplot(privacy_compare) + 
  xlab('Reference: mOTUs') + ylab('Reference: Freeze11') +
  ggtitle('Difference between self and closest other', subtitle = '(includes all species)')+
  geom_point(aes(x=Value.m,y=Value.f,colour = Sign), alpha=.8) + 
  facet_wrap(~Super.Site.m) + theme_bw() +  theme(legend.position = 'bottom') +
  scale_color_manual(values=c('#f6c47d','#5085b6'), name = 'Privacy')

p2b = ggplot(privacy_compare) + 
  ggtitle('Conservation of privacy between Freeze11 and mOTUs')+
  geom_histogram(aes(x=Sign,colour = Sign,fill=Sign), stat='count') + 
  facet_wrap(~Super.Site.m) + theme_bw() + theme(axis.title = element_blank(), legend.position = 'none') +
  scale_color_manual(values=c('#f6c47d','#5085b6'), name = 'Privacy') +
  scale_fill_manual(values=c('#f6c47d','#5085b6'), name = 'Privacy')

p3 = ggplot(privacy_m)+geom_histogram(aes(x=Privacy,fill=Privacy,col=Privacy), stat = 'count') +
  facet_wrap(~Super.Site) + ggtitle('Privacy of each available comparisons', subtitle='Reference: mOTUs') +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position = 'none') +
  ylab('Number of Comparisons') +
  scale_fill_manual(values=c('#f6c47d','#5085b6')) +
  scale_color_manual(values=c('#f6c47d','#5085b6')) 

ggsave('Results/privacy.2.0-1A.pdf', p1a, height = 5, width = 5)
ggsave('Results/privacy.2.0-1B.pdf', p1b, height = 5, width = 5)
ggsave('Results/privacy.2.0-1C.pdf', p1c, height = 5, width = 5)
ggsave('Results/privacy.2.0-2A.pdf', p2a, height = 5, width = 4)
ggsave('Results/privacy.2.0-2B.pdf', p2b, height = 5, width = 5)
ggsave('Results/privacy.2.0-3.pdf', p3, height = 5, width = 5)
