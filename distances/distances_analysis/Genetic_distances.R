# =======================================================
# metaSNV part of the mOTUs paper
# Lucas Paoli, lucas.paoli@ens.fr ; lucas.paoli@gmail.com
# =======================================================

# ===== Libraries =====
library(tidyverse)
library(reshape2)
library(readxl)

# ===== Variables =====
rm(list=ls())
setwd("/Users/Lucas/Documents/Publications/2018_mOTUs/")
id = read_tsv('distances/map_genomes_2_ref-motus.tsv', col_names = F)
id$X2 = gsub('\\..*', '', id$X2)
dir_tara_motu = 'distances/tara.motu.distances-m5-d10-b80-p0.9'
dir_hmp_motu = 'distances/hmp.motu.distances-m20-d10-b80-p0.9'
dir_tara_f11 = 'distances/tara.freeze11.distances-m5-d10-b40-p0.9'
dir_hmp_f11 = 'distances/hmp.freeze11.distances-m20-d10-b40-p0.9'
tara_snps = "distances/tara.final.motu.n.SNPs-m5.txt"
hmp_snps = "distances/hmp.final.motu.n.SNPs.txt"
subspecies = read_xlsx('inline-supplementary-material-4.xlsx', sheet = 2)

# ===== Functions =====
get.sp.distances <- function(sp, folder.freeze, folder.motus, num.snps, plot.folder, id_map, patf11 = '_G.*', patmot = '_G.*', dataset = 'Tara'){
  
  # ----- Load & format data -----
  snps = read_delim(num.snps, col_names = F, delim = ' ')
  snps$X2 = gsub('.filt.*','',snps$X2)
  snps$X1 = snps$X1 - 1
  
  freeze11 = read.delim(file = paste0(folder.freeze, '/', id_map[which(id_map$X1 == sp),]$X2, '.filtered.mann.dist'), sep = '\t', check.names = F, row.names = 1, stringsAsFactors = F)
  motus = read.delim(file = paste0(folder.motus, '/', sp, '.filtered.mann.dist'), sep = '\t', check.names = F, row.names = 1, stringsAsFactors = F)
  
  # if (dataset == 'HMP'){
  #   freeze11 = freeze11[grepl('[0-9]{9}-[a-z]', names(freeze11)), grepl('[0-9]{9}-[a-z]', names(freeze11))]
  #   motus = motus[grepl('[0-9]{9}-[a-z]', names(motus)), grepl('[0-9]{9}-[a-z]', names(motus))]
  #   if (length(nrow(freeze11)) == 0 | length(nrow(motus)) == 0){
  #     return(NA)
  #   }
  #   if (nrow(freeze11) < 20 | nrow(motus) < 20){
  #     return(NA)
  #   }
  # }
  
  names(freeze11) = gsub(patf11, '_G', names(freeze11)) ; row.names(freeze11) = gsub(patf11, '_G', row.names(freeze11))
  names(motus) = gsub(patmot, '_G', names(motus)) ; row.names(motus) = gsub(patmot, '_G', row.names(motus))
  names.subset = intersect(names(freeze11), names(motus))
  
  freeze11 = freeze11[names.subset,names.subset] ; motus = motus[names.subset,names.subset]
  freeze11[upper.tri(freeze11, diag = T)] = NA ; motus[upper.tri(motus, diag = T)] = NA
  
  freeze11.m = melt(cbind(freeze11, id = row.names(freeze11)), id = 'id', na.rm = T)
  motus.m = melt(cbind(motus, id = row.names(motus)), id = 'id', na.rm = T)
  
  # ----- Correlation & Plot -----
  corr.value = data.frame(Pearson = cor.test(freeze11.m$value, motus.m$value)$estimate,
                          Spearman = cor.test(freeze11.m$value, motus.m$value, method = 'spearman')$estimate, 
                          Species = id_map[which(id_map$X1 == sp),]$X2,
                          mOTUs = sp,
                          snps = snps[which(snps$X2 == sp),]$X1)
  
  p = ggplot() +
    geom_point(aes(x=(freeze11.m$value), y=(motus.m$value)), shape=1) +
    geom_abline(slope=1) + xlab('Freeze11 distances') + ylab('mOTUs distances') +
    theme_grey() + theme(legend.position='none') +
    ggtitle(paste0('qqPlot, species : ', sp, ' / mOTU : ', id[which(id$X2 == sp),]$X1),
            subtitle = paste0('Correlation : coeff = ',
                              round(cor.test(freeze11.m$value, motus.m$value)$estimate, digits = 4),
                              ', p-value = ',
                              round(cor.test(freeze11.m$value, motus.m$value)$p.value, digits = 3),
                              ' ; Number of SNVs = ', snps[which(snps$X2 == id[which(id$X2 == sp),]$X1),]$X1))
  
  ggsave(paste0(plot.folder, '/', sp, '.correlation.plot.pdf'),p)
  
  return(corr.value)
}
get.body.part <- function(sp, folder.freeze, body.parts){
  freeze11 = read.delim(file = paste0(folder.freeze, '/', sp, '.filtered.mann.dist'), sep = '\t', check.names = F, row.names = 1, stringsAsFactors = F)
  subseting = body.parts$SampleID %in% gsub('.freeze.*','',names(freeze11))
  body.parts = body.parts[subseting,]
  res = ifelse(nrow(body.parts) == 0, 'ST', names(table(body.parts$Site)[1]))
  return(res)
}

# ===== Tara =====
summary_tara = NULL
tara_species_motus = gsub('.filt.*', '', list.files(dir_tara_motu)[grepl('ref_.*mann',list.files(dir_tara_motu))])
tara_species_freeze = gsub('.filt.*', '', list.files(dir_tara_f11)[grepl('mann',list.files(dir_tara_f11))])
id_tara = id[id$X2 %in% tara_species_freeze & id$X1 %in% tara_species_motus,]
id_tara = id_tara[!duplicated(id_tara$X1),]

for (spe in id_tara$X1){
  corr.value = get.sp.distances(spe,
                                dir_tara_f11,
                                dir_tara_motu,
                                tara_snps,
                                'Distance_correlations/hmp.correlations',
                                id_tara)
  summary_tara = rbind(summary_tara, corr.value)
}

p_tara = ggplot(summary_tara, aes(x=Pearson, y=Spearman)) + theme_grey() +
  stat_ellipse(aes(x=Pearson, y=Spearman), inherit.aes = F) +
  geom_point(aes(size = snps, color = Species), alpha = .7) +
  ylim(0,NA)+xlim(0,NA)+
  #geom_label(aes(label = Species))+
  theme(legend.position = 'none')+
  ggtitle('Correlation for mOTUs and freeze11 genome distances', subtitle = 'Allele prevalence of 90%')

ggsave('Distance_correlations/Tara.correlations.summary.pdf',p_tara)


# ===== HMP =====
summary_hmp = NULL
hmp_species_motus = gsub('.filt.*', '', list.files(dir_hmp_motu)[grepl('ref_.*mann',list.files(dir_hmp_motu))])
hmp_species_freeze = gsub('.filt.*', '', list.files(dir_hmp_f11)[grepl('mann',list.files(dir_hmp_f11))])
id_hmp = id[id$X2 %in% hmp_species_freeze & id$X1 %in% hmp_species_motus,]
id_hmp = id_hmp[!duplicated(id_hmp$X1),]

for (spe in id_hmp$X1){
  corr.value = get.sp.distances(spe,
                                dir_hmp_f11,
                                dir_hmp_motu,
                                hmp_snps,
                                'Distance_correlations/hmp.correlations',
                                id_hmp,
                                patf11 = '.freeze11.*',
                                patmot = '.snv.bam',
                                dataset = 'HMP')
  summary_hmp = rbind(summary_hmp, corr.value)
}

p_hmp = ggplot(summary_hmp, aes(x=Pearson, y=Spearman)) + theme_grey() +
  stat_ellipse(aes(x=Pearson, y=Spearman), inherit.aes = F) +
  geom_point(aes(size = snps, color = Species), alpha = .7) +
  ylim(0,NA)+xlim(0,NA)+
  #geom_label(aes(label = Species))+
  theme(legend.position = 'none')+
  ggtitle('Correlation for mOTUs and freeze11 genome distances', subtitle = 'Allele prevalence of 90%')

ggsave('Distance_correlations/hmp.correlations.summary.pdf',p_hmp)

# ===== Summary =====
summary_table = as.tibble(rbind(data.frame(summary_hmp, Dataset = 'HMP'),
                                data.frame(summary_tara, Dataset = 'Tara')))
summary_table = summary_table[!is.na(summary_table$Pearson),]
super.sites = read_tsv("Body_parts_privacy/mOTUv2.insert.scaled.floor.rm0.site.map.tsv")
names(super.sites) = c('SampleID', 'Site')

summary_table$Site = NA
summary_table$Species = as.character(summary_table$Species)

for (s in subset(summary_table, Dataset == 'HMP')$Species){
  site = get.body.part(s, dir_hmp_f11, super.sites)
  summary_table[which(summary_table$Species == s), 'Site'] = site
}
summary_table[is.na(summary_table$Site), 'Site'] = 'Tara'
summary_table$Site = factor(summary_table$Site,labels=c(
  paste0('Oral\n(', nrow(subset(summary_table, Site %in% c('OR','Oral'))), ' Species)'),
  paste0('Skin\n(', nrow(subset(summary_table, Site %in% c('SK','Skin'))), ' Species)'),
  paste0('Stool\n(', nrow(subset(summary_table, Site %in% c('ST','Stool'))), ' Species)'),
  paste0('Ocean\n(', nrow(subset(summary_table, Site %in% c('Tara','Ocean'))), ' Species)')))

p_summary.site = ggplot(summary_table, aes(x=Pearson, y=Spearman)) + theme_grey() +
  stat_ellipse(aes(x=Pearson, y=Spearman, color = Site), inherit.aes = F) +
  geom_point(aes(size = snps, color = Site), alpha = .7) +
  ylim(0,NA)+xlim(0,NA)+
  theme_bw()+
  #geom_label(aes(label = Species))+
  ggtitle('Correlation of mOTUs and freeze11 genetic distances',
          subtitle = 'HMP : min 20 Samples, Tara : min 5 Samples') +
  scale_color_manual(values=c("#98be57","#91b5b5","#c7624f","#9555b4","#514143")) +
  scale_size_continuous(name = 'Num. of SNVs')


ggsave('correlations.summary.sitewise.pdf',p_summary.site)

p_boxplot.site = ggplot(summary_table, aes(x = Site, y=Pearson, fill=Site)) + theme_grey() +
  geom_boxplot() +
  theme_bw()+
  ylab('Freeze11/mOTUs Pearson coeff.') +
  #geom_label(aes(label = Species))+
  ggtitle('Correlation of mOTUs and freeze11 genetic distances') +
  #scale_fill_manual(values=c("#98be57","#91b5b5","#c7624f","#9555b4","#514143"))
  scale_fill_manual(values=c("#9ac161","#738386","#b65a49","#8e50ab")) +
  theme(legend.position = 'none', axis.title.x = element_blank())

ggsave('correlations.summary.boxplot.pdf',p_boxplot.site)


p_pars = ggplot(summary_table) +
  geom_smooth(aes(x=snps,y=Pearson), color = 'black') +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = mean(summary_table$Pearson), linetype = 2) +
  geom_point(aes(x=snps,y=Pearson, color = Site),size=2.5,alpha=.9) +
  ylim(0, NA)+
  scale_x_log10() +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  xlab('Number of SNVs') +
  ylab('Freeze11/mOTUs Pearson coeff.') +
  ggtitle('Correlation coeff. depending on log(Number of SNVs)') +
  scale_color_manual(values=c("#9ac161","#738386","#b65a49","#8e50ab"))

ggsave('correlations.pars.summary.pdf',p_pars)


summary_bad_corr = subset(summary_table, Pearson < .7)
summary_good_corr = subset(summary_table, Pearson >= .7)
duplicated(summary_bad_corr$Species)
summary_bad_corr$Species %in% summary_good_corr$Species
summary_bad_corr[summary_bad_corr$Species %in% summary_good_corr$Species,]

subsp_plot = subspecies[,c('Ref.genome (NCBI)', 'Nr. Subspecies')]
subsp_plot$Subspecies = subsp_plot$`Nr. Subspecies` >= 1
names(subsp_plot) = c('Species', 'Nr.Sub', 'Subspecies')
summary_table_subsp = left_join(summary_table, subsp_plot)
summary_table_subsp[!(summary_table_subsp$Subspecies)&!is.na(summary_table_subsp$Subspecies),]$Subspecies = NA
summary_table_subsp$bad_corr = summary_table_subsp$Pearson < .4
summary_table_subsp[summary_table_subsp$bad_corr,]

p_pars_sub = ggplot(summary_table_subsp) +
  geom_smooth(aes(x=snps,y=Pearson), color = 'black') +
  geom_hline(yintercept = 1) +
  geom_hline(yintercept = mean(summary_table$Pearson), linetype = 2) +
  geom_point(aes(x=snps,y=Pearson, color = Site),size=2.5,alpha=.9) +
  geom_point(aes(x=snps,y=Pearson,shape=Subspecies), alpha=.6) +
  ylim(0, NA)+
  scale_x_log10() +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  xlab('Number of SNVs') +
  ylab('Freeze11/mOTUs Pearson coeff.') +
  ggtitle('Correlation coeff. depending on log(Number of SNVs)') +
  scale_color_manual(values=c("#9ac161","#738386","#b65a49","#8e50ab")) +
  scale_shape_manual(values=3)

ggsave('correlations.pars.subsp.pdf',p_pars_sub)
write_tsv(x = summary_table_subsp, 'correlation_values.tsv')
