# =======================================================
# metaSNV part of the mOTUs paper
# Lucas Paoli, lucas.paoli@ens.fr ; lucas.paoli@gmail.com
# Privacy based on a multispecies approach
# =======================================================

# ===== Libraries =====
library(tidyverse)

# ===== Variables =====
rm(list=ls())
setwd("/Users/Lucas/Documents/Publications/2018_mOTUs/")
#source("Analysis/Load_data_for_privacy")
summary_hmp_m = read_tsv('Data/summary_hmp_mOTUs.tsv')
summary_hmp_m_sub = summary_hmp_m[summary_hmp_m$Part1 == summary_hmp_m$Part2,]
privacy_m = read_tsv('Data/privacy_mOTUs.tsv')
hmp_snps = read_delim("Data/Distances-final/hmp.final.motu.n.SNPs.txt", col_names = F, delim = ' ')
hmp_snps$X2 = gsub('.filt.*','',hmp_snps$X2) ; hmp_snps$X1 = hmp_snps$X1 - 1
names(hmp_snps) = c('N.snps', 'Species')
super.sites = read_tsv("Data/mOTUv2.insert.scaled.floor.rm0.site.map.tsv")
id = read_tsv('Data/Distances-final//map_genomes_2_ref-motus.tsv', col_names = F) ; id$X2 = gsub('\\..*', '', id$X2)
subspecies = read_xlsx('Biblio/inline-supplementary-material-4.xlsx', sheet = 2)

# ===== Functions =====
get_cosmopolitan_score <- function(V1,V2,V3,V4,V5, df = select.species.sub, n = -1, rate = .6, f.name = 'samples.names.txt'){
  x = c(V1,V2,V3,V4,V5)
  r = round((1-rate)*length(x))
  df.sub = df[, names(df) %in% x]
  df.bin = df.sub %>%
    rowwise() %>% 
    is.na()
  if (sum(rowSums(df.bin) <= r) == n){
    spl.names = row.names(df.bin)[which(rowSums(df.bin) <= r)]
    write(spl.names, file = f.name)
  }
  return(sum(rowSums(df.bin) <= r))
}
get_pangenome_privacy <- function(df, body.part=p){
  
  priv = data.frame(Indiv = unique(c(df$Indiv1, df$Indiv2)),
                    Super.Site = body.part,
                    Value = NA)
  df = df %>%
    rowwise() %>%
    mutate(order.comp = paste(sort(c(id,variable)),collapse=':'))
  
  df = df %>%
    group_by(order.comp) %>%
    mutate(n.obs.comp = length(value),
           cat.value = sum(value*N.snps)/sum(N.snps))
  
  df.cat = df[!duplicated(df$order.comp),]
  for (i in unique(c(df.cat$Indiv1, df.cat$Indiv2))){
    df.i = subset(df.cat, Indiv1 == i | Indiv2 == i)
    df.intra = subset(df.i, Indiv1 == i & Indiv2 == i)
    df.inter = subset(df.i, Indiv1 != i | Indiv2 != i)
    priv[which(priv$Indiv == i),'Value'] = ifelse(nrow(df.intra) > 0,
                                                  max(df.intra$cat.value) < min(df.inter$cat.value),
                                                  NA)
  }
  return(priv)
}


# ===== Cosmopolitan privacy =====

# Select best species for a 'cosmopolitan pangenome'
parts_freq = table(summary_hmp_m_sub[!duplicated(summary_hmp_m_sub$Species),]$SuperSite1)
parts = names(parts_freq[parts_freq>=5])

cosmopolitan_privacy = NULL
for (p in parts){
  print(p)
  df.p = subset(summary_hmp_m_sub, SuperSite1 == p)
  
  # Take samples with at least 5 species, then take species present in at least x% of samples
  # Create df with samples as rows, species as cols
  select.species = data.frame(Samples = unique(c(df.p$id, df.p$variable)))
  for (s in unique(df.p$Species)){
    df.s = subset(df.p, Species == s)
    species.i = data.frame(Samples = unique(c(df.s$id, df.s$variable)),
                           n.samples = unique(c(df.s$id, df.s$variable)))
    names(species.i)[2] = s
    row.names(species.i) = unique(c(df.s$id, df.s$variable))
    select.species = left_join(select.species, species.i)
  }
  
  # Get the best combinaison of n species to select as many samples
  # present in at least n-1 species
  n = 5
  select.species.sub = select.species[, 2:ncol(select.species)]
  # To speed up computation, remove species with more than 20% of NAs
  select.species.sub.speed = select.species.sub[, colSums(!is.na(select.species.sub)) > round(nrow(select.species.sub)/5)]
  row.names(select.species.sub.speed) = select.species$Samples
  df.combn = as.tibble(t(combn(names(select.species.sub.speed), n)))
  res.combn = df.combn %>% 
    rowwise() %>% 
    mutate(Score = get_cosmopolitan_score(V1,V2,V3,V4,V5, rate = .8))
  if (max(res.combn$Score) >= round(nrow(select.species.sub)/5)) {
    print('Speed up hypothesis is true -- All good!')
    print(paste0('Cosmopolitan analysis can include ',
                round(max(res.combn$Score)/nrow(select.species.sub)*100),
                '% of the ', nrow(select.species.sub), ' samples'))
  } else {
    print('Speed up hypothesis not true -- may be worth to double check')
  }
  
  ##
  ## Work In Progress
  ##
  
  # for (z in length(which(res.combn$Score == max(res.combn$Score)))) ... ?
  
  which(res.combn$Score == max(res.combn$Score))
  
  # rerun to id samples
  res.combn = df.combn %>% 
    rowwise() %>% 
    mutate(Score = get_cosmopolitan_score(V1,V2,V3,V4,V5, n = max(res.combn$Score), rate = .8))
  
  samples.list = read_csv('samples.names.txt', col_names = F)
  species.list = as.character(res.combn[which(res.combn$Score == max(res.combn$Score))[1],1:5])
  
  #Subset for the n 'cosmopolitan' species and corresponsing samples
  df.s = subset(df.p, Species %in% species.list)
  df.s = subset(df.s, id %in% samples.list$X1 & variable %in% samples.list$X1)
  df.s = left_join(df.s, snps, by = 'Species')
  cosmopolitan_privacy = rbind(cosmopolitan_privacy,get.privacy(df.s))
}

# ===== Pangenome privacy =====

Privacy.all = NULL
for (p in c('OR', 'ST', 'SK', 'VG')){
  df.p = subset(summary_hmp.m_sub, SuperSite1 == p)
  df.p = left_join(df.p, snps, by = 'Species')
  Privacy.all = rbind(Privacy.all,get.privacy(df.p))
}

Privacy.all = Privacy.all[!is.na(Privacy.all$Value),]
Privacy.all$Value = ifelse(Privacy.all$Value, 'Private', 'Not Private')
Privacy.all$Super.Site = factor(Privacy.all$Super.Site,labels=c('Oral','Skin','Stool','Vagina'))

# ===== Figures =====

p1 = ggplot(Privacy.all)+geom_histogram(aes(x=Value,fill=Value,col=Value), stat = 'count') +
  facet_wrap(~Super.Site) + ggtitle('Privacy of individual microbiome based on all available data') +
  theme_bw() + theme(axis.title.x = element_blank(), legend.position = 'none') +
  ylab('Number of Comparisons\n(two time points for an individual)') +
  scale_fill_manual(values=c('#f6c47d','#5085b6')) +
  scale_color_manual(values=c('#f6c47d','#5085b6')) 

ggsave('private.whole.microbiome.motu.m20.pdf', p1, height = 6, width = 5)

Privacy = Privacy[!is.na(Privacy$Value),]
Privacy$Value = ifelse(Privacy$Value, 'Private', 'Not Private')
Privacy$Super.Site = factor(Privacy$Super.Site,labels=c('Oral','Stool'))

p2 = ggplot(Privacy)+geom_histogram(aes(x=Value,fill=Value,col=Value), stat = 'count') +
  facet_wrap(~Super.Site) + 
  theme_bw() + theme(axis.title.x = element_blank(), legend.position = 'none') +
  ylab('Number of Comparisons\n(two time points for an individual)') +
  scale_fill_manual(values=c('#f6c47d','#5085b6')) +
  scale_color_manual(values=c('#f6c47d','#5085b6')) +
  ggtitle('Privacy of individual microbiome based on cosmopolitan species',
          subtitle = '5 species combinaison with the most samples detected in at least 4 of them')
ggsave('private.sub.microbiome.motu.m20.pdf', p2, height = 4, width = 6)

