# =======================================================
# metaSNV part of the mOTUs paper
# Lucas Paoli, lucas.paoli@ens.fr ; lucas.paoli@gmail.com
# Load Data
# =======================================================

# ===== Libraries =====
library(tidyverse)
library(reshape2)
library(readxl)

# ===== Variables =====
rm(list=ls())
setwd("/Users/Lucas/Documents/Publications/2018_mOTUs/")
dir_hmp_motu = 'Data/Distances-final/hmp.motu.distances-m20-d10-b80-p0.9'
dir_hmp_f11 = 'Data/Distances-final/hmp.freeze11.distances-m20-d10-b40-p0.9'
super.sites = read_tsv("Data/mOTUv2.insert.scaled.floor.rm0.site.map.tsv")
id = read_tsv('Data/Distances-final//map_genomes_2_ref-motus.tsv', col_names = F) ; id$X2 = gsub('\\..*', '', id$X2)
subspecies = read_xlsx('Biblio/inline-supplementary-material-4.xlsx', sheet = 2)

# ===== Functions =====
load_distances <- function(file, folder, pat){
  # ----- Load data -----
  dist = read.delim(file = paste0(folder, '/', file), sep = '\t', check.names = F, row.names = 1, stringsAsFactors = F)
  names(dist) = gsub(pat, '', names(dist)) ; row.names(dist) = gsub(pat, '', row.names(dist))
  
  # ----- Select only HMP1 -----
  dist = dist[grepl('[0-9]{9}-[a-z]', names(dist)), grepl('[0-9]{9}-[a-z]', names(dist))]
  if (length(nrow(dist)) == 0){
    return(0)
  } else if (nrow(dist) < 20) {
    return(0)
  }
  
  # ----- Format data -----
  dist[upper.tri(dist, diag = T)] = NA
  
  dist.m = melt(cbind(dist, id = row.names(dist)), id = 'id', na.rm = T)
  dist.m$Species = gsub('.filt.*','',file)
  
  dist.m$Indiv1 = gsub('-.*', '', dist.m$id)
  dist.m$Indiv2 = gsub('-.*', '', dist.m$variable)
  
  dist.m$Part1 = gsub('[0-9]|-', '', dist.m$id)
  dist.m$Part2 = gsub('[0-9]|-', '', dist.m$variable)
  
  dist.m$Sample1 = gsub('.*-[a-z]*', '', dist.m$id)
  dist.m$Sample2 = gsub('.*-[a-z]*', '', dist.m$variable)
  
  return(dist.m)
}
create_summary <- function(dir, pattern, sites = super.sites){
  files = list.files(dir)[grep('mann',list.files(dir))]
  res = NULL
  for (f in files){
    dist = load_distances(file = f,
                          folder = dir,
                          pat = pattern)
    if (length(dist) != 1){
      res = rbind(res, dist)
    }
  }
  
  res$Indiv = ifelse(res$Indiv1 == res$Indiv2, 'Same', 'Different')
  names(sites) = c('id', 'SuperSite1')
  res = left_join(res, sites, by = 'id')
  names(sites) = c('variable', 'SuperSite2')
  res = left_join(res, sites, by = 'variable')
  
  return(res)
}
get_privacy <- function(df, body.part){
  
  df = df %>%
    rowwise() %>%
    mutate(order.comp = paste0(paste(sort(c(id,variable)),collapse=':'), ':', Species),
           IndivPart1 = paste0(Indiv1, ':', Part1, ':', Species),
           IndivPart2 = paste0(Indiv2, ':', Part2, ':', Species))
  
  priv = data.frame(Comp = unique(c(df$IndivPart1, df$IndivPart2)),
                    Super.Site = body.part,
                    Value = NA)
  
  for (i in priv$Comp){
    df.i = subset(df, IndivPart1 == i | IndivPart2 == i)
    df.intra = subset(df.i, Indiv == 'Same')
    df.inter = subset(df.i, Indiv == 'Different')
    priv[which(priv$Comp == i),'Value'] = ifelse(nrow(df.intra) > 0 & nrow(df.inter) > 0,
                                                 max(df.intra$value) - min(df.inter$value),
                                                 NA)
  }
  return(priv)
}
privacy_summary <- function(summary, motus = T){
  res = NULL
  for (p in c('OR', 'ST', 'SK', 'VG')){
    df.p = subset(summary, SuperSite1 == p)
    res = rbind(res,get_privacy(df.p, body.part = p))
  }
  
  res = res[!is.na(res$Value),]
  res$Super.Site = factor(res$Super.Site,labels=c('Oral','Stool','Skin','Vagina'))
  res$Privacy = ifelse(res$Value<0, 'Private','Not Private')
  
  if (motus){
    res$Species.mOTUs = gsub('.*:.*:','',res$Comp)
  } else if (!motus) {
    res$Species = gsub('.*:.*:','',res$Comp)
  }
  res$Indiv = gsub(':.*','',res$Comp)
  res$Part = gsub('[0-9]|:|meta_.*|ref_.*','',res$Comp)
  
  return(res)
  
}

# ===== mOTUs =====

# ----- Create Summary -----
# Get master table with distances, samples, species...
summary_hmp_m = create_summary(dir = dir_hmp_motu, pattern = '.snv.*')
summary_hmp_m_sub = summary_hmp_m[summary_hmp_m$Part1 == summary_hmp_m$Part2,]

# ----- Get Privacy -----
privacy_m = privacy_summary(summary_hmp_m_sub)

# ===== Freeze 11 =====

# ----- Create Summary -----
# Get master table with distances, samples, species...
summary_hmp_f = create_summary(dir = dir_hmp_f11, pattern = '.freeze11.*')
summary_hmp_f_sub = summary_hmp_f[summary_hmp_f$Part1 == summary_hmp_f$Part2,]

# ----- Get Privacy -----
privacy_f = privacy_summary(summary_hmp_f_sub, motus = F)

# ===== Intersection =====

privacy_compare = subset(privacy_m, Species.mOTUs %in% id$X1)
df_compare_f = subset(privacy_f, Species %in% id$X2)
df_compare_f$Species.mOTUs = NA

for (species in unique(df_compare_f$Species)){
  df_compare_f[which(df_compare_f$Species == species),]$Species.mOTUs = id[which(id$X2 == species)[1],]$X1
}
df_compare_f$Ref = paste0(df_compare_f$Indiv,':',df_compare_f$Part,':',df_compare_f$Species.mOTUs)
privacy_compare$Ref = privacy_compare$Comp

privacy_compare = left_join(privacy_compare,df_compare_f, by='Ref', suffix = c(".m", ".f"))
privacy_compare = privacy_compare[!is.na(privacy_compare$Value.f),]
privacy_compare$Sign = ifelse(sign(privacy_compare$Value.m) == sign(privacy_compare$Value.f),'Remained','Changed')

# ===== Write tables =====

write_tsv(summary_hmp_m, 'summary_hmp_mOTUs.tsv')
write_tsv(summary_hmp_f, 'summary_hmp_freeze11.tsv')
write_tsv(privacy_m, 'privacy_mOTUs.tsv')
write_tsv(privacy_f, 'privacy_freeze11.tsv')
write_tsv(privacy_compare, 'privacy_compare.tsv')

