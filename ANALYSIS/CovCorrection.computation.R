####################################
# Correct for full genome division #
####################################

# Requires the Pi tables (nucleotide diversity within and between samples)

# Methodology based on Schloissnig et al. 2013 (Nature)
# Written by lucas paoli (lucas.paoli@ens.fr ; lucas.paoli@gmail.com)


############
# Arguments
wdirectory = "~/GUT_POPULATION_GENETICS/DATA/GUT_metaSNV/GUT_2_metaSNV/"
plot.dir = "~/GUT_POPULATION_GENETICS/RESULTS"
#setwd(wdirectory)
setwd(plot.dir)
setwd("/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/GUT_POPULATION_GENETICS/RESULTS/")
###########

###########
# Libraries
library(tidyverse)
library(data.table)
library(gridExtra)
library(cowplot)
###########

files = list.files("genetic_final.1/pi/")
species = gsub(".genetic.distance.tab","",files)

# Load coverage information
coverage.table = read.delim(paste0("statistics_gut2_3_17_04/all.cov.stat.tab"),check.names = F)
coverage.table$sample=gsub('.*/','',coverage.table$sample)

for (i in species){
  ### Classical Pi
  nucl.dist = read.delim(paste0("genetic_final.1/pi/",i,".genetic.distance.tab"),row.names=1,check.names = F)
  # Subset for the species of interest
  coverage.table.s = coverage.table[which(coverage.table$species==i),]
  # Creating a table to store the min coverage for each paiwise comparison
  coverage.correction = nucl.dist
  coverage.correction[,] = NA
  
  for (sample in names(coverage.correction)){
    cov = coverage.table.s[which(coverage.table.s$sample==paste0(sample,'.cov.summary')),'avg.cov']
    cov = cov[!duplicated(cov)]
    cov = ifelse(cov >= 5, cov, NA)
    corr = cov/(cov-1)
    coverage.correction[sample,sample] = corr
  }
  coverage.correction[lower.tri(coverage.correction)]=1
  nucl.dist = nucl.dist * coverage.correction
  write.table(nucl.dist,paste0("genetic_final.1/pi.cov-1/",i,".genetic.distance.tab"),sep="\t")
}

# for (i in species){
#   ### Classical Pi
#   nucl.dist = read.delim(paste0("genetic_final.1/pi/",i,".genetic.distance.tab"),row.names=1,check.names = F)
#   # Subset for the species of interest
#   coverage.table.s = coverage.table[which(coverage.table$species==i),]
#   # Creating a table to store the min coverage for each paiwise comparison
#   coverage.correction = nucl.dist
#   coverage.correction[,] = NA
#   
#   for (sample1 in names(coverage.correction)){
#     pos.sample1 = which(names(coverage.correction)==sample1)
#     for (sample2 in names(coverage.correction)[1:pos.sample1]){
#       coverage.correction[sample1,sample2]=min(
#         coverage.table.s[which(coverage.table.s$sample==paste0(sample1,'.cov.summary')),'hcov.1X'],
#         coverage.table.s[which(coverage.table.s$sample==paste0(sample2,'.cov.summary')),'hcov.1X'])/100
#     }
#   }
#   nucl.dist = nucl.dist * coverage.correction
#   write.table(nucl.dist,paste0("genetic_final.1/pi.full.genome/",i,".genetic.distance.tab"),sep="\t")
# }
