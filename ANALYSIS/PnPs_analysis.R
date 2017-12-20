
setwd('/Users/Lucas/Documents/M2_IMaLiS-EBE/STAGE/Paper_metaSNV')

pnps = read.delim('summary.pnps')
names(pnps) = gsub('TARA_|.filtered.screened.tara.adapters.on.freeze11.marine.repgenomes.solexaqa.allbest.l45.p97.unique.sorted.bam','',names(pnps))
names(pnps) = gsub('_0.22.3_|_0.22.1.6_','_',names(pnps))
names(pnps)

names = unique(gsub('_G|_T','',names(pnps)))

proper.pnps = data.frame(matrix(ncol = length(names),nrow=nrow(pnps)))
names(proper.pnps) = names
proper.pnps$X = pnps$X

for (i in names[2:length(names)]){
  if (paste0(i,'_T') %in% names(pnps) & paste0(i,'_G') %in% names(pnps)){
    proper.pnps[,i] = pnps[,paste0(i,'_T')]/pnps[,paste0(i,'_G')]
  }
}

max(proper.pnps[,2:ncol(proper.pnps)],na.rm=T)
min(proper.pnps[,2:ncol(proper.pnps)],na.rm=T)
mean(proper.pnps[,2:ncol(proper.pnps)],na.rm=T)

