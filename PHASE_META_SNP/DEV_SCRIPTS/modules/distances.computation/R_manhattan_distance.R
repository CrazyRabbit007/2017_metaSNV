# Rscript aiming to compare the distance calculation of metaSNP and R classical function

i = "1002672"

snp.freq = read.delim(paste0("/science/paolil/TARA_metaSNP/tara_4_metaSNP/filtered_universal/pop/",i,".core.filtered.freq"), row.names=1)
names(snp.freq) = gsub("X","",names(snp.freq))

head(snp.freq)

snp.freq.all=snp.freq

snp.freq.0=snp.freq
snp.freq.Ø=snp.freq.0[snp.freq.0==-1]=0

snp.freq[snp.freq==-1]=NA

head(snp.freq)

man.dist.all=as.data.frame(as.matrix(dist(t(snp.freq.all),method="manhattan",upper=T,diag=T)))
man.dist.na=as.data.frame(as.matrix(dist(t(snp.freq),method="manhattan",upper=T,diag=T)))
man.dist.0=as.data.frame(as.matrix(dist(t(snp.freq.0),method="manhattan",upper=T,diag=T)))
names(man.dis.all)=names(snp.freq)
names(man.dis.na)=names(snp.freq)
names(man.dis.0)=names(snp.freq)
row.names(man.dis.all)=names(snp.freq)
row.names(man.dis.na)=names(snp.freq)
row.names(man.dis.0)=names(snp.freq)

head(man.dist.all)
head(man.dist.na)
head(man.dist.0)

write.table(man.dist.all,file=paste0(i,".all.mann.R.dist.tab"),row.names=T,quote=F,sep="\t")
write.table(man.dist.na,file=paste0(i,".na.mann.R.dist.tab"),row.names=T,quote=F,sep="\t")
write.table(man.dist.0,file=paste0(i,".0.mann.R.dist.tab"),row.names=T,quote=F,sep="\t")


