#Get split for ref mOTUs

dict = read.table('map_genomes_2_ref-motus.tsv',header=F)
motu <- list.files('clust','ref.*mann.*tab')

nrClust = NULL

for (mSplit in motu) {#Ok, we have some ref based split.
  #Find relevant specI
  motu_number = strsplit(mSplit,'\\.')[[1]][1]
  specs = as.character(subset(dict,as.character(V1)==motu_number)$V2)
  specs = sapply(strsplit(specs,'\\.'),'[[',1)
  ff = NULL
  for (s in specs) {
    f = paste0('clust/',s,'.filtered.mann.dist_clustering.tab')
    if (file.exists(f)) {
      ff <- c(ff,f)
    }
  }
  if (length(ff)==1) {
    mm = read.table(paste0('clust/',mSplit),header=T,row.names=1,check.names=F)
    ss = read.table(ff[1],header=T,row.names=1,check.names=F)
    rownames(mm) <- gsub(".snv.bam",'',rownames(mm))
    rownames(ss) <- gsub(".freeze11.snv.mod.fix.bam",'',rownames(ss))
    common = rownames(mm)[which(rownames(mm) %in% rownames(ss))]
    mmc <- mm[common,]
    ssc <- ss[common,]
    nrClust = rbind(nrClust,c(motu_number,max(mm$clust),ff[1],max(ss$clust),sum(mmc==ssc)/length(mmc)))
  }
}

nrClust <- data.frame(nrClust,stringsAsFactors=F)
all <- nrClust
nrClust <- subset(nrClust,X2!=1 | X4 != 1)

correct <- subset(nrClust,X2==X4)
print(dim(correct))
wrong <- subset(nrClust,X2!=X4)

dim(subset(wrong,X5<0.9))
