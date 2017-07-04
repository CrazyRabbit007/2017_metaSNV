
######################### Distances ############################

allel.dist <- function(inMatrix, ...) {
	 matrixColSize <- length(colnames(inMatrix))
	 matrixRowSize <- length(rownames(inMatrix))
	 colnames <- colnames(inMatrix)
	 resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
       
        dyn.load(paste(dist_location,'freq_euc.so',sep=''))
        r = .C('allel_dist',as.double(as.matrix(inMatrix)),as.integer(dim(inMatrix)),as.double(resultsMatrix),as.integer(60))
	resultsMatrix <- t(matrix(r[[3]],matrixColSize,matrixColSize))
	
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
}

mann.dist <- function(inMatrix,...) {
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
        dyn.load(paste(dist_location,'freq_euc.so',sep=''))
        r = .C('manhattan_dist',as.double(as.matrix(inMatrix)),as.integer(dim(inMatrix)),as.double(resultsMatrix))
	resultsMatrix <- t(matrix(r[[3]],matrixColSize,matrixColSize))
	
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
  
}

na.dist <- function(inMatrix,...) {
	matrixColSize <- length(colnames(inMatrix))
	matrixRowSize <- length(rownames(inMatrix))
	colnames <- colnames(inMatrix)
	resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
        dyn.load(paste(dist_location,'freq_euc.so',sep=''))
        r = .C('NA_dist',as.double(as.matrix(inMatrix)),as.integer(dim(inMatrix)),as.double(resultsMatrix))
	resultsMatrix <- t(matrix(r[[3]],matrixColSize,matrixColSize))
	
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
  
}

euc.dist <- function(inMatrix, ...) {
	 matrixColSize <- length(colnames(inMatrix))
	 matrixRowSize <- length(rownames(inMatrix))
	 colnames <- colnames(inMatrix)
	 resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
        
        dyn.load(paste(dist_location,'freq_euc.so',sep=''))
        r = .C('euc_dist',as.double(as.matrix(inMatrix)),as.integer(dim(inMatrix)),as.double(resultsMatrix))
	resultsMatrix <- t(matrix(r[[3]],matrixColSize,matrixColSize))
	
	colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
	as.dist(resultsMatrix)->resultsMatrix
	attr(resultsMatrix, "method") <- "dist"
	return(resultsMatrix) 
}

kNeighborDistance <-function(dist1,dist2,k=10) {
  dist1 <- as.matrix(dist1)
  dist2 <- as.matrix(dist2)
  
  if (nrow(dist1) != nrow(dist2)) {
    print('Cannot compare distances on different sets')
    return(-1)
  }
  
  same = 0
  for (i in 1:nrow(dist1)) {
    a1 <- dist1[i,-i]
    a2 <- dist2[i,-i]
    a1 <- order(a1)[1:k]
    a2 <- order(a2)[1:k]
    same = same + length(which(a1 %in% a2))
  }
  
  same <- same/(nrow(dist1)*k)
  return(same)
  
}


################################################################
######################### Clustering ###########################

checkMetaHitEnrichment <- function(strain) {
  pvals <- NULL
  features <- NULL
  plots <- NULL
  strainMeta$strain <- strain[rownames(strainMeta),]
  
  for (i in 1:(ncol(strainMeta)-1)) {
    set1 <- subset(strainMeta,strain==1)[,i]
    set2 <- subset(strainMeta,strain==2)[,i]
    if (is.numeric(set1)) {
      w <- wilcox.test(set1,set2)
      if (!is.na(w[["p.value"]])) {
	pvals <- c(pvals,w[["p.value"]])
	features <- c(features,colnames(strainMeta)[i])
      }
      if (!is.na(w[["p.value"]]) & w[["p.value"]] < 0.05) {
	  print (colnames(strainMeta)[i])
	  plots <- c(plots,list(ggplot(strainMeta,aes_string(x="strain",y=colnames(strainMeta)[i])) + geom_boxplot() + ggtitle(as.character(w[["p.value"]]))))
	}
    } else {
      m <- subset(strainMeta,!is.na(strain))
      m <- subset(m,!is.na(m[,i]))
      
      if (nrow(m) > 1) {
	print(fisher.test(m$strain,m[,i]))
	print(paste('Fisher tested: ',colnames(m)[i]))
	print(m$strain)
	print(m[,i])
      }
    }
  }

  pvals <- p.adjust(pvals,method='fdr')
  #print(features)
  #print(pvals)
}

getClustering_tree <- function(euc_dist) {
  h <- hclust(euc_dist)
  l <- NULL
  for (k in 2:10) {
    c <- cutree(h,k)
    if (min(table(c))<=5) {
      break
    }
    for (i in 1:max(c)) {
      fDist <- (abs( rowMeans(data[,which(c==i)],na.rm=T) - rowMeans(data[,which(c!=i)],na.rm=T) ) )
      oList <- which(fDist >= 60)
      l <- rbind(l,c(length(oList),i,k))
    }
  }
  if (is.null(l)) {
    c <- rep(1,ncol(data))
  } else if (nrow(l) < 1) {#The is a major outlier! Skip for now
    c <- rep(1,ncol(data))
  } else {
    l <- data.frame(l)
    a <- aggregate(X1 ~ X3,l,mean)$X1  
    #Get clustering
    #cluster <- pamk(euc_dist,krange=1:10,criterion="asw")
    c <- cutree(h,(which.max(a)+1))
  }
  return(c)
}

library(cluster)

#Cluster the distance
pam.clustering=function(x,k) { # x is a distance matrix and k the number of cluster
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE))
  return(cluster)
}

getClustering_kMeans1 <- function(data,euc_dist) {
      
  ddd <- data
  #Total sum of squares
  sst <- apply(ddd,1,function(x) {x <- x[x!=-1]; sum((x-mean(x))^2)})  
  
  #Replace missing with mean
  ddd1 <- t(apply(ddd,1,function(x) {x1 <- x[x!=-1]; m <- mean(x1); x[x==-1] <- m; return(x)}))
  
  nsilhouette=NULL
  var = NULL

  nsilhouetteK=NULL
  varK = NULL
  
  for (k in 1:10) { 
      if (k==1) {
	nsilhouette[k]=0
	var[k] = 0
      }

      else {
	cluster=pam.clustering(euc_dist, k)$clustering
	nsilhouette[k]=mean(silhouette(cluster, euc_dist)[,3])
	#Sum of squares between
	ssb <- rep(0,length(sst))
	classes = unique(cluster)
	for (cl in classes) {
	  d1 <- ddd1[,which(cluster==cl),drop=FALSE]
	  ssb <- ssb + sum(cluster==cl)*((apply(d1,1,mean)-apply(ddd1,1,mean))^2)
	}
	var[k] = sum(ssb)/sum(sst)*100
	print(var)
	print(nsilhouette)
	
	cluster=kmeans(euc_dist,k)$cluster
	nsilhouetteK[k]=mean(silhouette(cluster, euc_dist)[,3])
	#Sum of squares between
	ssb <- rep(0,length(sst))
	classes = unique(cluster)
	for (cl in classes) {
	  d1 <- ddd1[,which(cluster==cl),drop=FALSE]
	  ssb <- ssb + sum(cluster==cl)*((apply(d1,1,mean)-apply(ddd1,1,mean))^2)
	}
	varK[k] = sum(ssb)/sum(sst)*100
	
      }
  }
  
  df <- data.frame(sil=nsilhouette,var=var,silK=nsilhouetteK,varK=varK)
  write.table(df,paste(species,'_variation_explained.tab',sep=''),quote=F,row.names=F,sep='\t')
}

getClustering_Pamk <- function(euc_dist) {
  #Get clustering
  cluster <- pamk(euc_dist,krange=1:10,criterion="asw")
  c <- cluster[["pamobject"]][["clustering"]]
  return(list(c,cluster$crit[cluster$nc]))
}

library(clue)

getClustering_kMeans_boot <- function(euc_dist) {
  sil = c(0)
  for (k in 2:10) {
    bb <- cl_boot(euc_dist,50,k)
    c1 <- cl_consensus(bb)
    cluster <- as.numeric(cl_class_ids(c1))
    if (length(unique(cluster)) < k) {#The consensus is that this clustering is terrible
      sil <- c(sil,0)
    } else {
      sil <- c(sil,mean(silhouette(cluster, euc_dist)[,3]))
    }
  }
  print(sil)
  k <- which.max(sil)
  bb <- cl_boot(euc_dist,50,k)
  c1 <- cl_consensus(bb)
  return(list(as.numeric(cl_class_ids(c1)),sil[k]))
}

################################################################
######################### Plots ################################

Family_plot <- function(dist,df) {
  dist <- as.matrix(dist)
  g <- grep('genome',colnames(dist))
  if (length(g) > 0) {
    dist <- dist[-g,-g]
  }
  #Remove donald1-11-30-0
  g <- grep('donald1-11-30-0',colnames(dist))
  if (length(g) > 0) {
    dist <- dist[-g,-g]
  }
  
  dist2 <- melt(dist)[melt(upper.tri(dist))$value,]
  
  #Fix the labels which we know to be flipped:
  # mickey3-11-0-0 with mickey5-11-30-0 
  
  r1 <- which(dist2$X1=='mickey3-11-0-0')
  r2 <- which(dist2$X1=='mickey5-11-30-0')
  if (length(r1) > 0) {
    dist2$X1[r1] <- 'mickey5-11-30-0'
  }
  if (length(r2) > 0) {
    dist2$X1[r2] <- 'mickey3-11-0-0'
  }
  
  r1 <- which(dist2$X2=='mickey3-11-0-0')
  r2 <- which(dist2$X2=='mickey5-11-30-0')
  if (length(r1) > 0) {
    dist2$X2[r1] <- 'mickey5-11-30-0'
  }
  if (length(r2) > 0) {
    dist2$X2[r2] <- 'mickey3-11-0-0'
  }
  
  
  dist2$Fam1 <- fam_data[as.character(dist2$X1),]$family
  dist2$Fam2 <- fam_data[as.character(dist2$X2),]$family
  dist2$Sample1 <- df[as.character(dist2$X1),]$sample
  dist2$Sample2 <- df[as.character(dist2$X2),]$sample
  
  dist2$Relation <- 'Strangers'
  dist2$Relation[which(dist2[["Fam1"]]==dist2[["Fam2"]])] <- 'Family'
  dist2$Relation[which(dist2[["Sample1"]]==dist2[["Sample2"]])] <- 'Individual'
  
  fam_set <- subset(dist2,Relation=='Family')
  fam_count <- table(fam_set[["Fam1"]])
  relevant_fam <- names(which(fam_count > 4))
  
  plots <- list()
  all_df <- NULL
  
  for (f in relevant_fam) {#We're in business
    #Get individual distances, keep the replicates
    indiv <- subset(dist2, Fam1==f & Relation=='Individual')
    fSet <- subset(dist2,(Fam1==f | Fam2==f) & Relation=='Family')
    fSet <- cbind(c(as.character(fSet$X1),as.character(fSet$X2)),c(fSet$Sample1,fSet$Sample2))
    colnames(fSet) <- c('X1','Sample')
    toKeep <- aggregate(X1 ~ Sample,fSet,function(x) as.character(x)[1])
    
    ff <- subset(dist2,(Fam1==f | Fam2==f) & Relation=='Family')
    #Remove from dist2 anything that is samples1 but not first replicate
    w <- which((ff[["X1"]] %in% toKeep[["X1"]] & (ff[["X2"]] %in% toKeep[["X1"]])))
    indiv <- rbind(indiv,ff[w,])
    
    ff <- subset(dist2,(Fam1==f | Fam2==f) & Relation=='Strangers')
    #Remove from dist2 anything that is samples1 but not first replicate
    #w <- which((ff[["X1"]] %in% toKeep[["X1"]] | (ff[["X2"]] %in% toKeep[["X1"]])))
    indiv <- rbind(indiv,ff)
    indiv[["Relation"]] <- factor(indiv[["Relation"]],levels=c("Individual","Family","Strangers"),ordered=T)
    
    s <- max(indiv[["value"]])
    lines_df <- structure(list(x = c(1, 1, 2, 2, 2, 3, 3, 3, 1), y = c(s, s+2, s, s+5, s+ 7, s+ 5, s+10, s+12, s+10), xend = c(1, 2, 2, 2, 3 ,3, 3, 1, 1), yend = c(s+2, s+2, s+2, s+7, s+7, s+7, s+12, s+12, s+12)), .Names = c("x", "y", "xend", "yend"),
            row.names = c(NA, -9L), class = "data.frame")
    astpos_df <- structure(list(x = c(1.5, 2.5, 2), y = c(s+3, s+8, s+13)), .Names = c("x", "y"), row.names = c(NA, -3L), class = "data.frame" )
    
    w <- pairwise.wilcox.test(indiv[["value"]],indiv[["Relation"]])[['p.value']]
    a <- melt(w)[melt(lower.tri(w,diag=T))$value,]
    
    astpos_df$value <- a$value[c(1,3,2)]
    indiv <- indiv[order(indiv[["value"]]),]
    all_df <- rbind(all_df,cbind(indiv,f))
    g <- ggplot(indiv,aes(x=Relation,y=value)) + geom_boxplot() + geom_point() + ggtitle(paste('Family',f)) + geom_segment(data = lines_df, size = .5, aes(x=x, y=y, xend=xend, yend=yend)) + geom_text(data = astpos_df, aes(x=x, y=y,label=sprintf('%3.6f',value)), size = 4)
    plots <- c(plots,list(g))
  }
  
  if (length(relevant_fam)>0) {
    do.call(grid.arrange,plots)
    write.table(all_df,paste(species,'_family_cons.tab',sep=''),quote=F,row.names=F,sep='\t')
  }
  
}

Replicate_plot <- function(dist,df) {
  dist <- as.matrix(dist)
  dist2 <- melt(dist)
  dist2 <- dist2[order(dist2$X2,dist2$value),]
  
  withReplicates <- names(which(table(df[["sample"]])>1))
  samples_to_keep <- subset(df,sample %in% withReplicates)
  
  dist2 <- subset(dist2,X2 %in% rownames(samples_to_keep))
  prox <- NULL
  for (s in rownames(samples_to_keep)) {
    dd <- subset(dist2,X2==s)[-1,]
    closest_self <- min(which(df[as.character(dd$X1),]$sample == samples_to_keep[s,]$sample))
    closest_non_self <- min(which(df[as.character(dd$X1),]$sample != samples_to_keep[s,]$sample))
    prox <- rbind(prox,c(s,dd[["value"]][closest_self],'self'))
    prox <- rbind(prox,c(s,dd[["value"]][closest_non_self],'non-self'))
  }
  prox <- data.frame(prox)
  prox$X2 <- as.numeric(as.character(prox$X2))
  g <- ggplot(prox,aes(x=X1,y=X2,color=X3,group=X3)) + geom_point(size=3) + ylab('Distance') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_line()
  return(g)
}

Replicate_plot_US_only <- function(dist,df,spec) {
  dist <- as.matrix(dist)
  s <- which(rownames(dist) %in% samples_252)
  dist <- dist[s,s]
  
  dist2 <- melt(dist)
  dist2 <- dist2[order(dist2$X2,dist2$value),]
  
  withReplicates <- names(which(table(df[["sample"]])>1))
  samples_to_keep <- subset(df,sample %in% withReplicates)
  
  keep <- grep('stool',rownames(samples_to_keep))
  if (length(keep)==0) {
    return()
  }
  samples_to_keep <- samples_to_keep[keep,]
  
  dist2 <- subset(dist2,X2 %in% rownames(samples_to_keep))
  prox <- NULL
  for (s in rownames(samples_to_keep)) {
    dd <- subset(dist2,X2==s)[-1,]
    closest_self <- min(which(df[as.character(dd$X1),]$sample == samples_to_keep[s,]$sample))
    closest_non_self <- min(which(df[as.character(dd$X1),]$sample != samples_to_keep[s,]$sample))
    prox <- rbind(prox,c(s,dd[["value"]][closest_self],'self'))
    prox <- rbind(prox,c(s,dd[["value"]][closest_non_self],'non-self'))
  }
  prox <- data.frame(prox)
  prox$X2 <- as.numeric(as.character(prox$X2))
  write.table(prox,paste(spec,'_Timepoint_proximity_252.tab',sep=''),quote=F,row.names=F,col.names=F)
  g <- ggplot(prox,aes(x=X1,y=X2,color=X3,group=X3)) + geom_point(size=3) + ylab('Distance') + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + geom_line()
  return(g)
}

Replicate_plot_Comp <- function(dist,nDist,df) {
  dist <- as.matrix(dist)
  dist2 <- melt(dist)[melt(upper.tri(dist))$value,]
  o <- order(dist2$X2,dist2$value)
  dist2 <- dist2[o,]
  
  nDist <- as.matrix(nDist)
  ndist <- melt(nDist)[melt(upper.tri(nDist))$value,]
  ndist <- ndist[o,]
  
  all <- cbind(dist2,ndist[,3])
  colnames(all) <- c('S1','S2','Dist','Bases.compared')
  all$same_sample <- df[as.character(all$S1),]$sample==df[as.character(all$S2),]$sample
 
  non_same <- which(all$same_sample==FALSE)
  same <- which(all$same_sample==TRUE)
  all <- all[c(same,sample(non_same,100)),]
  all$sample <- df[as.character(all$S1),]$sample
  all$sample[which(all$same_sample==FALSE)] <- NA
 
  g <- ggplot(all,aes(x=Dist,y=log10(Bases.compared),color=same_sample)) + geom_point() + scale_color_brewer(palette='Set1') + geom_text(aes(label=paste(S1,S2,sep='\n')))
  return(g)
}

NN_plot <- function(dist,df) {
  dist <- as.matrix(dist)
  dist2 <- melt(dist)
  dist2 <- dist2[order(dist2$X2,dist2$value),]
  #Keep only 10 for each "thing"
  k = 1:10
  l = dim(dist)[1]
  for (i in 1:l) {
    k = c(k,(l*i+1:10))
  }
  dist2 <- dist2[k,]
  dist2$from <- df[as.character(dist2$X1),]$from
  dist2$order <- 1:10
  
  g <- ggplot(dist2,aes(x=order,y=as.factor(X2),fill=from)) + geom_tile() + scale_fill_manual(values=col)
  return(g)
}


