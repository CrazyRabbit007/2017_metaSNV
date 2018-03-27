#setwd('/g/bork5/costea/Strain_bigger/SNP_data')
#setwd('/home/costea/g5/Strain_bigger/SNP_data')

source('util.R')

#This is an fpc package function, which needed a bit of fixing. 
classifdist <- function (cdist, clustering, method = "averagedist", centroids = NULL, 
    nnk = 1) 
{
    k <- max(clustering)
    n <- nrow(data)
    cdist <- as.matrix(cdist)
    topredict <- clustering < 0
    if (method == "averagedist") {
        prmatrix <- matrix(0, ncol = k, nrow = sum(topredict))
        for (j in 1:k) prmatrix[, j] <- rowMeans(as.matrix(cdist[topredict,clustering == j]))
        clpred <- apply(prmatrix, 1, which.min)
        clustering[topredict] <- clpred
    }
    if (method == "centroid") 
        clustering[topredict] <- apply(cdist[topredict, centroids,drop=F], 
            1, which.min)
    if (method == "knn") {
        cdist[topredict, topredict] <- max(cdist) + 1
        if (nnk == 1) {
            bestobs <- apply(cdist[topredict, ], 1, which.min)
            clustering[topredict] <- clustering[bestobs]
        }
        else {
            for (i in (1:n)[topredict]) {
                bestobs <- order(cdist[i, ])[1:k]
                clasum <- numeric(0)
                for (j in 1:k) clasum[j] <- sum(clustering[bestobs] == 
                  j)
                clustering[i] <- which.max(clasum)
            }
        }
    }
    clustering
}

cleanTimepoints <- function(data) {
  removed <- NULL
  #Only the American and the German samples have multiple time-points. So, clean that up
  americans <- grep('American',as.character(meta[colnames(data),,drop=F]$Nationality))
  am_names <- colnames(data)[americans]
  am_indiv <- sapply(strsplit(am_names,'-'),'[[',1)
  
  if (anyDuplicated(am_indiv)) {
    data <- data[-americans[which(duplicated(am_indiv))],-americans[which(duplicated(am_indiv))]]
    removed <- c(removed,am_names[which(duplicated(am_indiv))])
  }
  
  #Now handle the germans
  germans <- grep('German',as.character(meta[colnames(data),,drop=F]$Nationality))
  gr_names <- colnames(data)[germans]
  gr_indiv <- sapply(strsplit(gr_names,'-'),'[[',1)
  if (anyDuplicated(gr_indiv)) {
    data <- data[-germans[which(duplicated(gr_indiv))],-germans[which(duplicated(gr_indiv))]]
    removed <- c(removed,gr_names[which(duplicated(gr_indiv))])
  }
  
  #Now handle the Kazakhstan
  kazak <- grep('Kazakhstan',as.character(meta[colnames(data),,drop=F]$Nationality))
  kz_names <- colnames(data)[kazak]
  kz_indiv <- sapply(strsplit(kz_names,'-'),'[[',1)
  if (anyDuplicated(kz_indiv)) {
    data <- data[-kazak[which(duplicated(kz_indiv))],-kazak[which(duplicated(kz_indiv))]]
    removed <- c(removed,kz_names[which(duplicated(kz_indiv))])
  }
  
  return(list(data,removed))
}

#Get prediction strength. This is modified to work directly on the distance matrix
pred <- function (distance, Gmin = 2, Gmax = 10, M = 50, 
    classification = "centroid", cutoff = 0.75, nnk = 1, ...) 
{
    require(cluster)
    require(class)
    dist <- as.matrix(distance)
    n <- nrow(dist)
    nf <- c(floor(n*0.5), n - floor(n*0.5))
    indvec <- clcenters <- clusterings <- jclusterings <- classifications <- list()
    prederr <- list()
    
    for (k in Gmin:Gmax) {
        prederr[[k]] <- numeric(0)
        for (l in 1:M) {
            nperm <- sample(n, n)
            indvec[[l]] <- list()
            indvec[[l]][[1]] <- nperm[1:nf[1]]
            indvec[[l]][[2]] <- nperm[(nf[1] + 1):n]
            for (i in 1:2) {
                clusterings[[i]] <- as.vector(pam(as.dist(dist[indvec[[l]][[i]],indvec[[l]][[i]]]), k, diss=TRUE))
                jclusterings[[i]] <- rep(-1, n)
                jclusterings[[i]][indvec[[l]][[i]]] <- clusterings[[i]]$clustering
		centroids <- clusterings[[i]]$medoids
                j <- 3 - i
                classifications[[j]] <- classifdist(distance, jclusterings[[i]], 
                  method = classification, centroids = centroids, 
                  nnk = nnk)[indvec[[l]][[j]]]
            }
            
            ps_f <- matrix(0, nrow = 2, ncol = k)
            
            for (i in 1:2) {
                for (kk in 1:k) {
                  nik <- sum(clusterings[[i]]$clustering == kk)
                  if (nik > 1) {
		         a <- which(clusterings[[i]]$clustering[1:(nf[i] - 1)] == kk)
			       ps_f[i,kk] <- sum(outer(classifications[[i]][a],classifications[[i]][a],'=='))-length(a)
                   ps_f[i,kk] <- ps_f[i, kk]/(nik * (nik - 1))
                  }
                }
            }
            
            ps <- ps_f
            
            prederr[[k]][l] <- mean(c(min(ps[1, ]), min(ps[2, 
                ])))
        }
    }
    mean.pred <- numeric(0)
    if (Gmin > 1) 
        mean.pred <- c(1)
    if (Gmin > 2) 
        mean.pred <- c(mean.pred, rep(NA, Gmin - 2))
    for (k in Gmin:Gmax) mean.pred <- c(mean.pred, mean(prederr[[k]]))
    optimalk <- max(which(mean.pred > cutoff))
    out <- list(predcorr = prederr, mean.pred = mean.pred, optimalk = optimalk, 
        cutoff = cutoff, method = clusterings[[1]]$clustermethod, 
        Gmax = Gmax, M = M)
    class(out) <- "predstr"
    out
}

library(fpc)
library(ape)
source('Strain_functions.R')

meta <- read.table('All_samples_metadata',header=T,row.names=1,sep='\t')

args <- commandArgs(trailingOnly = TRUE)
#Get file to process
file <- args[1]
file_name <- strsplit(file,"/")[[1]][2]

dist <- read.table(file,header=T,row.names=1,check.names=F)
clean <- cleanTimepoints(dist) 

if (nrow(dist) < 50) {
	quit()
}
#Get predictions strength
set.seed(135649685)
res <- pred(dist)
#Get bootstrapped k-means clustering for optimal k
if (res[["optimalk"]] == 1) {#There's no clustering here
  df <- data.frame(row.names=colnames(dist),clust=rep(1,ncol(dist)))
} else {
  #Get clustering for the non-replicate containing set
  clustering <- pam(dist, res[["optimalk"]], diss=TRUE)
  #Now, get cluster assignments for all the samples that were left out
  assign <- rep(-1,ncol(dist))
  names(assign) <- colnames(dist)
  assign[names(clustering$clustering)] <- clustering$clustering
  centroids <- clustering$medoids
  #Compute centroid based assignment
  clust <- classifdist(dist, assign, method = "centroid", centroids = centroids)
  
  df <- data.frame(clust)
}

write.table(res[["mean.pred"]],paste('clust/',file_name,'.ps',sep=''),sep='\t',quote=F)
write.table(df,paste('clust/',file_name,'_clustering.tab',sep=''),sep='\t',quote=F)

#Get PCOA projection
#pca <- pcoa(dist)
#eig <- pca[["values"]][["Eigenvalues"]]
#eig[eig<0] <- 0
#eig <- eig/sum(eig)*100

#pcoa_df <- data.frame(pca[["vectors"]][,1:3])
#pcoa_df$f <- freq[rownames(pcoa_df),]$freq_data_sample_10_90
#pcoa_df$clust <- df[rownames(pcoa_df),]

#write.table(pcoa_df,paste('clust/',species,'_mann_pcoa_proj.tab',sep=''),sep='\t',quote=F)

