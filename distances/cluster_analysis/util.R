##### Keep all "common" functions in here.
##### If something is defined here, it'd better not be anywhere else.

#Remove duplicate time-points
cleanTimepoints <- function(names,meta) {#Returns a vector on unique time-points

  americans <- grep('American',as.character(meta[names,,drop=F]$Nationality))
  am_names <- names[americans]
  am_indiv <- sapply(strsplit(am_names,'-'),'[[',1)
  if (anyDuplicated(am_indiv)) {
    names <- names[-americans[which(duplicated(am_indiv))]]
  }
  
  #Now handle the germans
  germans <- grep('German',as.character(meta[names,,drop=F]$Nationality))
  gr_names <- names[germans]
  gr_indiv <- sapply(strsplit(gr_names,'-'),'[[',1)
  if (anyDuplicated(gr_indiv)) {
    names <- names[-germans[which(duplicated(gr_indiv))]]
  }
  
  #Now handle the Kazakhstan
  kazak <- grep('Kazakhstan',as.character(meta[names,,drop=F]$Nationality))
  kz_names <- names[kazak]
  kz_indiv <- sapply(strsplit(kz_names,'-'),'[[',1)
  if (anyDuplicated(kz_indiv)) {
    names <- names[-kazak[which(duplicated(kz_indiv))]]
  }
  
  return(names)
}

getSamplesWithMultipleTimepoints <- function(names,meta) {

  multi <- NULL

  americans <- grep('American',as.character(meta[names,,drop=F]$Nationality))
  if (length(americans) > 0) {
    am_names <- names[americans]
    am_indiv <- sapply(strsplit(am_names,'-'),'[[',1)
    #Get duplicated
    multiple <- names(which(table(am_indiv)>1))
    if (length(multiple) > 0) {
      multi <- c(multi,names[americans[which(am_indiv %in% multiple)]])
    }
  }
  #Now handle the germans
  germans <- grep('German',as.character(meta[names,,drop=F]$Nationality))
  if (length(germans) > 0) {
    gr_names <- names[germans]
    gr_indiv <- sapply(strsplit(gr_names,'-'),'[[',1)
    #Get duplicated
    multiple <- names(which(table(gr_indiv)>1))
    if (length(multiple) > 0) {
      multi <- c(multi,names[germans[which(gr_indiv %in% multiple)]])
    }
  }
  #Now handle the Kazakhstan
  kazak <- grep('Kazakhstan',as.character(meta[names,,drop=F]$Nationality))
  if (length(kazak) > 0) {
    kz_names <- names[kazak]
    kz_indiv <- sapply(strsplit(kz_names,'-'),'[[',1)
    #Get duplicated
    multiple <- names(which(table(kz_indiv)>1))
    if (length(multiple) > 0) {
      multi <- c(multi,names[kazak[which(kz_indiv %in% multiple)]])
    }
  }
  return(multi)
}
