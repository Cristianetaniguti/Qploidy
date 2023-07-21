#' Edit Axiom Summary file
#'
#' Remove consecutive A allele probe from the summary file
#'
#' @param summary_df data.frame with A and B probes intensities
#'
#' @return list with A and B probes
#'
#' @export
clean_summary <- function(summary_df){
  summaries_filt <- summary_df
  if(length(grep("Contig", summary_df$probeset_id)) > 0) summaries_filt <- summary_df[-grep("Contig", summaries_filt$probeset_id),]
  if(length(grep("contig", summaries_filt$probeset_id)) > 0) summaries_filt <- summaries_filt[-grep("contig", summaries_filt$probeset_id),]
  if(length(grep("comp", summaries_filt$probeset_id)) > 0) summaries_filt <- summaries_filt[-grep("comp", summaries_filt$probeset_id),]

  A_probes <- summaries_filt[seq(1,dim(summaries_filt)[1], 2),]
  B_probes <- summaries_filt[seq(2,dim(summaries_filt)[1], 2),]
  return(list(A_probes=A_probes, B_probes = B_probes))
}

#' Get R and theta values from summary file with option to do standard normalization by plate and markers
#'
#' @param cleaned_summary summary object from clean_summary function. List with A and B intensities.
#' @param ind.names data.frame with individual names correspondence
#' @param sd.normalization logical. If TRUE performs standard normalization
#' @param atan logical. If TRUE calculates the theta using atan2
#'
#' @import tidyr
#' @importFrom dplyr mutate
#'
#' @export
get_R_theta <- function(cleaned_summary, ind.names, sd.normalization = FALSE, atan = FALSE){
  R_all <- as.data.frame(cleaned_summary$A_probes[,-1] + cleaned_summary$B_probes[,-1])
  if(atan){
    theta_all <- as.data.frame((atan2(as.matrix(cleaned_summary$B_probes[,-1]), as.matrix(cleaned_summary$A_probes[,-1])))/(pi/2))
  }  else {
    theta_all <- as.data.frame(cleaned_summary$B_probes[,-1]/(cleaned_summary$B_probes[,-1] + cleaned_summary$A_probes[,-1]))
  }

  probes_names <-  cleaned_summary$A_probes[,1]
  probes_names <- gsub("-A","", as.vector(probes_names$probeset_id))

  R_all <- cbind(MarkerName = probes_names, R_all)
  theta_all <- cbind(MarkerName = probes_names, theta_all)

  array.id <- colnames(R_all)
  ind.names <- as.data.frame(ind.names)
  ids <- ind.names[,2][match(colnames(R_all), ind.names[,1])]
  colnames(R_all) <- ids

  ids <- ind.names[,2][match(colnames(theta_all), ind.names[,1])]
  colnames(theta_all) <- ids

  colnames(theta_all)[1] <- colnames(R_all)[1] <- "MarkerName"

  # R Quantile normalization
  if(sd.normalization){
    R_melt <- pivot_longer(R_all, cols = -1, names_to = "ind", values_to = "R")

    exp <- sapply(strsplit(array.id[-1], "-"), "[[", 4)
    exp <- sapply(strsplit(exp, "_"), "[[", 1)

    array_ids <- data.frame(array=exp, inds = ids[-1])
    R_melt$array <- array_ids$array[match(R_melt$ind, array_ids$inds)]

    R_z <- R_melt %>% group_by(MarkerName, array) %>% mutate(zscore = (R-mean(R))/sd(R))
    R_all <- pivot_wider(R_z[,c(1,2,5)], names_from = ind, values_from = zscore)
  }

  return(list(R_all=R_all, theta_all = theta_all))
}

#' Convert Summary file to fitpoly format
#'
#' @param R_all data.frame with R values
#' @param theta_all data.frame with theta values
#' @param ind.names data.frame first column with the plate name and second with sample name
#'
#' @importFrom tidyr pivot_longer
#' @export
summary_to_fitpoly <- function(R_all, theta_all){
  Y_all <- R_all[,-1]*theta_all[,-1]
  X_all <- R_all[,-1] - Y_all

  X_all <- cbind(MarkerName= R_all[,1], X_all)
  Y_all <- cbind(MarkerName = R_all[,1], Y_all)

  R <- pivot_longer(R_all, cols = 2:ncol(R_all), values_to = "R", names_to = "SampleName")
  theta <- pivot_longer(theta_all, cols = 2:ncol(theta_all), values_to = "ratio", names_to = "SampleName")
  X <- pivot_longer(X_all, cols = 2:ncol(X_all), values_to = "X", names_to = "SampleName")
  Y <- pivot_longer(Y_all, cols = 2:ncol(Y_all), values_to = "Y", names_to = "SampleName")

  fitpoly_input <- cbind(X, Y = Y[,3], R = R[,3], ratio = theta[,3])

  return(fitpoly_input)
}

#' Extracts genotype probabilities information from polyOrigin polyancentry output file
#' Code adapted from diaQTL
#'
#' @param filename character with path to polyancentry CSV file
#'
#' @export
read_polyOrigin <- function(filename){
  con <- file(filename,"r")
  temp <- readLines(con)
  ix <- grep("PolyOrigin",temp)

  k <- grep("offspringinfo",temp[ix])
  x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
  ped <- data.frame(id=sapply(x,function(x){x[1]}),pop=sapply(x,function(x){x[2]}),stringsAsFactors = F)

  k <- grep("ancestralgenotype",temp[ix])
  x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
  parents <- apply(cbind(sapply(x,function(x){x[1]}),sapply(x,function(x){x[3]})),1,paste,collapse="|")
  parents <- unique(parents)
  y <- strsplit(parents,split="|",fixed=T)
  parents <- data.frame(pop=sapply(y,function(y){y[1]}),parent1=sapply(y,function(y){y[2]}),parent2=sapply(y,function(y){y[length(y)]}),stringsAsFactors = F)
  ped <- merge(ped,parents)

  k <- grep("parentgeno",temp[ix])
  header <- setdiff(strsplit(temp[ix[k]+1],split=",",fixed=T)[[1]],"")
  parents <- header[-(1:3)]
  header[1:3] <- c("marker","chrom","cM")
  x <- strsplit(temp[c((ix[k]+2):(ix[k+1]-1))],split=",")
  m <- length(x)
  pp <- matrix("",nrow=m,ncol=length(header))
  for (i in 1:3) {
    pp[,i] <- sapply(x,function(z){z[i]})
  }
  for (i in 4:length(header)) {
    v <- sapply(x,function(z){z[i]})
    v <- gsub("1","0",v)
    pp[,i] <- gsub("2","1",v)
  }
  colnames(pp) <- header

  k <- grep("genoprob",temp[ix])
  id <- temp[ix[k]+1]
  id <- strsplit(id,split=",",fixed=T)[[1]][-(1:3)]
  keep <- 1:length(id)
  n <- length(keep)
  genoprob <- matrix("",nrow=m,ncol=n+1)
  colnames(genoprob) <- c("marker",id[keep])
  for (j in 1:m) {
    i <- ix[k]+1+j
    x <- strsplit(temp[i],split=",",fixed=T)[[1]]
    genoprob[j,] <- x[c(1,keep+3)]
  }
  df <- cbind(pp[,1:3],genoprob[match(pp[,1],genoprob[,1]),-1])
  close(con)
  return(df)
}

#' Converts polyOrigin genotype probabilities format to mappoly homoprob format
#'
#' @param df data.frame resulted from read_polyOrigin function
#' @param f1.codes data.frame with polyOrigin stages codification
#' @param ind character defining the individual to be evalueated
#' @param ploidy integer defining the ploidy
#'
#' @importFrom tidyr pivot_longer
#' @import parallel
#' @export
get_probs_polyorigin <- function(df, f1.codes, ind = NULL, ploidy = 4, n.cores= 1){
  df_prob <- data.frame(df)
  if(!is.null(ind)) df_prob <- df_prob[,c(1:3,which(colnames(df_prob) %in% ind))]
  melt_df <- pivot_longer(df_prob, cols = 4:ncol(df_prob), names_to = "ind", values_to = "geno")

  melt_df_lst <- split(melt_df, melt_df$ind)

  clust <- makeCluster(n.cores)
  clusterExport(clust, c("get_probs_sing"))
  homoprob.lst <- parLapply(clust, melt_df_lst, function(x) get_probs_sing(melt_df = x,
                                                                           f1.codes = f1.codes,
                                                                           ploidy = ploidy))
  stopCluster(clust)

  homoprob <- do.call(rbind, homoprob.lst)

  structure(list(info = list(ploidy = ploidy,
                             n.ind = length(unique(homoprob$individual))) ,
                 homoprob = homoprob), class = "mappoly.homoprob")
}

#' Converts polyOrigin genotype probabilities format to mappoly homoprob format for a single individual
#'
#' @param melt_df long format of the data.frame resulted from read_polyOrigin function
#' @param f1.codes data.frame with polyOrigin stages codification
#' @param ploidy integer defining the ploidy
#'
#' @importFrom tidyr pivot_longer
#'
#' @export
get_probs_sing <- function(melt_df, f1.codes, ploidy){
  df_sep <- strsplit(melt_df$geno, "=>")
  states <- sapply(df_sep, "[[", 1)
  states <- strsplit(states, "[|]")
  states <- lapply(states, function(x) f1.codes$State[match(x, f1.codes$Code)])
  probs <- sapply(df_sep, "[[", 2)
  probs <- strsplit(probs, "[|]")
  probs <- lapply(probs, function(x) as.numeric(x)/sum(as.numeric(x)))
  #probs <- lapply(probs, function(x) as.numeric(x))

  p.probs <- list()
  for(i in 1:(2*ploidy)){
    p.probs[[i]] <- mapply(function(states, probs){
      sum(probs[grep(i, states)])
    }, states, probs)
  }

  # Detect double reduction
  p.probs <- do.call(cbind,p.probs)
  p.probs[,1:ploidy] <- t(apply(p.probs[,1:ploidy], 1, function(x) (2*x)/sum(x)))
  p.probs[,(ploidy+1):(2*ploidy)] <- t(apply(p.probs[,(ploidy+1):(2*ploidy)], 1, function(x) (2*x)/sum(x)))


  colnames(p.probs) <- letters[1:(ploidy*2)]

  df_all <- cbind(melt_df[,-ncol(melt_df)],p.probs)
  homoprob <- pivot_longer(df_all, cols = 5:ncol(df_all), names_to = "homolog", values_to = "probability")
  homoprob <- data.frame(marker = homoprob$marker, homolog = homoprob$homolog,
                         individual = homoprob$ind, probability = homoprob$probability,
                         map.position = as.numeric(homoprob$cM), LG = homoprob$chrom)
  print(paste0("Individual ", unique(homoprob$ind), " done!"))
  return(homoprob)
}
