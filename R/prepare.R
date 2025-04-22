globalVariables(c("ind", "zscore", "chr", "write.csv"))



#' Convert VCF File to Qploidy Data
#'
#' This function converts a VCF file into a format compatible with Qploidy analysis. It extracts genotype and allele depth information and formats it into a data frame.
#'
#' @param vcf_file Path to the VCF file.
#' @param geno Logical. If TRUE, the output columns will include MarkerName, SampleName, geno, and prob. If FALSE, the output will include MarkerName, SampleName, X, Y, R, and ratio.
#' @param geno.pos Logical. If TRUE, the output will include MarkerName, Chromosome, and Position columns.
#' @return A data frame containing the processed VCF data.
#' @export
#' @import vcfR
#' @importFrom tidyr pivot_longer
qploidy_read_vcf <- function(vcf_file, geno = FALSE, geno.pos = FALSE) {

  vcf <- read.vcfR(vcf_file)

  if(!(geno | geno.pos)){
    DP <- extract.gt(vcf, "AD")

    mknames <- pivot_longer(data.frame(mks = rownames(DP), DP), cols = 2:(ncol(DP)+1))

    dp_split <- strsplit(mknames$value, ",")

    # remove markers without two alleles
    dp_split_filt <- lapply(dp_split, function(x){
      if(length(x) != 2) return(c(NA,NA)) else return(x)
    })

    ref <- as.numeric(sapply(dp_split_filt, "[[", 1))
    alt <- as.numeric(sapply(dp_split_filt, "[[", 2))

    data_qploidy <- data.frame(MarkerName = mknames$mks,
                               SampleName = mknames$name,
                               X= ref,
                               Y= alt,
                               R = alt+ref,
                               ratio = alt/(ref+alt))
  } else if(geno){
    GT <- extract.gt(vcf, "GT")
    GT <- pivot_longer(data.frame(mks = rownames(GT), GT), cols = 2:(ncol(GT)+1))
    GT$value <- stringr::str_count(GT$value, "1")
    colnames(GT) <- c("MarkerName","SampleName","geno")

    prob <- extract.gt(vcf, "MPP")
    prob <- pivot_longer(data.frame(mks = rownames(prob), prob), cols = 2:(ncol(prob)+1))
    colnames(prob) <- c("MarkerName","SampleName","prob")

    data_qploidy <- merge.data.frame(GT, prob, by = c("MarkerName", "SampleName"))
    data_qploidy$prob <- as.numeric(data_qploidy$prob)


  } else if(geno.pos){
    data_qploidy <- data.frame("MarkerName" = vcf@fix[,3],
                               "Chromosome" = vcf@fix[,1],
                               "Position" = vcf@fix[,2])


  }

  return(data_qploidy)
}

#' Read Illumina Array Files
#'
#' This function reads Illumina array files and processes them into a format suitable for Qploidy analysis. It adds a suffix to sample IDs if multiple files are provided.
#'
#' @param ... One or more Illumina array filenames.
#' @return A data frame containing the processed Illumina array data.
#' @export
#' @importFrom tidyr pivot_longer
read_illumina_array <- function(...){
  files <- list(...)

  raw_all <- vector()
  for(i in 1:length(files)){
    skipn <- find_header_line(files[[i]], word = "SNP Name")
    raw <- vroom(files[[i]], delim = "\t", skip = skipn-1)
    raw1 <- raw[,c("SNP Name","Sample ID","X","Y")]
    if(length(files) > 1) raw1$`Sample ID` <- paste0(raw1$`Sample ID`, "_",i)
    raw_all <- rbind(raw_all, raw1)
  }

  fitpoly_potato_input <- cbind(raw_all, R = raw_all$X + raw_all$Y, ratio = raw_all$Y/(raw_all$X + raw_all$Y))
  colnames(fitpoly_potato_input)[1:2] <- c("MarkerName", "SampleName")
  return(fitpoly_potato_input)
}


#' Find the Header Line in a File
#'
#' This function scans a file to locate the first line containing a specific keyword, such as 'probeset_id'.
#' It is useful for identifying the starting point of data in files with headers or metadata.
#'
#' @param summary_file The path to the file to be scanned.
#' @param max_lines The maximum number of lines to scan. Default is 6000.
#' @param word The keyword to search for in the first column. Default is "probeset_id".
#' @return The line number where the keyword is found.
#' 
#' @export
find_header_line <- function(summary_file, word = "probeset_id", max_lines = 6000) {
  con <- file(summary_file, open = "r")
  on.exit(close(con))
  
  for (i in 1:max_lines) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) break  # EOF
    if (startsWith(line, word)) return(i)
  }
  
  stop(substitute(word)," not found in the first column within the first ", max_lines, " lines.")
}

#' Convert Axiom Array Summary File to Qploidy Input
#'
#' This function processes an Axiom array summary file and converts it into a format compatible with Qploidy and fitpoly analysis.
#'
#' @param summary_file Path to the Axiom summary file.
#' @param ind_names Optional. A file with two columns: Plate_name (sample IDs in the summary file) and Sample_Name (desired sample names).
#' @param atan Logical. If TRUE, calculates theta using atan2.
#' @return A data frame formatted for Qploidy analysis.
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom tidyr pivot_wider
read_axiom <- function(summary_file, ind_names = NULL, atan = FALSE){

  skip_line <- find_header_line(summary_file)
  raw <- vroom(summary_file, skip = as.numeric(skip_line) -1)

  # Guarantee that every A has a B
  summary_filt <- clean_summary(raw)

  if(!is.null(ind_names)){
    ind.names <- vroom(ind_names)

    new.names <- ind.names$Sample_Name[match(colnames(summary_filt$A_probes), ind.names$Plate_Name)]
    colnames(summary_filt$A_probes)[which(!is.na(new.names))] <- new.names[which(!is.na(new.names))]
    colnames(summary_filt$B_probes)[which(!is.na(new.names))] <- new.names[which(!is.na(new.names))]
  }

  R_theta <- get_R_theta(cleaned_summary = summary_filt, atan)

  fitpoly_input <- summary_to_fitpoly(R_theta$R_all, R_theta$theta_all)

  return(fitpoly_input)
}


#' Clean Axiom Summary File
#'
#' This function removes consecutive A allele probes from an Axiom summary file.
#'
#' @param summary_df A data frame containing A and B probe intensities.
#' @return A list with cleaned A and B probes.
#' @examples NULL
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

#' Get R and Theta Values from Summary File
#'
#' This function calculates R and theta values from a cleaned summary file. It optionally performs standard normalization by plate and markers.
#'
#' @param cleaned_summary A summary object from the clean_summary function.
#' @param atan Logical. If TRUE, calculates theta using atan2.
#' @return A list containing R and theta values.
#' @export
#' @importFrom dplyr mutate
get_R_theta <- function(cleaned_summary, atan = FALSE){
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

  colnames(theta_all)[1] <- colnames(R_all)[1] <- "MarkerName"

  return(list(R_all=R_all, theta_all = theta_all))
}

#' Convert Summary file to fitpoly format
#'
#' @param R_all data.frame with R values
#' @param theta_all data.frame with theta values
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


