globalVariables(c("theta", "R", "geno", "Var1", "array.id",
                  "ids", "sd", "median", "z", "lm", "sd",
                  "snp", "centers_theta", "pnorm"))


##' Standardize allelic ratio data and compute BAF and Z-scores
##'
##' Performs signal standardization of genotype data by aligning `theta` values (allelic ratios or normalized intensities)
##' to expected genotype clusters. Outputs standardized BAF (B-allele frequency) and Z-scores per sample and marker.
##'
##' Reference genotypes are used to estimate cluster centers either from dosage data (e.g., via `fitpoly` or `updog`)
##' or using an `updog` `multidog` object directly. Supports both array-based (intensity) and sequencing-based (count) data.
##'
##' The function applies marker and genotype-level quality filters, uses parallel computing to estimate BAF, and generates
##' a final annotated output suitable for CNV or dosage variation analyses.
##'
##' Filtering steps include:
##'   1. Genotype probability filter: datapoints with probability below `threshold.geno.prob` are set to missing.
##'   2. Marker missingness filter: Markers with fraction of missing datapoints above `threshold.missing.geno` are removed.
##'   3. Cluster number filter: Markers with fewer than `threshold.n.clusters` clusters are removed.
##'   4. Genomic info filter: Markers lacking chromosome or position info are removed.
##'
##' Merging logic:
##'   - The function merges filtered genotype and signal data, estimates cluster centers, computes BAFs in parallel,
##'     and calculates Z-scores. Results are merged into a single output data.frame containing BAF, Z-score, and genotype info.
##'
##' @param data A `data.frame` containing the full dataset with columns:
##'   - MarkerName: Marker identifiers
##'   - SampleName: Sample identifiers
##'   - X: Reference allele intensity or count
##'   - Y: Alternative allele intensity or count
##'   - R: Total signal intensity or read depth (X + Y)
##'   - ratio: Allelic ratio, typically Y / (X + Y)
##'
##' @param genos A `data.frame` containing genotype dosage information for the reference panel, with columns:
##'   - MarkerName: Marker identifiers
##'   - SampleName: Sample identifiers
##'   - geno: Estimated dosage (0, 1, 2, ...)
##'   - prob: Genotype call probability (used for filtering low-confidence datapoints)
##'
##' @param geno.pos A `data.frame` with marker position metadata, with columns:
##'   - MarkerName: Marker identifiers
##'   - Chromosome: Chromosome names
##'   - Position: Base-pair positions on the genome
##'
##' @param threshold.missing.geno Numeric (0–1). Maximum fraction of missing datapoints allowed per marker. Markers with a higher fraction will be removed.
##' @param threshold.geno.prob Numeric (0–1). Minimum genotype call probability threshold. Datapoints with lower probability will be treated as missing.
##' @param ploidy.standardization Integer. The ploidy level of the reference panel used for standardization.
##' @param threshold.n.clusters Integer. Minimum number of expected dosage clusters per marker. For diploid data, this is typically 3 (corresponding to dosages 0, 1, and 2). Cannot exceed `ploidy.standardization + 1`.
##' @param n.cores Integer. Number of cores to use in parallel computations (e.g., for cluster center estimation and BAF generation).
##' @param out_filename Optional. Path to save the final standardized dataset to disk as a CSV file (suitable for Qploidy).
##' @param type Character. Type of data used for clustering: "intensities" (array-based), "counts" (sequencing), or "updog" (set automatically if `multidog_obj` is provided).
##' @param multidog_obj Optional. An object of class `multidog` from the `updog` package, containing model fits and estimated biases. If provided, this will override the `type` parameter and use `updog`'s expected cluster positions.
##' @param parallel.type Character. Parallel backend to use ("FORK" or "PSOCK"). "FORK" is faster but only works on Unix-like systems.
##' @param verbose Logical. If TRUE, prints progress and filtering information to the console.
##' @param rm_outlier Logical. If TRUE, uses Bonferroni-Holm corrected residuals to remove outliers before estimating cluster centers.
##' @param cluster_median Logical. If TRUE, uses the median of theta values to estimate cluster centers. If FALSE, uses the mean.
##' @param filter_R Logical. If TRUE, calculates Z-scores only for markers that passed missing data filter (default: FALSE).
##'
##' @return An object of class `qploidy_standardization` (list) with the following components:
##'   - info: Named vector of standardization parameters
##'   - filters: Named vector summarizing how many markers were removed at each filtering step
##'   - data: A data.frame containing merged BAF, Z-score, and genotype information by marker and sample
##'
##' @importFrom dplyr group_by summarize filter mutate inner_join full_join
##' @importFrom tidyr pivot_wider pivot_longer
##' @importFrom vroom vroom_write
##' @importFrom parallel makeCluster stopCluster clusterExport parLapply clusterEvalQ
##' @importFrom stats sd
##'
##' @examples
##' # Example usage:
##' # data <- ... # see vignette for example data
##' # genos <- ...
##' # geno.pos <- ...
##' # result <- standardize(
##' #   data = data,
##' #   genos = genos,
##' #   geno.pos = geno.pos,
##' #   ploidy.standardization = 2,
##' #   threshold.n.clusters = 3,
##' #   n.cores = 2
##' # )
##'
##' @export
standardize <- function(data = NULL,
                        genos = NULL,
                        geno.pos = NULL,
                        threshold.missing.geno=0.9,
                        threshold.geno.prob=0.8,
                        ploidy.standardization = NULL,
                        threshold.n.clusters = NULL,
                        n.cores =1,
                        out_filename = NULL,
                        type = "intensities",
                        multidog_obj = NULL,
                        parallel.type = "PSOCK",
                        verbose = TRUE,
                        rm_outlier = TRUE,
                        cluster_median = TRUE,
                        filter_R = FALSE){

  if(is.null(data) || is.null(genos) || is.null(geno.pos) ||
     is.null(ploidy.standardization) || is.null(threshold.n.clusters)) {
    stop("Not all required inputs were defined.")
  }

  # Scalar type checks
  if (!is.data.frame(data))     stop("'data' must be a data.frame.")
  if (!is.data.frame(genos))    stop("'genos' must be a data.frame.")
  if (!is.data.frame(geno.pos)) stop("'geno.pos' must be a data.frame.")
  if (nrow(data) == 0)     stop("'data' has no rows.")
  if (nrow(genos) == 0)    stop("'genos' has no rows.")
  if (nrow(geno.pos) == 0) stop("'geno.pos' has no rows.")

  if (!is.numeric(threshold.missing.geno) || length(threshold.missing.geno) != 1 ||
      threshold.missing.geno < 0 || threshold.missing.geno > 1)
    stop("'threshold.missing.geno' must be a single numeric value in [0, 1].")

  if (!is.numeric(threshold.geno.prob) || length(threshold.geno.prob) != 1 ||
      threshold.geno.prob < 0 || threshold.geno.prob > 1)
    stop("'threshold.geno.prob' must be a single numeric value in [0, 1].")

  if (!is.numeric(ploidy.standardization) || length(ploidy.standardization) != 1 ||
      ploidy.standardization < 1 || ploidy.standardization != as.integer(ploidy.standardization))
    stop("'ploidy.standardization' must be a single positive integer.")

  if (!is.null(threshold.n.clusters) &&
      (!is.numeric(threshold.n.clusters) || length(threshold.n.clusters) != 1 ||
       threshold.n.clusters < 1 || threshold.n.clusters != as.integer(threshold.n.clusters)))
    stop("'threshold.n.clusters' must be a single positive integer.")

  if (!is.null(threshold.n.clusters) && threshold.n.clusters > ploidy.standardization + 1)
    stop(sprintf("'threshold.n.clusters' cannot be higher than ploidy.standardization + 1 (%d).",
                 ploidy.standardization + 1L))

  if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1 || n.cores != as.integer(n.cores))
    stop("'n.cores' must be a single positive integer.")

  if (!is.null(out_filename) && (!is.character(out_filename) || length(out_filename) != 1))
    stop("'out_filename' must be a single character string or NULL.")

  allowed_types <- c("intensities", "counts", "updog")
  if (!is.character(type) || length(type) != 1 || !type %in% allowed_types)
    stop(sprintf("'type' must be one of: %s.", paste(allowed_types, collapse = ", ")))

  allowed_parallel <- c("FORK", "PSOCK")
  if (!is.character(parallel.type) || length(parallel.type) != 1 || !parallel.type %in% allowed_parallel)
    stop(sprintf("'parallel.type' must be one of: %s.", paste(allowed_parallel, collapse = ", ")))

  if (!is.logical(verbose)       || length(verbose) != 1)       stop("'verbose' must be a single logical value.")
  if (!is.logical(rm_outlier)    || length(rm_outlier) != 1)    stop("'rm_outlier' must be a single logical value.")
  if (!is.logical(cluster_median)|| length(cluster_median) != 1) stop("'cluster_median' must be a single logical value.")
  if (!is.logical(filter_R)      || length(filter_R) != 1)      stop("'filter_R' must be a single logical value.")

  if (!is.null(multidog_obj) && !inherits(multidog_obj, "multidog"))
    stop("'multidog_obj' must be an object of class 'multidog' (from the updog package) or NULL.")

  # Check column names and order for data
  required_data_cols <- c("MarkerName", "SampleName", "X", "Y", "R", "ratio")
  if (!identical(colnames(data), required_data_cols)) {
    stop(sprintf(
      "Column names of the provided data object must be exactly: %s (in this order).",
      paste(required_data_cols, collapse = ", ")
    ))
  }
  # Check column names and order for genos
  allowed_genos_cols <- list(
    c("MarkerName", "SampleName", "geno"),
    c("MarkerName", "SampleName", "geno", "prob")
  )
  genos_colnames <- colnames(genos)
  if (!any(sapply(allowed_genos_cols, function(cols) identical(genos_colnames, cols)))) {
    stop(sprintf(
      "Column names of the provided genos object must be exactly: %s (in this order), or: %s (in this order).",
      paste(allowed_genos_cols[[1]], collapse = ", "),
      paste(allowed_genos_cols[[2]], collapse = ", ")
    ))
  }
  has_prob_col <- identical(genos_colnames, allowed_genos_cols[[2]])

  # Check column names and order for geno.pos
  required_genopos_cols <- c("MarkerName", "Chromosome", "Position")
  if (!identical(colnames(geno.pos), required_genopos_cols)) {
    stop(sprintf(
      "Column names of the provided geno.pos object must be exactly: %s (in this order).",
      paste(required_genopos_cols, collapse = ", ")
    ))
  }

  dose <- max(genos$geno, na.rm = T)
  if(dose != ploidy.standardization) stop("Ploidy of the provided reference samples do not match with the one defined in the ploidy.standardization parameter.")

  vmsg("Generating standardized BAFs", verbose = verbose, level = 0, type = ">>")
  if(is.null(threshold.n.clusters)) threshold.n.clusters <- ploidy.standardization + 1

  ## Filter by prob (only if prob column is present)
  if (has_prob_col) {
    idx <- genos$prob < threshold.geno.prob
    prob.rm <- table(idx)
    if(!is.na(prob.rm["TRUE"])) prob.rm <- round(as.numeric(prob.rm["TRUE"]/sum(prob.rm)*100),2) else prob.rm <- 0
    idx <- which(idx)
    if(length(idx) > 0) genos$geno[idx] <- NA
    vmsg("Percentage of datapoints turned into missing data because of low genotype probability: %s", verbose = verbose, level = 2, type = ">>", prob.rm)
  } else {
    prob.rm <- 0
    vmsg("No 'prob' column in genos: skipping genotype probability filtering.", verbose = verbose, level = 1, type = ">>")
  }

  ## Filter by missing data
  n.na <- genos %>% group_by(MarkerName) %>% summarize(n.na = (sum(is.na(geno))/length(geno)))
  rm.mks <- n.na$MarkerName[which(n.na$n.na > threshold.missing.geno)]
  mis.rm <- length(rm.mks)

  vmsg("Markers removed because of excess of missing data: %s", verbose = verbose, level = 2, type = ">>", mis.rm)

  if(length(rm.mks) > 0){
    genos_filt <- genos[which(!(genos$MarkerName %in% rm.mks)),]
    data_filt <- data[which(!(data$MarkerName %in% rm.mks)),]
  } else {
    genos_filt <- genos
    data_filt <- data
  }

  # Organize data
  data_standardization <- inner_join(data_filt[,-c(3,4,5)], genos_filt[,-4], by = c("MarkerName", "SampleName"))

  if(dim(data_standardization)[1] == 0) stop("Individuals in `data` and `genos` don't have same ID.")

  colnames(data_standardization)[3] <- "theta"

  rm.na <- which(is.na(data_standardization$geno))
  if(length(rm.na) > 0){
    data_standardization$theta[rm.na] <- NA
  }

  if(is.null(multidog_obj)){
    lst_standardization <- split(data_standardization, data_standardization$MarkerName)

    vmsg("Going to parallel mode", verbose = verbose, level = 1, type = ">>")
    clust <- makeCluster(n.cores, type = parallel.type)
    clusterExport(clust, c("get_centers", "rm_outlier"), envir = .GlobalEnv)
    clusterEvalQ(clust, { library(magrittr); library(dplyr)}) # load required packages and Qploidy on workers
    clusters <- parLapply(clust, lst_standardization, get_centers,
                          ploidy= ploidy.standardization,
                          n.clusters.thr = threshold.n.clusters,
                          type = type,
                          rm_outlier = rm_outlier,
                          cluster_median = cluster_median)

    stopCluster(clust)

    vmsg("Back to single core usage", verbose = verbose, level = 1, type = ">>")

    gc(verbose = FALSE)

  } else { # centers defined using updog bias
    clusters <- updog_centers(multidog_obj, threshold.n.clusters = threshold.n.clusters, rm.mks = rm.mks)
  }

  # Filter by number of clusters
  rm.mks <- sapply(clusters, function(x) x$rm != 0)
  clusters.rm <- sum(rm.mks)

  vmsg("Markers removed because of smaller number of clusters than set threshold: %s", verbose = verbose, level = 2, type = ">>", clusters.rm)

  if(length(which(rm.mks)) > 0)  clusters_filt <- clusters[-which(rm.mks)] else clusters_filt <- clusters

  if(length(clusters_filt) == 0) stop("All markers were filtered, adapt filtering parameters.")

  keep.mks <- sapply(clusters_filt, function(x) x$MarkerName)

  # Getting BAF for complete dataset
  theta_filt <- pivot_wider(data[which(data$MarkerName %in% keep.mks),-c(3:5)], names_from = "SampleName", values_from = "ratio")
  centers_filt <- t(sapply(clusters_filt, function(x) x$centers_theta))
  centers_filt <- data.frame(mks=keep.mks, centers_filt)
  centers_filt <- centers_filt[match(theta_filt$MarkerName, centers_filt$mks),]

  par <- rep(seq_len(n.cores),
             each = round((nrow(theta_filt) / n.cores) + 1, 0))[seq_len(nrow(theta_filt))]

  par_theta <- split.data.frame(theta_filt, par)
  par_clusters_filt <- split(centers_filt, par)

  if(length(par_theta) < n.cores) n.cores <- length(par_theta)
  par_all <- list()
  for(i in 1:n.cores){
    par_all[[i]] <- list()
    par_all[[i]][[1]] <- as.data.frame(par_theta[[i]])
    par_all[[i]][[2]] <- as.data.frame(par_clusters_filt[[i]])
  }

  # Get BAF
  vmsg("Going to parallel mode", verbose = verbose, level = 1, type = ">>")

  clust <- makeCluster(n.cores, type = parallel.type)
  clusterExport(clust, c("get_baf_par", "get_baf"))
  bafs <- parLapply(clust, par_all, get_baf_par, ploidy = ploidy.standardization)

  stopCluster(clust)

  vmsg("Back to single core usage", verbose = verbose, level = 1, type = ">>")

  gc()

  bafs_lt <- unlist(bafs, recursive = F)
  bafs_m <- do.call(rbind, bafs_lt)
  rownames(bafs_m) <- theta_filt$MarkerName
  colnames(bafs_m) <- colnames(theta_filt)[-1]
  bafs_df <- as.data.frame(bafs_m)
  bafs_df <- cbind(mks = theta_filt$MarkerName, bafs_df)

  # Add chr and pos info
  chr <- geno.pos$Chromosome[match(bafs_df$mks,geno.pos$MarkerName)]
  pos <- geno.pos$Position[match(bafs_df$mks,geno.pos$MarkerName)]

  bafs_join <- cbind(MarkerName=bafs_df$mks, Chr = chr, Position = as.numeric(pos), bafs_df[,-1])

  no.geno.info <- length(which(is.na(bafs_join$Chr)))
  if(no.geno.info > 0)
    bafs_join <- bafs_join[-which(is.na(bafs_join$Chr)),]

  baf_melt <- pivot_longer(bafs_join, cols = 4:ncol(bafs_join), names_to = "SampleName", values_to = "baf")
  vmsg("BAFs ready!", verbose = verbose, level = 1, type = ">>")

  # Z score
  vmsg("Generating z scores", verbose = verbose, level = 0, type = ">>")
  if(filter_R)
    zscore <- get_zscore(data_filt, geno.pos) # z score is calculated only for markers that passed missing data filter
  else zscore <- get_zscore(data, geno.pos)

  vmsg("Z scores ready!", verbose = verbose, level = 1, type = ">>")

  vmsg("Preparing outputs", verbose = verbose, level = 0, type = ">>")

  vmsg("Merging results into qploidy_standardization object", verbose = verbose, level = 1, type = ">>")
  qploidy_data <- full_join(data, data_standardization[,-3], c("MarkerName", "SampleName"))
  qploidy_data <- full_join(qploidy_data,baf_melt[,-c(2,3)], c("MarkerName", "SampleName"))
  qploidy_data <- full_join(qploidy_data, zscore[,-c(2,3)], c("MarkerName", "SampleName"))

  # Fill Chr and Position columns in qploidy_data using geno.pos
  marker_idx <- match(qploidy_data$MarkerName, geno.pos$MarkerName)
  qploidy_data$Chr <- geno.pos$Chromosome[marker_idx]
  qploidy_data$Position <- as.numeric(geno.pos$Position[marker_idx])

  result <- structure(list(info = c(threshold.missing.geno = threshold.missing.geno,
                                    threshold.geno.prob = threshold.geno.prob,
                                    ploidy.standardization = ploidy.standardization,
                                    threshold.n.clusters = threshold.n.clusters,
                                    out_filename = out_filename,
                                    type = if(!is.null(multidog_obj)) "updog" else type),
                           filters = c(n.markers.start = length(unique(data$MarkerName)),
                                       geno.prob.rm = prob.rm,
                                       miss.rm = mis.rm,
                                       clusters.rm= clusters.rm,
                                       no.geno.info.rm = no.geno.info,
                                       n.markers.end = length(unique(baf_melt$MarkerName))),
                           data = qploidy_data), class = "qploidy_standardization")

  if(!is.null(out_filename)) {
    vmsg("Writting QploidyApp input file: %s", verbose = verbose, level = 1, type = ">>", out_filename)
    write_qploidy_standardization(result, out_filename)
  }

  vmsg("Done!", verbose = verbose, level = 1, type = ">>")
  return(result)
}


#' Print method for object of class 'qploidy_standardization'
#'
#' @param x object of class 'qploidy_standardization'
#' @param ... print parameters
#'
#' @return printed information about Qploidy standardization process
#'
#' @method print qploidy_standardization
#'
#' @export
print.qploidy_standardization <- function(x, ...){

  info <- data.frame(c1 = c("Standardization type:", "Ploidy:",
                            "Min # of heterozygous classes (clusters) present:",
                            "Max proportion of missing genotype by marker:",
                            "Min genotype probability:"),
                     c2 = c(x$info["type"], x$info["ploidy.standardization"],
                            x$info["threshold.n.clusters"],
                            as.numeric(x$info["threshold.missing.geno"]),
                            x$info["threshold.geno.prob"]))

  format.df <- data.frame(c1 = c("# markers at raw data:",
                                 "% datapoints filtered by low probability:",
                                 "# markers filtered by missing data:",
                                 "# markers filtered by min number of clusters:",
                                 "# markers filtered by lack of genomic information:",
                                 "# markers remaining with estimated BAF:"),
                          c2 = c(x$filters["n.markers.start"],
                                 "-",
                                 x$filters["miss.rm"],
                                 x$filters["clusters.rm"],
                                 x$filters["no.geno.info.rm"],
                                 x$filters["n.markers.end"]),
                          c3 = c("(100%)",
                                 paste0("(",x$filters["geno.prob.rm"], " %)"),
                                 paste0("(",round(x$filters["miss.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["clusters.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["no.geno.info.rm"]/x$filters["n.markers.start"]*100,2)," %)"),
                                 paste0("(", round(x$filters["n.markers.end"]/x$filters["n.markers.start"]*100,2)," %)")))

  colnames(info) <- rownames(info) <- colnames(format.df) <- rownames(format.df) <- NULL

  cat("This is an object of class 'qploidy_standardization'\n")
  cat("--------------------------------------------------------------------\n")
  cat("Parameters\n")
  print(format(info, justify="left", digit = 2))
  cat("--------------------------------------------------------------------\n")
  cat("Filters\n")
  print(format(format.df, justify="left", digit = 2))
}

#' Read a Qploidy Standardized Dataset (.tsv or .tsv.gz) with format checks
#'
#' Expected layout (tab-separated):
#'   1) "info" header line
#'   2) one row of info values
#'   3) "filters" header line
#'   4) one row of numeric filter values
#'   5) "data" header line
#'   6+) data rows (must include at least SampleName and Chr)
#'
#' @param qploidy_standardization_file Path to the TSV/TSV.GZ file
#' @return An object of class "qploidy_standardization" with fields:
#'   - info:    named character vector (single row)
#'   - filters: named numeric vector (single row, typically >= 6 values)
#'   - data:    data.frame/tibble with at least SampleName and Chr
#' @export
read_qploidy_standardization <- function(qploidy_standardization_file) {
  if (is.null(qploidy_standardization_file) || !nzchar(qploidy_standardization_file)) {
    stop("Path is empty or NULL.", call. = FALSE)
  }
  if (!file.exists(qploidy_standardization_file)) {
    stop("File does not exist: ", qploidy_standardization_file, call. = FALSE)
  }

  # Soft extension check (warn if unexpected)
  if (!grepl("\\.(tsv|txt)(\\.gz)?$", tolower(qploidy_standardization_file))) {
    warning("File extension is not .tsv[.gz]; expecting a tab-separated file. Attempting to read anyway.")
  }

  # Read sections (force tab delimiter; suppress noisy messages)
  suppressMessages({
    info <- vroom::vroom(
      qploidy_standardization_file,
      delim = "\t", n_max = 1, progress = FALSE, show_col_types = FALSE
    )
  })
  # If we only got 1 column, it's likely not tab-separated
  if (ncol(info) <= 1) {
    stop("File must be tab-separated (.tsv or .tsv.gz). Detected a single-column first section.", call. = FALSE)
  }
  if (nrow(info) != 1L) {
    stop("File format error: expected the 'info' section to have exactly one data row (line 2).", call. = FALSE)
  }
  info_v <- as.character(info[1, , drop = TRUE])
  names(info_v) <- names(info)

  suppressMessages({
    filters <- vroom::vroom(
      qploidy_standardization_file,
      delim = "\t", skip = 2, n_max = 1, progress = FALSE, show_col_types = FALSE
    )
  })
  if (ncol(filters) < 1L || nrow(filters) != 1L) {
    stop("File format error: expected a 'filters' header (line 3) followed by exactly one row of numeric values (line 4).", call. = FALSE)
  }
  filters_v <- suppressWarnings(as.numeric(filters[1, , drop = TRUE]))
  names(filters_v) <- names(filters)
  if (any(is.na(filters_v))) {
    stop("File format error: non-numeric values found in the 'filters' row (line 4).", call. = FALSE)
  }
  if (length(filters_v) < 6L) {
    warning("Filters row has fewer than 6 values; downstream UI may not show all counters.")
  }

  suppressMessages({
    data <- vroom::vroom(
      qploidy_standardization_file,
      delim = "\t", skip = 4, progress = FALSE, show_col_types = FALSE
    )
  })
  if (!nrow(data)) {
    stop("File format error: data section is empty (no rows after line 5).", call. = FALSE)
  }

  # Normalize common name variants for downstream compatibility
  nm <- names(data)
  if (!("Chr" %in% nm) && ("Chromosome" %in% nm)) {
    names(data)[nm == "Chromosome"] <- "Chr"
  }
  nm <- names(data)
  if (!("Pos" %in% nm) && ("Position" %in% nm)) {
    names(data)[nm == "Position"] <- "Position"
  }

  data$Position <- as.numeric(data$Position)

  # Minimal required columns used by the module/UI
  required_cols <- c("SampleName", "Chr")
  missing <- setdiff(required_cols, names(data))
  if (length(missing)) {
    stop(
      sprintf("File format error: data section must contain column(s): %s.",
              paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  structure(
    list(info = info_v, filters = filters_v, data = data),
    class = "qploidy_standardization"
  )
}

#' Write Qploidy Standardization Object to File
#'
#' This function writes a `qploidy_standardization` object to a specified file.
#' The output file includes metadata, filtering information, and the standardized dataset.
#'
#' @param qploidy_standardization_object An object of class `qploidy_standardization` to be written to file.
#' @param out_filename A string specifying the path to the output file where the data will be saved. Specify in the file name the desired extension (e.g., `my_output.csv`).
#' @return None. The function writes the data to the specified file.
#' @import vroom
#' @export
write_qploidy_standardization <- function(qploidy_standardization_object, out_filename){
  if(!inherits(qploidy_standardization_object, "qploidy_standardization")) stop("The provided object is not of class 'qploidy_standardization'.")

  info = data.frame(t(qploidy_standardization_object$info))
  filters = data.frame(t(qploidy_standardization_object$filters))
  vroom_write(info, file = out_filename, col_names = T)
  vroom_write(filters, file = out_filename, append = TRUE, col_names = T)
  vroom_write(qploidy_standardization_object$data, file = out_filename, append = TRUE, col_names = T)
}


##' Re-standardize a Qploidy HMM object
##'
##' Re-runs the Qploidy standardization pipeline on a HMM object (from `hmm_estimate_CN_multi`),
##' generating a new `qploidy_standardization` object for downstream analysis or Shiny app input.
##' This is useful for re-standardizing after model selection or parameter tuning.
##'
##' The function extracts marker-level data and genotype dosages for a selected model, applies quality filters,
##' estimates cluster centers, computes BAF (B-allele frequency) and Z-scores, and merges results into a standardized data frame.
##'
##' @param data A `data.frame` containing the full dataset with columns:
##'   - MarkerName: Marker identifiers
##'   - SampleName: Sample identifiers
##'   - X: Reference allele intensity or count
##'   - Y: Alternative allele intensity or count
##'   - R: Total signal intensity or read depth (X + Y)
##'   - ratio: Allelic ratio, typically Y / (X + Y)
##'
##' @param geno.pos A `data.frame` with marker genomic positions, with columns:
##'   - MarkerName: Marker identifiers
##'   - Chromosome: Chromosome identifier
##'   - Position: Genomic position (base-pair coordinate)
##'
##' @param hmm_CN_multi An object of class `hmm_CN` as returned by `hmm_estimate_CN_multi`.
##' @param selected_model Character. The name of the model to use from `hmm_CN_multi$params_samples` (required).
##' @param ploidy.standardization Integer. Ploidy level to use for standardization. If NULL, uses the mode of CN_call in dosages.
##' @param threshold.n.clusters Integer. Minimum number of expected dosage clusters per marker. Defaults to ploidy + 1.
##' @param n.cores Integer. Number of cores for parallel computation. Default is 1.
##' @param threshold.geno.prob Numeric (0–1). Minimum genotype call probability. Default is 0.5.
##' @param threshold.missing.geno Numeric (0–1). Maximum fraction of missing datapoints per marker. Default is 0.90.
##' @param out_filename Optional. Path to save the standardized dataset (CSV/TSV).
##' @param type Character. Data type for clustering: "intensities" (default), "counts", or "updog".
##' @param multidog_obj Optional. updog multidog object for cluster center estimation.
##' @param parallel.type Character. Parallel backend: "FORK" or "PSOCK". Default is "PSOCK".
##' @param rm_outlier Logical. Remove outliers before estimating cluster centers. Default is TRUE.
##' @param cluster_median Logical. Use median (TRUE) or mean (FALSE) for cluster centers. Default is TRUE.
##' @param verbose Logical. Print progress messages. Default is TRUE.
##'
##' @return An object of class `qploidy_standardization` (list) with elements:
##'   - info: Named vector of standardization parameters
##'   - filters: Named vector summarizing filtering steps
##'   - data: Data frame with merged BAF, Z-score, and genotype info by marker and sample
##'
##' @details
##' This function is intended for re-standardizing a Qploidy HMM object after model selection or parameter changes.
##' It extracts the relevant marker-level and dosage data, applies quality filters, estimates cluster centers,
##' computes BAF and Z-scores, and merges all results into a standardized data frame suitable for downstream analysis or Shiny app input.
##'
##' Filtering steps include:
##'   1. Genotype probability filter: datapoints with probability below `threshold.geno.prob` are set to missing.
##'   2. Marker missingness filter: Markers with fraction of missing datapoints above `threshold.missing.geno` are removed.
##'   3. Cluster number filter: Markers with fewer than `threshold.n.clusters` clusters are removed.
##'   4. Genomic info filter: Markers lacking chromosome or position info are removed.
##'
##' @examples
##' # Example usage:
##' # hmm_CN_multi <- hmm_estimate_CN_multi(...)
##' # std <- re_standardize(
##' #   data = my_data,
##' #   geno.pos = my_geno_pos,
##' #   hmm_CN_multi = hmm_CN_multi,
##' #   selected_model = "model1",
##' #   ploidy.standardization = 4
##' # )
##'
##' @export
re_standardize <- function(data = NULL,
                           geno.pos = NULL,
                           hmm_CN_multi,
                           selected_model = NULL,
                           ploidy.standardization = NULL,
                           threshold.n.clusters = NULL,
                           n.cores = 1,
                           threshold.geno.prob = 0.5,
                           threshold.missing.geno = 0.90,
                           out_filename = NULL,
                           type = "intensities",
                           multidog_obj = NULL,
                           parallel.type = "PSOCK",
                           rm_outlier = TRUE,
                           cluster_median = TRUE,
                           verbose = TRUE) {

  # Check input object
  # check if hmm_CN_multi is hmm_CN class
  if (!inherits(hmm_CN_multi, "hmm_CN")) stop("Input object must be of class 'hmm_CN'")
  # check if hmm_CN_multi has params_samples element
  if (!any(names(hmm_CN_multi) == "params_samples")) stop("Input object must come from hmm_estimate_CN_multi function")

  # make selected_model required
  if (is.null(selected_model)) stop("selected_model must be provided")

  # Call dosages for the selected model
  dosages <- call_hmm_dosages(hmm_CN = hmm_CN_multi, selected_model)

  # Set ploidy.standardization if not provided
  if(is.null(ploidy.standardization)) {
    ploidy.standardization <- mode( hmm_CN_multi$by_marker$CN_call)
    if (verbose) cat("ploidy.standardization not provided, using:", ploidy.standardization, "\n")
  }

  dosages[which(dosages$CN_call != ploidy.standardization), "dosage"] <- NA

  genos <- dosages[, c("MarkerName", "SampleName", "dosage", "post_max_dosage")]
  colnames(genos)[3:4] <- c("geno", "prob")
  genos <- genos[order(genos$MarkerName, genos$SampleName),]

  if (is.null(threshold.n.clusters)) {
    threshold.n.clusters <- ploidy.standardization + 1
    if (verbose) cat("threshold.n.clusters not provided, using:", threshold.n.clusters, "\n")
  }

  re_qploidy_standardization <- standardize(
    data = data,
    genos = genos,
    geno.pos = geno.pos,
    ploidy.standardization = ploidy.standardization,
    threshold.n.clusters = threshold.n.clusters,
    n.cores = n.cores,
    threshold.geno.prob = threshold.geno.prob,
    threshold.missing.geno = threshold.missing.geno,
    out_filename = out_filename,
    type = type,
    multidog_obj = multidog_obj,
    parallel.type = parallel.type,
    rm_outlier = rm_outlier,
    cluster_median = cluster_median,
    verbose = verbose
  )

  return(re_qploidy_standardization)
}
