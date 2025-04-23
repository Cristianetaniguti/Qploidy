#' Simulate a VCF file with GT, DP, and AD format fields for 2 chromosomes
#'
#' @param file_path The path where the simulated VCF file will be saved.
#' @param seed The seed for random number generation to ensure reproducibility.
#' @param n_tetraploid Number of tetraploid samples. Default is 35.
#' @param n_diploid Number of diploid samples. Default is 5.
#' @param n_triploid Number of triploid samples. Default is 10.
#' @return None. The function writes the simulated VCF content to the specified file.
#'
#' @export
#' @importFrom stats rbinom
simulate_vcf <- function(
    file_path, seed,
    n_tetraploid = 35,
    n_diploid = 5,
    n_triploid = 10,
    n_markers = 100) {

  set.seed(seed) # For reproducibility

  # Define chromosomes and markers
  chromosomes <- c("chr1", "chr2")
  chr_lengths <- c(50e6, 50e6) # 50Mb each
  n_samples <- sum(n_tetraploid, n_diploid, n_triploid)

  # Define samples and ploidy
  samples <- c(
    paste0("Tetraploid", seq_len(n_tetraploid)),
    paste0("Diploid", seq_len(n_diploid)),
    paste0("Triploid", seq_len(n_triploid))
  )
  ploidies <- c(rep(4, n_tetraploid), rep(2, n_diploid), rep(3, n_triploid))

  # Generate markers for each chromosome
  vcf_content <- "##fileformat=VCFv4.2\n##source=SimulatedVCF\n##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depths\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
  vcf_content <- paste0(vcf_content, paste(samples, collapse = "\t"))

  for (chr in chromosomes) {
    positions <- sort(sample(1:chr_lengths[1], n_markers))
    for (i in seq_along(positions)) {
      ref <- sample(c("A", "T", "G", "C"), 1)
      alt <- sample(setdiff(c("A", "T", "G", "C"), ref), 1)
      qual <- sample(50:100, 1)
      filter <- "PASS"
      info <- paste0("DP=", sample(100:200, 1))

      # Generate genotype data for each sample
      genotypes <- sapply(ploidies, function(ploidy) {
        gt <- paste0(sort(sample(0:1, ploidy, replace = TRUE)), collapse = "/")
        ref_count <- sum(rbinom(ploidy, size = n_samples, prob = stringr::str_count(gt, "0") / ploidy))
        alt_count <- sum(rbinom(ploidy, size = n_samples, prob = stringr::str_count(gt, "1") / ploidy))
        dp <- ref_count + alt_count
        ad <- paste(ref_count, alt_count, sep = ",")
        paste(gt, dp, ad, sep = ":")
      })

      vcf_content <- paste0(
        vcf_content, "\n", chr, "\t", positions[i], "\t", chr, "_mk", i, "\t", ref, "\t", alt, "\t",
        qual, "\t", filter, "\t", info, "\tGT:DP:AD\t",
        paste(genotypes, collapse = "\t")
      )
    }
  }

  writeLines(vcf_content, file_path)
}

#' Simulate an Axiom array summary file
#'
#' This function generates a simulated Axiom array summary file with probe IDs ending in `-A` or `-B` and sample intensities. The intensities are simulated based on the genotype of the sample: homozygous for A, homozygous for B, or heterozygous.
#'
#' @param file_path The path where the simulated summary file will be saved.
#' @param n_probes Number of probes to simulate. Default is 100.
#' @param n_samples Number of samples to simulate. Default is 10.
#' @param seed The seed for random number generation to ensure reproducibility.
#' @return None. The function writes the simulated summary content to the specified file.
#'
#' @export
#' @importFrom utils write.table
simulate_axiom_summary <- function(file_path, n_probes = 100, n_samples = 10, seed) {
  set.seed(seed) # For reproducibility

  # Generate probe IDs
  probe_ids <- paste0("AX-", sample(10000000:99999999, n_probes), c("-A", "-B"))

  # Generate sample names
  sample_names <- paste0("Sample", 1:n_samples)

  # Initialize the data frame
  summary_data <- data.frame(probeset_id = probe_ids)

  # Simulate intensities for each sample
  for (sample in sample_names) {
    intensities <- sapply(probe_ids, function(probe) {
      if (grepl("-A$", probe)) {
        # Intensity for A probe
        if (runif(1) < 0.5) {
          # Homozygous for A
          rnorm(1, mean = 100, sd = 10)
        } else if (runif(1) < 0.5) {
          # Heterozygous
          rnorm(1, mean = 50, sd = 10)
        } else {
          # Homozygous for B
          rnorm(1, mean = 10, sd = 5)
        }
      } else {
        # Intensity for B probe
        if (runif(1) < 0.5) {
          # Homozygous for A
          rnorm(1, mean = 10, sd = 5)
        } else if (runif(1) < 0.5) {
          # Heterozygous
          rnorm(1, mean = 50, sd = 10)
        } else {
          # Homozygous for B
          rnorm(1, mean = 100, sd = 10)
        }
      }
    })
    summary_data[[sample]] <- intensities
  }

  # Write to file
  write.table(summary_data, file = file_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

#' Simulate an Illumina File
#'
#' This function generates a simulated Illumina file with SNP data for a specified number of SNPs and samples.
#' The file includes a header section and a data section with fields such as SNP Name, Sample ID, GC Score, Theta, X, Y, X Raw, Y Raw, and Log R Ratio.
#'
#' @param filepath The path where the simulated Illumina file will be saved. Default is "simulated_summary.txt".
#' @param num_snps The number of SNPs to simulate. Default is 10.
#' @param num_samples The number of samples to simulate. Default is 1.
#' @param sample_id_prefix The prefix for sample IDs. Default is "SAMP".
#' @param mk_id The prefix for marker IDs. Default is "MK-".
#' @param seed The seed for random number generation to ensure reproducibility. Default is 123.
#' @return None. The function writes the simulated Illumina file to the specified path.
#' @details The simulated data includes random values for GC Score, Theta, X, Y, X Raw, Y Raw, and Log R Ratio. The header section provides metadata about the file, including the number of SNPs and samples.
#'
#' @export
simulate_illumina_file <- function(
    filepath,
    num_snps = 10,
    num_samples = 1,
    sample_id_prefix = "SAMP",
    mk_id = "MK-",
    seed = 123) {
  set.seed(seed)

  # --- Create header ---
  header <- c(
    "[Header]",
    "GSGT Version\t2.0.5",
    paste("Processing Date", format(Sys.time(), "%m/%d/%Y %I:%M %p"), sep = "\t"),
    "Content\tSpecies.bpm",
    paste("Num SNPs", num_snps, sep = "\t"),
    paste("Total SNPs", num_snps, sep = "\t"),
    paste("Num Samples", num_samples, sep = "\t"),
    paste("Total Samples", num_samples, sep = "\t"),
    "[Data]"
  )

  # --- Simulate data ---
  snp_names <- paste0(mk_id, seq_len(num_snps))
  sample_ids <- paste0(sample_id_prefix, sprintf("%05d", seq_len(num_samples)))

  data_list <- lapply(sample_ids, function(sample_id) {
    data.frame(
      `SNP Name` = snp_names,
      `Sample ID` = sample_id,
      `GC Score` = round(runif(num_snps, 0, 1), 4),
      `Theta` = round(runif(num_snps, 0, 1), 3),
      `X` = round(runif(num_snps, 0.001, 0.05), 3),
      `Y` = round(runif(num_snps, 0.001, 0.05), 3),
      `X Raw` = sample(300:700, num_snps, replace = TRUE),
      `Y Raw` = sample(10:200, num_snps, replace = TRUE),
      `Log R Ratio` = round(rnorm(num_snps, 0, 1), 7)
    )
  })

  data_df <- do.call(rbind, data_list)
  rownames(data_df) <- NULL
  colnames(data_df) <- c(
    "SNP Name", "Sample ID", "GC Score", "Theta", "X", "Y",
    "X Raw", "Y Raw", "Log R Ratio"
  )

  # --- Write to file ---
  writeLines(header, filepath)
  write.table(
    data_df,
    filepath,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  message("File written to: ", normalizePath(filepath))
}

#' Simulate Genotyping Data with Flexible Ploidy
#'
#' Generates synthetic genotyping and signal intensity data for a given ploidy level.
#' Returns a structured list containing input data suitable for standardization analysis.
#'
#' @param n_markers Integer. Number of markers to simulate (default: 10).
#' @param n_samples Integer. Number of individuals/samples to simulate (default: 5).
#' @param ploidy Integer. Ploidy level of the organism (e.g., 2 for diploid, 4 for tetraploid).
#' @param seed Integer. Random seed for reproducibility (default: 2025).
#'
#' @return A named list with:
#' \describe{
#'   \item{sample_data}{Allelic signal intensities (X, Y, R, ratio).}
#'   \item{geno_data}{Genotype dosage and probability data.}
#'   \item{geno_pos}{Genomic coordinates for each marker.}
#'   \item{standardization_input}{Merged input data with theta and genotype.}
#' }
#'
#' @importFrom stats cor rnorm runif
#'
#' @export
simulate_standardization_input <- function(n_markers = 10,
                                           n_samples = 5,
                                           ploidy = 2,
                                           seed = 2025) {
  set.seed(seed)

  marker_ids <- paste0("m", 1:n_markers)
  sample_ids <- paste0("S", 1:n_samples)

  data_grid <- expand.grid(MarkerName = marker_ids, SampleName = sample_ids)

  final_dataset <- data_grid %>%
    rowwise() %>%
    mutate(
      geno = sample(0:ploidy, 1),
      # Simulate allele signal as interpolation between "pure" profiles
      X = round(rnorm(1, mean = 280 - geno * ((280 - 20) / ploidy), sd = 10)),
      Y = round(rnorm(1, mean = 20 + geno * ((280 - 20) / ploidy), sd = 10)),
      R = X + Y,
      ratio = round(Y / R, 4),
      prob = round(runif(1, 0.85, 1.00), 2)
    ) %>%
    ungroup()

  sample_data <- final_dataset %>%
    select(MarkerName, SampleName, X, Y, R, ratio)

  geno_data <- final_dataset %>%
    select(MarkerName, SampleName, geno, prob)

  geno_pos <- data.frame(
    MarkerName = marker_ids,
    Chromosome = rep(c("1", "2", "3"), length.out = n_markers),
    Position = sort(sample(1e5:1e6, n_markers)),
    stringsAsFactors = FALSE
  )

  standardization_input <- inner_join(
    sample_data %>% select(MarkerName, SampleName, ratio) %>% rename(theta = ratio),
    geno_data %>% select(MarkerName, SampleName, geno),
    by = c("MarkerName", "SampleName")
  )

  list(
    sample_data = sample_data,
    geno_data = geno_data,
    geno_pos = geno_pos,
    standardization_input = standardization_input
  )
}
