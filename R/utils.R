globalVariables(c("rnom", "X", "Y", "runif", "prob"))

#' Pascal Triangle for Expected Peaks Calculation
#'
#' This function generates the Pascal triangle for a given ploidy value. The Pascal triangle is used to define the expected peaks for each ploidy level, which can be useful in various genetic analyses.
#'
#' @param h An integer representing the ploidy value.
#'
#' @return A list where each element corresponds to a row of the Pascal triangle, up to the specified ploidy value.
#'
pascalTriangle <- function(h) {
  lapply(0:h, function(i) choose(i, 0:i))
}

#' Calculate the Statistical Mode
#'
#' This function returns the most frequent (modal) value in a vector.
#' If there are multiple values with the same highest frequency,
#' it returns the first one encountered.
#'
#' @param x A vector of numeric, character, or factor values.
#'
#' @return A single value representing the mode of the input vector.
#'
#'
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#' Adjust allele dosage using sequencing error and allele bias - function adapted from updog
#'
#' Adjusts allele dosage `p` by the sequencing error rate `eps`
#' and the allele bias `h`. This reproduces the behavior of the
#' `xi_fun()` from the updog package.
#'
#' @param p A numeric vector of allele dosages (e.g., 0 to 1).
#' @param eps A numeric vector (or scalar) of sequencing error rates.
#' @param h A numeric vector (or scalar) of allele bias values.
#'
#' @return A numeric vector of adjusted probabilities of the reference read.
#'
#' @author David Gerard
#' @keywords internal
#' @noRd
xi_fun <- function(p, eps, h) {
  n <- length(p)

  if (!(length(eps) %in% c(1, n))) {
    stop("xi_fun: eps must either have length 1 or the same length as p.")
  }

  if (!(length(h) %in% c(1, n))) {
    stop("xi_fun: h must either have length 1 or the same length as p.")
  }

  # Vector recycling
  eps <- rep(eps, length.out = n)
  h <- rep(h, length.out = n)

  # Compute adjusted probabilities
  xi <- p * (1 - h) + (1 - p) * h * (1 - eps)

  return(xi)
}

#' Perform a Sanity Check on a VCF File
#'
#' This function performs a series of checks on a VCF file to ensure its validity and integrity. It verifies the presence of required headers, columns, and data fields, and checks for common issues such as missing or malformed data.
#'
#' @param vcf_path A character string specifying the path to the VCF file. The file can be plain text or gzipped.
#' @param n_data_lines An integer specifying the number of data lines to sample for detailed checks. Default is 100.
#' @param max_markers An integer specifying the maximum number of markers allowed in the VCF file. Default is 10,000.
#' @param verbose A logical value indicating whether to print detailed messages during the checks. Default is FALSE.
#'
#' @return A list containing:
#' - `checks`: A named vector indicating the results of each check (TRUE or FALSE).
#' - `messages`: A data frame containing messages for each check, indicating success or failure.
#' - `duplicates`: A list containing any duplicated sample or marker IDs found in the VCF file.
#' - `ploidy_max`: The maximum ploidy detected from the genotype field, if applicable.
#'
#' @details The function performs the following checks:
#' - **VCF_header**: Verifies the presence of the `##fileformat` header.
#' - **VCF_columns**: Ensures required columns (`#CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`) are present.
#' - **max_markers**: Checks if the total number of markers exceeds the specified limit.
#' - **GT**: Verifies the presence of the `GT` (genotype) field in the FORMAT column.
#' - **allele_counts**: Checks for allele-level count fields (e.g., `AD`, `RA`, `AO`, `RO`).
#' - **samples**: Ensures sample/genotype columns are present.
#' - **chrom_info** and **pos_info**: Verifies the presence of `CHROM` and `POS` columns.
#' - **ref_alt**: Ensures `REF` and `ALT` fields contain valid nucleotide codes.
#' - **multiallelics**: Identifies multiallelic sites (ALT field with commas).
#' - **phased_GT**: Checks for phased genotypes (presence of `|` in the `GT` field).
#' - **duplicated_samples**: Checks for duplicated sample IDs.
#' - **duplicated_markers**: Checks for duplicated marker IDs.
#'
#' @importFrom stats setNames
#'
#' @export
vcf_sanity_check <- function(
    vcf_path,
    n_data_lines = 100,
    max_markers = 10000,
    verbose = FALSE) {
  if (!file.exists(vcf_path)) stop("File does not exist.")

  is_gz <- grepl("\\.gz$", vcf_path)
  con <- if (is_gz) gzfile(vcf_path, open = "rt") else file(vcf_path, open = "r")
  lines <- readLines(con, warn = FALSE)
  close(con)

  # --- Prepare result vector ---
  checks_names <- c(
    "VCF_header",
    "VCF_columns",
    "max_markers",
    "GT",
    "allele_counts",
    "samples",
    "chrom_info",
    "pos_info",
    "ref_alt",
    "multiallelics",
    "phased_GT",
    "duplicated_samples",
    "duplicated_markers",
    "mixed_ploidies"
  )
  checks <- setNames(rep(NA, length(checks_names)), checks_names)


  # Container for duplicated IDs
  duplicates <- list(
    duplicated_samples = character(0),
    duplicated_markers = character(0)
  )

  # --- Header checks ---
  header_lines <- grep("^##", lines, value = TRUE)
  if (!any(grepl("^##fileformat=VCFv", header_lines))) {
    checks["VCF_header"] <- FALSE
    if (verbose) warning("Missing ##fileformat header.")
  } else {
    checks["VCF_header"] <- TRUE
    if (verbose) cat("VCF header is present.\n")
  }

  # --- Column header line ---
  column_header_line <- grep("^#CHROM", lines, value = TRUE)
  if (length(column_header_line) != 1) stop("Missing or multiple #CHROM lines.")
  column_names <- unlist(strsplit(column_header_line, "\t"))
  has_genotypes <- length(column_names) > 8

  required_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  checks["VCF_columns"] <- all(required_columns %in% column_names[1:8])
  if (checks["VCF_columns"]) {
    if (verbose) cat("Required VCF columns are present.\n")
  } else {
    if (verbose) warning("Missing one or more required VCF columns.")
  }

  # --- Total marker count ---
  data_line_indices <- grep("^[^#]", lines)
  total_markers <- length(data_line_indices)
  if (verbose) cat(sprintf("Total markers (data rows): %d\n", total_markers))

  checks["max_markers"] <- total_markers <= max_markers
  if (!checks["max_markers"]) {
    warning(sprintf("More than %d markers found. Consider subsampling.", max_markers))
  }

  # --- Check for duplicated marker IDs ---
  if (total_markers > 0) {
    marker_ids <- sapply(lines[data_line_indices], function(line) {
      fields <- strsplit(line, "\t")[[1]]
      if (length(fields) >= 3) fields[3] else NA
    })
    marker_ids <- marker_ids[!is.na(marker_ids)]
    duplicated_markers <- marker_ids[duplicated(marker_ids)]
    duplicates$duplicated_markers <- duplicated_markers
    if (length(duplicated_markers) > 0) {
      checks["duplicated_markers"] <- TRUE
      if (verbose) warning("Duplicated marker IDs found: ", paste(head(duplicated_markers, 10), collapse = ", "), "...")
    } else {
      if (verbose) cat("No duplicated marker IDs.\n")
      checks["duplicated_markers"] <- FALSE
    }
  }

  # --- FORMAT field checks (GT, AD, etc.) ---
  if (has_genotypes) {
    checks["samples"] <- TRUE
    sample_indices <- head(data_line_indices, n_data_lines)
    format_fields <- character()

    for (i in seq_along(sample_indices)) {
      fields <- strsplit(lines[sample_indices[i]], "\t")[[1]]
      if (length(fields) >= 9) {
        format_keys <- unlist(strsplit(fields[9], ":"))
        format_fields <- union(format_fields, format_keys)
      }
    }

    # GT check
    checks["GT"] <- "GT" %in% format_fields
    if (verbose) {
      if (checks["GT"]) cat("FORMAT field 'GT' (genotype) is present.\n") else warning("FORMAT field 'GT' is missing.")
    }
    # Allele counts check
    support_fields <- c("AD", "RA", "AO", "RO", "NR", "NV", "SB", "F1R2", "F2R1")
    checks["allele_counts"] <- any(support_fields %in% format_fields)
    if (checks["allele_counts"]) {
      if (verbose) {
        cat(sprintf(
          "Allele count FORMAT field(s) found: %s\n",
          paste(intersect(support_fields, format_fields), collapse = ", ")
        ))
      }
    } else {
      warning(" No allele-level count fields found (e.g., AD, RA, AO, RO).")
    }

    # Optional: phased GT (presence of '|' instead of '/' in genotypes)
    phased_lines <- grep("\\|", lines[sample_indices])
    checks["phased_GT"] <- length(phased_lines) > 0

    # --- Check for duplicated sample names ---
    sample_names <- column_names[10:length(column_names)]
    duplicated_samples <- sample_names[duplicated(sample_names)]
    duplicates$duplicated_samples <- duplicated_samples
    if (length(duplicated_samples) > 0) {
      if (verbose) warning("Duplicated sample names found: ", paste(duplicated_samples, collapse = ", "))
      checks["duplicated_samples"] <- TRUE
    } else {
      if (verbose) cat("No duplicated sample names.\n")
      checks["duplicated_samples"] <- FALSE
    }
  } else {
    checks["samples"] <- FALSE
    checks["GT"] <- FALSE
    checks["allele_counts"] <- FALSE
    checks["phased_GT"] <- FALSE
    checks["duplicated_samples"] <- FALSE

    warning("No sample/genotype columns found.")
  }

  # --- Ploidy inference (based on GT) ---
  ploidy_max <- NA # default if GT not found or no valid genotypes

  if (checks["GT"]) {
    ploidy_values <- c()

    for (i in seq_along(sample_indices)) {
      fields <- strsplit(lines[sample_indices[i]], "\t")[[1]]

      # Skip if not enough fields
      if (length(fields) < 10) next

      format_keys <- unlist(strsplit(fields[9], ":"))
      gt_index <- which(format_keys == "GT")

      # Loop through samples (from column 10 onward)
      for (sample_field in fields[10:length(fields)]) {
        sample_fields <- unlist(strsplit(sample_field, ":"))
        if (length(sample_fields) >= gt_index) {
          gt_raw <- sample_fields[gt_index]
          if (grepl("[/|]", gt_raw)) {
            ploidy <- length(unlist(strsplit(gt_raw, "[/|]")))
            ploidy_values <- c(ploidy_values, ploidy)
          }
        }
      }
    }

    if (length(ploidy_values) > 0) {
      ploidy_max <- max(ploidy_values, na.rm = TRUE)
      if (verbose) cat("Highest ploidy detected from GT field:", ploidy_max, "\n")
    }
    if (length(ploidy_values) > 1) {
      checks["mixed_ploidies"] <- TRUE
      if (verbose) cat("Mixed ploidies\n")
    } else {
      checks["mixed_ploidies"] <- FALSE
    }
  }

  # --- CHROM and POS column checks ---
  checks["chrom_info"] <- "#CHROM" %in% column_names
  checks["pos_info"] <- "POS" %in% column_names

  if (checks["chrom_info"] && checks["pos_info"]) {
    if (verbose) cat("Both CHROM and POS columns are present.\n")
  } else {
    if (!checks["chrom_info"]) warning(" Column 'CHROM' is missing.")
    if (!checks["pos_info"]) warning(" Column 'POS' is missing.")
  }

  # --- REF/ALT basic check on sample rows ---
  sample_lines <- lines[head(data_line_indices, n_data_lines)]
  ref_alt_valid <- sapply(sample_lines, function(line) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) >= 5) {
      ref <- fields[4]
      alt <- fields[5]
      grepl("^[ACGTN]+$", ref) && grepl("^[ACGTN.,<>]+$", alt)
    } else {
      FALSE
    }
  })
  checks["ref_alt"] <- all(ref_alt_valid)

  # --- Multiallelic site check (ALT with ',' separator) ---
  multiallelic_flags <- grepl(",", sapply(sample_lines, function(line) strsplit(line, "\t")[[1]][5]))
  checks["multiallelics"] <- any(multiallelic_flags)

  # --- Compile messages ---

  # Messages in case of failure
  messages <- data.frame(
    "VCF_header" = c(
      "VCF header is missing. Please check the file format.",
      "VCF header is present."
    ),
    "VCF_columns" = c(
      "Required VCF columns are missing. Please check the file format.",
      "Required VCF columns are present."
    ),
    "max_markers" = c(
      "More than 10,000 markers found. Consider subsampling or running in HPC.",
      "Less than maximum number of markers found."
    ),
    "GT" = c(
      "Genotype information is not available in the VCF file.",
      "Genotype information is available in the VCF file."
    ),
    "allele_counts" = c(
      "Allele counts are not available in the VCF file.",
      "Allele counts are available in the VCF file."
    ),
    "samples" = c(
      "Sample information is not available in the VCF file.",
      "Sample information is available in the VCF file."
    ),
    "chrom_info" = c(
      "Chromosome information is not available in the VCF file.",
      "Chromosome information is available in the VCF file."
    ),
    "pos_info" = c(
      "Position information is not available in the VCF file.",
      "Position information is available in the VCF file."
    ),
    "ref_alt" = c(
      "REF/ALT fields contain invalid nucleotide codes.",
      "REF/ALT fields are valid."
    ),
    "multiallelics" = c(
      "Multiallelic sites not found in the VCF file.",
      "Multiallelic sites found in the VCF file."
    ),
    "phased_GT" = c(
      "Phased genotypes (|) are not present in the VCF file.",
      "Phased genotypes (|) are present in the VCF file."
    ),
    "duplicated_samples" = c(
      "No duplicated sample IDs found.",
      paste("Duplicated sample IDs found: ", paste(duplicates$duplicated_samples, collapse = ", "))
    ),
    "duplicated_markers" = c(
      "No duplicated marker IDs found.",
      paste("Duplicated marker IDs found: ", paste(duplicates$duplicated_markers, collapse = ", "))
    ),
    "mixed_ploidies" = c(
      "Mixed ploidies detected.",
      "No mixed ploidies detected."
    )
  )
  rownames(messages) <- c("false", "true")

  # --- Done ---
  if (verbose) cat("Sanity check complete.\n")
  return(structure(
    list(
      checks = checks, messages = messages,
      duplicates = duplicates, ploidy_max = ploidy_max
    ),
    class = "vcf_sanity_check"
  ))
}
