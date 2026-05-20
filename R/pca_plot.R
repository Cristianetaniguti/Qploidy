#' PCA Plot for Qploidy Standardized Data
#'
#' Performs Principal Component Analysis (PCA) on a Qploidy standardized object
#' using either the "R", "z", or "ratio" column, and colors the samples based on a grouping
#' variable from a passport file (CSV or TXT).
#'
#' @param input Accepts one of:
#'   1. A data.frame with columns: MarkerName, SampleName, Chr, Position, R. Column with z is optional but required if col2use = "z".
#'   2. An object of class 'qploidy_standardization' (from Qploidy standardization functions)
#'   3. A character string: path to a .tsv[.gz] file with standardized Qploidy data
#'   The function will validate and process the input according to its type.
#' @param passport_file Either a data.frame with sample metadata (first column must contain sample IDs matching the Qploidy object), or a character path to a CSV or TXT file containing sample metadata.
#' @param geno.pos Optional. A data.frame with columns: MarkerName, Chromosome, Position. Required if input data.frame does not contain Chr and Position columns.
#' @param group_column Character. The name of the column in the passport file
#'   used to color/group the PCA plots.
#' @param sampleID_column character or numeric indicating the column name or index where the sample name is located. Sample name must match `input`
#' @param col2use Character. Either "R", "z", or "ratio". Which Qploidy column to use for
#'   the PCA. Defaults to c("R","z","ratio").
#' @param plot_title Optional. A character string appended to the plot title after the PC pair label
#'   (e.g., "PC1 vs PC2 - My Title"). If NULL (default), the title will include the column used
#'   (e.g., "PCA: PC1 vs PC2 | R").
#' @param samples Optional character vector of sample IDs to include in the
#'   PCA. Values must match entries in the `SampleName` column of the input
#'   data. If `NULL` (default), all samples are used. Any requested sample not
#'   found in the data triggers a warning.
#' @param palette Character string naming a built-in colour palette, or a
#'   character vector of colours to use directly. Built-in options:
#'   \describe{
#'     \item{`"auto"` (default)}{ggplot2 default discrete colour scale.}
#'     \item{`"alphabet"`}{26 highly-distinguishable colours (hand-tuned for
#'       large group counts).}
#'     \item{`"polychrome"`}{24 perceptually-distinct colours from the
#'       Polychrome palette, suitable for 20+ groups.}
#'     \item{`"okabe_ito"`}{8-colour Okabe-Ito palette — the recommended
#'       colorblind-friendly choice. Colours are recycled if there are more
#'       than 8 groups.}
#'   }
#'   When a named palette has fewer colours than the number of groups the
#'   colours are recycled with a warning.
#'
#' @return Prints three ggplot2 PCA plots: PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3.
#'   Returns a list of these plots invisibly.
#'
#' @examples
#' \dontrun{
#' pca_plot(
#'   input = "standardization.tsv.gz",
#'   passport_file = "pca_ids.csv",
#'   group_column = "Plate",
#'   col2use = "R",
#'   plot_title = "My Experiment"
#' )
#' }
#'
#' @importFrom dplyr left_join select
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs scale_color_hue scale_color_manual
#' @importFrom stats prcomp
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom scales percent_format
#'
#' @author Josue Chinchilla-Vargas
#'
#' @export
pca_plot <- function(
    input,
    geno.pos = NULL,
    passport_file,
    group_column,
    sampleID_column,
    col2use = c("R", "z", "ratio"),
    plot_title = NULL,
    palette = "auto",
    samples = NULL
) {
  #### Validate col2use ####
  col2use <- match.arg(col2use)
  message("Using column: ", col2use)

  # Resolve plot title suffix
  title_suffix <- if (!is.null(plot_title)) plot_title else col2use

  #### Read passport (data.frame, csv or txt) ####
  if (is.data.frame(passport_file)) {
    passport <- passport_file
  } else if (is.character(passport_file) && length(passport_file) == 1) {
    if (grepl("\\.csv$", passport_file, ignore.case = TRUE)) {
      passport <- read.csv(passport_file, stringsAsFactors = FALSE)
    } else if (grepl("\\.txt$", passport_file, ignore.case = TRUE)) {
      passport <- read.delim(passport_file, stringsAsFactors = FALSE)
    } else {
      stop("passport_file must be a data.frame, .csv or .txt")
    }
  } else {
    stop("passport_file must be a data.frame, .csv or .txt")
  }

  if (!group_column %in% colnames(passport)) {
    stop("group_column not found in passport file.")
  }

  #### Handle input type ####
  if (inherits(input, "qploidy_standardization")) {
    ttest_data <- input$data
  } else if (is.character(input) && length(input) == 1 && file.exists(input)) {
    ttest_data <- read_qploidy_standardization(input)$data
  } else if (is.data.frame(input)) {
    required_cols <- c("MarkerName", "SampleName", "R")
    missing_cols <- setdiff(required_cols, colnames(input))
    if (length(missing_cols) > 0) {
      stop(sprintf("Input data.frame is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    if (all(c("Chr", "Position") %in% colnames(input))) {
      ttest_data <- input
    } else {
      if (is.null(geno.pos)) {
        stop("Input data.frame does not contain Chr and Position columns, and geno.pos was not provided.")
      }
      if (!all(c("MarkerName", "Chromosome", "Position") %in% colnames(geno.pos))) {
        stop("geno.pos must have columns: MarkerName, Chromosome, Position.")
      }
      m_idx <- match(input$MarkerName, geno.pos$MarkerName)
      input$Chr      <- geno.pos$Chromosome[m_idx]
      input$Position <- geno.pos$Position[m_idx]
      ttest_data <- input
    }
  } else {
    stop("input must be a data.frame, a 'qploidy_standardization' object, or a valid file path.")
  }

  #### Validate chosen column exists ####
  if (!is.null(samples)) {
    missing_samples <- setdiff(samples, unique(ttest_data$SampleName))
    if (length(missing_samples) > 0)
      warning("The following requested samples were not found and will be ignored: ",
              paste(missing_samples, collapse = ", "))
    ttest_data <- ttest_data[ttest_data$SampleName %in% samples, , drop = FALSE]
    if (nrow(ttest_data) == 0)
      stop("No data remaining after filtering to the requested samples.")
  }

  if (col2use == "z" && !("z" %in% colnames(ttest_data))) {
    stop("Column 'z' not found in the data. Please choose 'R' or ensure 'z' is present.")
  }
  if (col2use == "ratio" && !("ratio" %in% colnames(ttest_data))) {
    stop("Column 'ratio' not found in the data. Please choose 'R' or ensure 'ratio' is present.")
  }

  #### Build matrix and run PCA ####
  if (col2use == "ratio") {

    message("Building genomic relationship matrix from ratio column via AGHmatrix::Gmatrix ...")

    test_matrix <- ttest_data %>%
      select(SampleName, MarkerName, ratio) %>%
      pivot_wider(
        names_from  = MarkerName,
        values_from = ratio,
        values_fill = list(ratio = NA)
      ) %>%
      column_to_rownames(var = "SampleName") %>%
      as.matrix()

    if (!requireNamespace("AGHmatrix", quietly = TRUE))
      stop("Package 'AGHmatrix' is required when col2use = 'ratio'. ",
           "Install it with: install.packages('AGHmatrix')")

    ratio_matrix <- AGHmatrix::Gmatrix(
      test_matrix,
      ratio        = TRUE,
      missingValue = NA,
      ratio.check  = TRUE
    )

    pca_res <- prcomp(ratio_matrix, center = TRUE, scale. = FALSE)

  } else {

    gt_dp <- select(ttest_data, MarkerName, SampleName, all_of(col2use)) %>%
      pivot_wider(
        names_from  = SampleName,
        values_from = all_of(col2use)
      ) %>%
      column_to_rownames(var = "MarkerName") %>%
      as.matrix()

    gt_dp_num <- t(matrix(
      as.numeric(gt_dp),
      nrow = nrow(gt_dp),
      ncol = ncol(gt_dp),
      dimnames = dimnames(gt_dp)
    ))

    marker_sd <- apply(gt_dp_num, 2, sd, na.rm = TRUE)
    # Use explicit TRUE/FALSE (not NA) to avoid NA columns being silently included
    gt_dp_num_filtered <- gt_dp_num[, !is.na(marker_sd) & marker_sd > 0, drop = FALSE]

    for (j in seq_len(ncol(gt_dp_num_filtered))) {
      nas <- is.na(gt_dp_num_filtered[, j])
      if (any(nas)) {
        gt_dp_num_filtered[nas, j] <- mean(gt_dp_num_filtered[, j], na.rm = TRUE)
      }
    }

    # Replace any non-finite values (NaN, Inf, -Inf) introduced by imputation
    gt_dp_num_filtered[!is.finite(gt_dp_num_filtered)] <- NA

    # Final comprehensive filter: remove columns with any NA or zero/non-finite sd
    keep <- apply(gt_dp_num_filtered, 2, function(col) {
      !anyNA(col) && is.finite(sd(col)) && sd(col) > 0
    })
    n_removed <- sum(!keep)
    if (n_removed > 0)
      message(n_removed, " marker(s) removed after imputation due to zero variance or non-finite values.")
    gt_dp_num_filtered <- gt_dp_num_filtered[, keep, drop = FALSE]

    pca_res <- prcomp(gt_dp_num_filtered, center = TRUE, scale. = TRUE)
  }

  #### Shared PCA summary ####
  pca_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

  pca_scores <- data.frame(
    ID  = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    stringsAsFactors = FALSE
  )

  #### Merge passport ####
  if (is.numeric(sampleID_column)) {
    id <- colnames(passport)[sampleID_column]
  } else {
    id <- sampleID_column
  }

  merged_df <- left_join(pca_scores, passport, by = c("ID" = id))
  merged_df[[group_column]] <- as.factor(merged_df[[group_column]])

  #### Resolve colour palette ####
  .palettes <- list(
    alphabet = c(
      "#AA0DFE", "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356",
      "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", "#DEA0FD",
      "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", "#90AD1C", "#F6222E",
      "#1CFFCE", "#2ED9FF", "#B10DA1", "#C075A6", "#FC1CBF", "#B00068",
      "#FBE426", "#FA0087"
    ),
    polychrome = c(
      "#5A5156", "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE",
      "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", "#DEA0FD",
      "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", "#1C8356", "#85660D",
      "#B10DA1", "#FBE426", "#1CBE4F", "#FA0087", "#FC1CBF", "#F7E1A0"
    ),
    okabe_ito = c(
      "#E69F00", "#56B4E9", "#009E73", "#F0E442",
      "#0072B2", "#D55E00", "#CC79A7", "#000000"
    )
  )

  n_groups <- nlevels(merged_df[[group_column]])

  if (identical(palette, "auto")) {
    colour_scale <- scale_color_hue()
  } else if (is.character(palette) && length(palette) == 1 && palette %in% names(.palettes)) {
    pal_cols <- .palettes[[palette]]
    if (n_groups > length(pal_cols)) {
      warning(sprintf(
        "Palette '%s' has %d colours but there are %d groups; colours will be recycled.",
        palette, length(pal_cols), n_groups
      ))
      pal_cols <- rep_len(pal_cols, n_groups)
    }
    colour_scale <- scale_color_manual(values = pal_cols)
  } else if (is.character(palette) && length(palette) > 1) {
    # Custom vector supplied directly
    if (n_groups > length(palette)) {
      warning(sprintf(
        "Supplied palette has %d colours but there are %d groups; colours will be recycled.",
        length(palette), n_groups
      ))
      palette <- rep_len(palette, n_groups)
    }
    colour_scale <- scale_color_manual(values = palette)
  } else {
    stop("'palette' must be one of 'auto', 'alphabet', 'polychrome', 'okabe_ito', or a character vector of colours.")
  }

  #### PCA Plots ####
  p1 <- ggplot(merged_df, aes(x = PC1, y = PC2, color = .data[[group_column]])) +
    geom_point(size = 3) +
    colour_scale +
    theme_minimal() +
    labs(
      x     = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y     = paste0("PC2 (", round(pca_var[2] * 100, 2), "%)"),
      color = group_column,
      title = paste0("PCA: PC1 vs PC2 | ", title_suffix)
    )

  p2 <- ggplot(merged_df, aes(x = PC1, y = PC3, color = .data[[group_column]])) +
    geom_point(size = 3) +
    colour_scale +
    theme_minimal() +
    labs(
      x     = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y     = paste0("PC3 (", round(pca_var[3] * 100, 2), "%)"),
      color = group_column,
      title = paste0("PCA: PC1 vs PC3 | ", title_suffix)
    )

  p3 <- ggplot(merged_df, aes(x = PC2, y = PC3, color = .data[[group_column]])) +
    geom_point(size = 3) +
    colour_scale +
    theme_minimal() +
    labs(
      x     = paste0("PC2 (", round(pca_var[2] * 100, 2), "%)"),
      y     = paste0("PC3 (", round(pca_var[3] * 100, 2), "%)"),
      color = group_column,
      title = paste0("PCA: PC2 vs PC3 | ", title_suffix)
    )

  return(list(PC1xPC2=p1, PC1xPC3 = p2, PC2xPC3 = p3))
}
