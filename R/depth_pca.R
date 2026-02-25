#' PCA Plot for Qploidy Standardized Data
#'
#' Performs Principal Component Analysis (PCA) on a Qploidy standardized object
#' using either the "R" or "z" column, and colors the samples based on a grouping
#' variable from a passport file (CSV or TXT).
#'
#' @param qploidy_object Character. Path to the Qploidy standardization file
#'   (e.g., "standardization.tsv.gz") or Qploidy object.
#' @param passport_file Character. Path to a CSV or TXT file containing sample
#'   metadata. The first column must contain the sample IDs matching the Qploidy
#'   object.
#' @param group_column Character. The name of the column in the passport file
#'   used to color/group the PCA plots.
#' @param col2use Character. Either "R" or "z". Which Qploidy column to use for
#'   the PCA. Defaults to c("R","z").
#'
#' @return Prints three ggplot2 PCA plots: PC1 vs PC2, PC1 vs PC3, and PC2 vs PC3.
#'   Returns a list of these plots invisibly.
#'
#' @examples
#' depth_pca_plot(
#'   qploidy_object = "standardization.tsv.gz",
#'   passport_file = "pca_ids.csv",
#'   group_column = "Plate",
#'   col2use = "R"
#' )
#'
#' @export
depth_pca_plot <- function(
    qploidy_object,
    passport_file,
    group_column,
    col2use = c("R","z")
) {

  #### Validate col2use ####
  col2use <- base::match.arg(col2use)
  base::message("Using column: ", col2use)

  #### Read passport (csv or txt) ####
  if (grepl("\\.csv$", passport_file, ignore.case = TRUE)) {
    passport <- read.csv(passport_file, stringsAsFactors = FALSE)
  } else if (grepl("\\.txt$", passport_file, ignore.case = TRUE)) {
    passport <- read.delim(passport_file, stringsAsFactors = FALSE)
  } else {
    stop("passport_file must be .csv or .txt")
  }

  if (!group_column %in% colnames(passport)) {
    stop("group_column not found in passport file.")
  }

  #### Read standardized data ####
  ttest_data <- Qploidy::read_qploidy_standardization(qploidy_object)$data
  ttest_data$SampleName <- base::sub("^X", "", ttest_data$SampleName)

  #### Build genotype matrix ####
  gt_dp <- dplyr::select(
    ttest_data,
    MarkerName,
    SampleName,
    dplyr::all_of(col2use)
  ) %>%
    tidyr::pivot_wider(
      names_from  = SampleName,
      values_from = dplyr::all_of(col2use)
    ) %>%
    tibble::column_to_rownames(var = "MarkerName") %>%
    base::as.matrix()

  gt_dp_num <- t(base::matrix(
    as.numeric(gt_dp),
    nrow = nrow(gt_dp),
    ncol = ncol(gt_dp),
    dimnames = dimnames(gt_dp)
  ))

  #### PCA ####

  # Remove zero variance markers
  marker_sd <- apply(gt_dp_num, 2, stats::sd, na.rm = TRUE)
  gt_dp_num_filtered <- gt_dp_num[, marker_sd > 0, drop = FALSE]

  # Mean impute missing values
  for (j in seq_len(ncol(gt_dp_num_filtered))) {
    nas <- is.na(gt_dp_num_filtered[, j])
    if (any(nas)) {
      gt_dp_num_filtered[nas, j] <- mean(gt_dp_num_filtered[, j], na.rm = TRUE)
    }
  }

  pca_res <- stats::prcomp(
    gt_dp_num_filtered,
    center = TRUE,
    scale. = TRUE
  )

  pca_var <- (pca_res$sdev^2) / sum(pca_res$sdev^2)

  pca_scores <- data.frame(
    ID  = rownames(pca_res$x),
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    PC3 = pca_res$x[, 3],
    stringsAsFactors = FALSE
  )

  #### Merge passport ####
  merged_df <- dplyr::left_join(
    pca_scores,
    passport,
    by = c("ID" = colnames(passport)[1])  # assumes first column is ID
  )

  # Ensure the grouping column is a factor for discrete coloring
  merged_df[[group_column]] <- as.factor(merged_df[[group_column]])

  #### PCA Plots ####

  p1 <- ggplot2::ggplot(
    merged_df,
    ggplot2::aes(x = PC1, y = PC2, color = .data[[group_column]])
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(pca_var[2] * 100, 2), "%)"),
      color = group_column,
      title = "PCA: PC1 vs PC2"
    )

  p2 <- ggplot2::ggplot(
    merged_df,
    ggplot2::aes(x = PC1, y = PC3, color = .data[[group_column]])
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y = paste0("PC3 (", round(pca_var[3] * 100, 2), "%)"),
      color = group_column,
      title = "PCA: PC1 vs PC3"
    )

  p3 <- ggplot2::ggplot(
    merged_df,
    ggplot2::aes(x = PC2, y = PC3, color = .data[[group_column]])
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0("PC2 (", round(pca_var[2] * 100, 2), "%)"),
      y = paste0("PC3 (", round(pca_var[3] * 100, 2), "%)"),
      color = group_column,
      title = "PCA: PC2 vs PC3"
    )

  #### Print plots ####
  print(p1)
  print(p2)
  print(p3)

  invisible(list(p1, p2, p3))
}
