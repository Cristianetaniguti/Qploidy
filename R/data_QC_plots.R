#' PCA Plot for Qploidy Standardized Data
#'
#' Performs Principal Component Analysis (PCA) on a Qploidy standardized object
#' using either the "R" or "z" column, and colors the samples based on a grouping
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
#' @importFrom dplyr left_join
#' @importFrom ggplot2 ggplot aes geom_point theme_minimal labs
#' @importFrom stats prcomp
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom scales percent_format
#'
#' @author Josue Chinchilla-Vargas
#'
#' @export
depth_pca_plot <- function(
    input,
    geno.pos = NULL,
    passport_file,
    group_column,
    sampleID_column,
    col2use = c("R","z")
) {

  #### Validate col2use ####
  col2use <- match.arg(col2use)
  message("Using column: ", col2use)


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

  # Handle input type
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
    # Check for Chr and Position
    if (all(c("Chr", "Position") %in% colnames(input))) {
      ttest_data <- input
    } else {
      # Need geno.pos to add Chr and Position
      if (is.null(geno.pos)) {
        stop("Input data.frame does not contain Chr and Position columns, and geno.pos was not provided.")
      }
      if (!all(c("MarkerName", "Chromosome", "Position") %in% colnames(geno.pos))) {
        stop("geno.pos must have columns: MarkerName, Chromosome, Position.")
      }
      m_idx <- match(input$MarkerName, geno.pos$MarkerName)
      input$Chr <- geno.pos$Chromosome[m_idx]
      input$Position <- geno.pos$Position[m_idx]
      ttest_data <- input
    }
  } else {
    stop("input must be a data.frame, a 'qploidy_standardization' object, or a valid file path.")
  }

  # Check if user chose col2use as "z" but it's not present in the data
  if (col2use == "z" && !("z" %in% colnames(ttest_data))) {
    stop("Column 'z' not found in the data. Please choose 'R' or ensure 'z' is present in the input data by running the standardization.")
  }

  #### Build genotype matrix ####
  gt_dp <- select(
    ttest_data,
    MarkerName,
    SampleName,
    all_of(col2use)
  ) %>%
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

  #### PCA ####

  # Remove zero variance markers
  marker_sd <- apply(gt_dp_num, 2, sd, na.rm = TRUE)
  gt_dp_num_filtered <- gt_dp_num[, marker_sd > 0, drop = FALSE]

  # Mean impute missing values
  for (j in seq_len(ncol(gt_dp_num_filtered))) {
    nas <- is.na(gt_dp_num_filtered[, j])
    if (any(nas)) {
      gt_dp_num_filtered[nas, j] <- mean(gt_dp_num_filtered[, j], na.rm = TRUE)
    }
  }

  pca_res <- prcomp(
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
  if(is.numeric(sampleID_column)) id <- colnames(passport)[sampleID_column] else id <- sampleID_column

  merged_df <- left_join(
    pca_scores,
    passport,
    by = c("ID" = id)  # assumes first column is ID
  )

  # Ensure the grouping column is a factor for discrete coloring
  merged_df[[group_column]] <- as.factor(merged_df[[group_column]])

  #### PCA Plots ####

  p1 <- ggplot(
    merged_df,
    aes(x = PC1, y = PC2, color = .data[[group_column]])
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      x = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y = paste0("PC2 (", round(pca_var[2] * 100, 2), "%)"),
      color = group_column,
      title = "PCA: PC1 vs PC2"
    )

  p2 <- ggplot(
    merged_df,
    aes(x = PC1, y = PC3, color = .data[[group_column]])
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
      x = paste0("PC1 (", round(pca_var[1] * 100, 2), "%)"),
      y = paste0("PC3 (", round(pca_var[3] * 100, 2), "%)"),
      color = group_column,
      title = "PCA: PC1 vs PC3"
    )

  p3 <- ggplot(
    merged_df,
    aes(x = PC2, y = PC3, color = .data[[group_column]])
  ) +
    geom_point(size = 3) +
    theme_minimal() +
    labs(
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


#' Heterozygosity Heatmap per Sample
#'
#' Computes the proportion of markers with BAF values within a heterozygosity
#' range (\code{het_range}) for each sample and displays the result as a
#' near-square heatmap. Hovering over a cell shows the sample name,
#' heterozygosity rate, and number of markers used (requires the \pkg{plotly}
#' package; a static \pkg{ggplot2} plot is returned if \pkg{plotly} is not
#' installed).
#'
#' @param input Accepts one of:
#'   1. A data.frame with at least columns \code{MarkerName}, \code{SampleName},
#'      and the column specified by \code{col2use}.
#'   2. An object of class \code{'qploidy_standardization'} (as returned by
#'      \code{standardize()}).
#'   3. A character string: path to a \code{.tsv[.gz]} file produced by
#'      \code{write_qploidy_standardization()}.
#' @param col2use Character. Column to use for computing heterozygosity.
#'   Either \code{"baf"} (default) or \code{"ratio"}.
#' @param het_range Numeric vector of length 2 giving the lower and upper
#'   bounds of \code{col2use} that define a heterozygous marker.
#'   Defaults to \code{c(0.2, 0.8)}.
#'
#' @return If \pkg{plotly} is available, a \code{plotly} htmlwidget with hover
#'   tooltips; otherwise a \code{ggplot2} object. Values are returned
#'   invisibly so the plot is still displayed when the function is called
#'   interactively.
#'
#' @examples
#' \dontrun{
#' plot_heterozygosity(my_standardization_object)
#' plot_heterozygosity("standardization.tsv.gz", het_range = c(0.3, 0.7))
#' plot_heterozygosity(my_standardization_object, col2use = "ratio")
#' }
#'
#' @importFrom dplyr group_by summarize
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn scale_x_continuous
#'   scale_y_continuous theme_minimal theme element_blank labs
#' @importFrom scales percent_format
#'
#' @export
plot_heterozygosity <- function(
    input,
    col2use = c("baf", "ratio"),
    het_range = c(0.2, 0.8)
) {
  ## --- validate arguments --------------------------------------------------
  col2use <- match.arg(col2use)

  if (!is.numeric(het_range) || length(het_range) != 2 ||
      het_range[1] >= het_range[2] ||
      het_range[1] < 0 || het_range[2] > 1) {
    stop("'het_range' must be a numeric vector of length 2 with 0 <= het_range[1] < het_range[2] <= 1.")
  }

  ## --- parse input (mirrors depth_pca_plot) --------------------------------
  if (inherits(input, "qploidy_standardization")) {
    tdata <- as.data.frame(input$data)
  } else if (is.character(input) && length(input) == 1 && file.exists(input)) {
    tdata <- as.data.frame(read_qploidy_standardization(input)$data)
  } else if (is.data.frame(input)) {
    tdata <- input
  } else {
    stop(paste(
      "'input' must be a data.frame with columns MarkerName, SampleName, and", col2use,
      "; a 'qploidy_standardization' object; or a valid file path."
    ))
  }

  required_cols <- c("MarkerName", "SampleName", col2use)
  missing_cols <- setdiff(required_cols, colnames(tdata))
  if (length(missing_cols) > 0)
    stop(sprintf(
      "Input data is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))

  ## --- compute per-sample heterozygosity -----------------------------------
  het_df <- tdata %>%
    group_by(SampleName) %>%
    summarize(
      het_ratio  = mean(.data[[col2use]] >= het_range[1] & .data[[col2use]] <= het_range[2], na.rm = TRUE),
      n_markers  = sum(!is.na(.data[[col2use]])),
      .groups    = "drop"
    )

  het_df <- het_df[order(het_df$het_ratio, decreasing = TRUE), ]

  ## --- lay out samples in a near-square grid -------------------------------
  n     <- nrow(het_df)
  ncols <- ceiling(sqrt(n))
  nrows <- ceiling(n / ncols)

  het_df$col_idx <- ((seq_len(n) - 1L) %% ncols) + 1L
  het_df$row_idx <- nrows - ((seq_len(n) - 1L) %/% ncols)  # row 1 = bottom (incomplete), highest row = top

  ## --- build ggplot heatmap ------------------------------------------------
  p <- ggplot(
    het_df,
    aes(
      x    = col_idx,
      y    = row_idx,
      fill = het_ratio,
      text = paste0(
        "Sample: ",          SampleName,
        "\nHeterozygosity: ", round(het_ratio * 100, 2), "%",
        "\nMarkers used: ",   n_markers
      )
    )
  ) +
    geom_tile(color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colours = c("#2166ac", "#fdae61", "#d7191c"),
      limits  = c(0, 1),
      name    = "Heterozygosity",
      labels  = percent_format(accuracy = 1)
    ) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      axis.title   = element_blank(),
      axis.text    = element_blank(),
      panel.grid   = element_blank()
    ) +
    labs(title = sprintf(
      "Heterozygosity per Sample  [%s in (%.2f, %.2f)]  n = %d",
      col2use, het_range[1], het_range[2], n
    ))

  ## --- add hover if plotly is available ------------------------------------
  if (requireNamespace("plotly", quietly = TRUE)) {
    out <- plotly::ggplotly(p, tooltip = "text")
    # xgap/ygap are only valid on heatmap traces; apply selectively to avoid warnings
    heatmap_idx <- which(vapply(out$x$data,
                                function(tr) isTRUE(tr$type == "heatmap"),
                                logical(1)))
    if (length(heatmap_idx) > 0)
      out <- plotly::style(out, xgap = 1.5, ygap = 1.5, traces = heatmap_idx)
  } else {
    message("Install the 'plotly' package for interactive hover tooltips.")
    out <- p
  }

  invisible(out)
}
