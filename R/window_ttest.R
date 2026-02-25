#' Window-level  Depth Comparison Plot
#'
#' Generates a scatter plot comparing normalized SNP depth of selected samples against the population mean in genomic windows.
#' Significance is calculated per window using unpaired t-tests and indicated in the legend.
#'
#' @param input Accepts one of:
#'   1. A data.frame with columns: MarkerName, SampleName, Chr, Position, R. Column with z is optional but required if col2test = "z".
#'   2. An object of class 'qploidy_standardization' (from Qploidy standardization functions)
#'   3. A character string: path to a .tsv[.gz] file with standardized Qploidy data
#'   The function will validate and process the input according to its type.
#' @param selected_samples Character vector of sample names to compare against the population.
#' @param col2test Column name to use for depth comparison. Options are "R" or "z".
#' @param window_size Size of genomic windows in base pairs (default 1e6).
#'
#' @details
#' The function reads standardized SNP depth data from the input, normalizes it per sample,
#' and aggregates values within windows across chromosomes. It performs unpaired t-tests per window
#' to compare the selected sample against the population mean. Significance is categorized by direction
#' (higher or lower than population) and adjusted p-value (FDR).
#'
#' The resulting ggplot2 scatter plot shows normalized SNP depth per window along the x-axis (megabase),
#' colored by significance, with facets for each sample and chromosome.
#'
#' @return A ggplot2 object representing window-level SNP depth comparison.
#'
#' @examples
#' window_level_test_plot_qploidy(
#'   qploidy_object = "standardization.tsv.gz",
#'   selected_samples = "J432",
#'   col2test = "z",
#'   window_size = 1e6
#' )
#'
#' @importFrom dplyr select distinct bind_rows case_when
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames
#' @importFrom ggplot2 ggplot aes geom_point facet_grid scale_color_manual scale_y_continuous theme_bw labs
#' @importFrom stats t.test p.adjust
#' @importFrom scales percent_format
#' @export
window_level_test_plot_qploidy <- function(
    input,
    selected_samples,
    col2test = c("R","z"),
    window_size = 1e6
) {

  # Validate col2test
  col2test <- match.arg(col2test)
  message("Using column: ", col2test)

 # Handle input type
  if (inherits(input, "qploidy_standardization")) {
    ttest_data <- input$data
  } else if (is.character(input) && length(input) == 1 && file.exists(input)) {
    ttest_data <- read_qploidy_standardization(input)$data
  } else if (is.data.frame(input)) {
    required_cols <- c("MarkerName", "SampleName", "Chr", "Position", "R")
    missing_cols <- setdiff(required_cols, colnames(input))
    if (length(missing_cols) > 0) {
      stop(sprintf("Input data.frame is missing required columns: %s", paste(missing_cols, collapse = ", ")))
    }
    ttest_data <- input
  } else {
    stop("input must be a data.frame, a 'qploidy_standardization' object, or a valid file path.")
  }

  # Check if user chose col2test as "z" but it's not present in the data
  if (col2test == "z" && !("z" %in% colnames(ttest_data))) {
    stop("Column 'z' not found in the data. Please choose 'R' or ensure 'z' is present in the input data by running the standardization.")
  }

  # Build genotype matrix
  gt_dp <- select(ttest_data, MarkerName, SampleName, all_of(col2test)) %>%
    pivot_wider(
      names_from  = SampleName,
      values_from = all_of(col2test)
    ) %>%
    column_to_rownames(var = "MarkerName") %>%
    as.matrix()

  gt_dp_num <- matrix(
    as.numeric(gt_dp),
    nrow = nrow(gt_dp),
    ncol = ncol(gt_dp),
    dimnames = dimnames(gt_dp)
  )

  # Normalization
  snp_prop <- sweep(
    gt_dp_num,
    2,
    colSums(gt_dp_num, na.rm = TRUE),
    "/"
  ) * 100

  all_samples <- colnames(snp_prop)

  # Window labels
  marker_info <- select(ttest_data, MarkerName, Chr, Position) %>%
    distinct(MarkerName, .keep_all = TRUE)

  chroms    <- marker_info$Chr
  positions <- marker_info$Position

  window_id    <- floor((positions - 1) / window_size) + 1
  window_label <- paste0(chroms, "_MB", window_id)
  unique_windows <- unique(window_label)

  stopifnot(length(window_label) == nrow(snp_prop))  # sanity check

  # Window-level test
  window_plot_df <- lapply(selected_samples, function(selected_sample) {

    other_samples <- setdiff(all_samples, selected_sample)

    # Population mean
    snp_pop_mean <- rowMeans(snp_prop[, other_samples], na.rm = TRUE)

    # Aggregate per window
    df <- data.frame(
      WINDOW     = unique_windows,
      CHROM      = sub("_MB.*", "", unique_windows),
      Megabase   = as.integer(sub(".*_MB", "", unique_windows)),
      Sample     = selected_sample,
      SampleProp = vapply(unique_windows, function(win) {
        mean(snp_prop[window_label == win, selected_sample], na.rm = TRUE)
      }, numeric(1)),
      PopAvgProp = vapply(unique_windows, function(win) {
        mean(snp_pop_mean[window_label == win], na.rm = TRUE)
      }, numeric(1))
    )

    # Unpaired t-test per window
    df$p_raw <- vapply(unique_windows, function(win) {
      idx <- window_label == win
      x_sel <- snp_prop[idx, selected_sample]
      x_pop <- snp_pop_mean[idx]
      ok <- is.finite(x_sel) & is.finite(x_pop)
      if (sum(ok) >= 5) t.test(x_sel[ok], x_pop[ok], paired = FALSE)$p.value else NA_real_
    }, numeric(1))

    df$p_adj <- p.adjust(df$p_raw, method = "fdr")
    df$delta <- df$SampleProp - df$PopAvgProp

    # Significance categories
    df$SigDir <- case_when(
      is.na(df$p_adj) | df$p_adj >= 0.05 ~ "Not significant",
      df$delta < 0 & df$p_adj < 0.001    ~ "Lower: P < 0.001",
      df$delta < 0 & df$p_adj < 0.01     ~ "Lower: P < 0.01",
      df$delta < 0                        ~ "Lower: P < 0.05",
      df$delta > 0 & df$p_adj < 0.001    ~ "Higher: P < 0.001",
      df$delta > 0 & df$p_adj < 0.01     ~ "Higher: P < 0.01",
      TRUE                               ~ "Higher: P < 0.05"
    )

    df
  }) %>% bind_rows()

  # Plotting
  window_plot_long <- pivot_longer(
    window_plot_df,
    cols = c(SampleProp, PopAvgProp),
    names_to = "ID",
    values_to = "DepthProp"
  )

  window_plot_long$LegendGroup <- case_when(
    window_plot_long$ID == "PopAvgProp" ~ "Population mean",
    window_plot_long$ID == "SampleProp" & window_plot_long$SigDir == "Not significant" ~ "Sample not significant",
    TRUE ~ window_plot_long$SigDir
  )

  legend_colors <- c(
    "Population mean"        = "grey70",
    "Sample not significant" = "gold",
    "Lower: P < 0.05"        = "#6baed6",
    "Lower: P < 0.01"        = "#3182bd",
    "Lower: P < 0.001"       = "#08519c",
    "Higher: P < 0.05"       = "#fc9272",
    "Higher: P < 0.01"       = "#fb6a4a",
    "Higher: P < 0.001"      = "#cb181d"
  )

  ggplot(
    window_plot_long,
    aes(x = Megabase, y = DepthProp, color = LegendGroup)
  ) +
    geom_point(size = 0.75, alpha = 0.85) +
    facet_grid(Sample ~ CHROM, scales = "free_x") +
    scale_color_manual(values = legend_colors, name = "Legend") +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    theme_bw() +
    labs(
      x = "Window (Mb)",
      y = "Normalized SNP depth (%)",
      title = "Window-level SNP depth: sample vs population"
    )
}
