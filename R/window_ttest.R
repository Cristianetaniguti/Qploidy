library(Qploidy)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

#### function ####
window_level_test_plot_qploidy <- function(
    qploidy_object,
    selected_samples,
    col2test = c("R","z"),
    window_size = 1e6
) {


  # Validate col2test
  col2test <- base::match.arg(col2test)
  base::message("Using column: ", col2test)


  # Read standardized data
  ttest_data <- Qploidy::read_qploidy_standardization(qploidy_object)$data
  ttest_data$SampleName <- base::sub("^X", "", ttest_data$SampleName)


  # Build genotype matrix
  gt_dp <- dplyr::select(ttest_data, MarkerName, SampleName, dplyr::all_of(col2test)) %>%
    tidyr::pivot_wider(
      names_from  = SampleName,
      values_from = dplyr::all_of(col2test)
    ) %>%
    tibble::column_to_rownames(var = "MarkerName") %>%
    base::as.matrix()

  gt_dp_num <- base::matrix(
    as.numeric(gt_dp),
    nrow = nrow(gt_dp),
    ncol = ncol(gt_dp),
    dimnames = dimnames(gt_dp)
  )

  # Normalization
  snp_prop <- base::sweep(
    gt_dp_num,
    2,
    base::colSums(gt_dp_num, na.rm = TRUE),
    "/"
  ) * 100

  all_samples <- colnames(snp_prop)


  # Window labels aligned with snp_prop rows
  marker_info <- dplyr::select(ttest_data, MarkerName, Chr, Position) %>%
    dplyr::distinct(MarkerName, .keep_all = TRUE)

  chroms    <- marker_info$Chr
  positions <- marker_info$Position

  window_id    <- base::floor((positions - 1) / window_size) + 1
  window_label <- base::paste0(chroms, "_MB", window_id)
  unique_windows <- base::unique(window_label)

  base::stopifnot(length(window_label) == nrow(snp_prop))  # sanity check


  # WINDOW-LEVEL TEST
  window_plot_df <- base::lapply(selected_samples, function(selected_sample) {

    other_samples <- base::setdiff(all_samples, selected_sample)

    # Population mean
    snp_pop_mean <- base::rowMeans(snp_prop[, other_samples], na.rm = TRUE)

    # Aggregate per window
    df <- base::data.frame(
      WINDOW     = unique_windows,
      CHROM      = base::sub("_MB.*", "", unique_windows),
      Megabase   = base::as.integer(base::sub(".*_MB", "", unique_windows)),
      Sample     = selected_sample,
      SampleProp = base::vapply(unique_windows, function(win) {
        base::mean(snp_prop[window_label == win, selected_sample], na.rm = TRUE)
      }, numeric(1)),
      PopAvgProp = base::vapply(unique_windows, function(win) {
        base::mean(snp_pop_mean[window_label == win], na.rm = TRUE)
      }, numeric(1))
    )

    # Unpaired t-test per window
    df$p_raw <- base::vapply(unique_windows, function(win) {
      idx <- window_label == win
      x_sel <- snp_prop[idx, selected_sample]
      x_pop <- snp_pop_mean[idx]
      ok <- base::is.finite(x_sel) & base::is.finite(x_pop)
      if (base::sum(ok) >= 5) stats::t.test(x_sel[ok], x_pop[ok], paired = FALSE)$p.value else NA_real_
    }, numeric(1))

    df$p_adj <- stats::p.adjust(df$p_raw, method = "fdr")
    df$delta <- df$SampleProp - df$PopAvgProp

    # Significance categories
    df$SigDir <- dplyr::case_when(
      base::is.na(df$p_adj) | df$p_adj >= 0.05 ~ "Not significant",
      df$delta < 0 & df$p_adj < 0.001    ~ "Lower: P < 0.001",
      df$delta < 0 & df$p_adj < 0.01     ~ "Lower: P < 0.01",
      df$delta < 0                        ~ "Lower: P < 0.05",
      df$delta > 0 & df$p_adj < 0.001    ~ "Higher: P < 0.001",
      df$delta > 0 & df$p_adj < 0.01     ~ "Higher: P < 0.01",
      TRUE                               ~ "Higher: P < 0.05"
    )

    df
  }) %>% dplyr::bind_rows()


  #plotting

  window_plot_long <- tidyr::pivot_longer(
    window_plot_df,
    cols = c(SampleProp, PopAvgProp),
    names_to = "ID",
    values_to = "DepthProp"
  )

  window_plot_long$LegendGroup <- dplyr::case_when(
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

  ggplot2::ggplot(
    window_plot_long,
    ggplot2::aes(x = Megabase, y = DepthProp, color = LegendGroup)
  ) +
    ggplot2::geom_point(size = 0.75, alpha = 0.85) +
    ggplot2::facet_grid(Sample ~ CHROM, scales = "free_x") +
    ggplot2::scale_color_manual(values = legend_colors, name = "Legend") +
    ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      x = "Window (Mb)",
      y = "Normalized SNP depth (%)",
      title = "Window-level SNP depth: sample vs population"
    )

}


#### example ####
window_level_test_plot_qploidy(qploidy_object = "standardization.tsv.gz",
                               selected_samples = "J432",
                               col2test = "z")
