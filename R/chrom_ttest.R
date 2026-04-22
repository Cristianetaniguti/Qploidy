#' Chromosome-level Depth Comparison Plot
#'
#' Generates a bar plot comparing the normalized SNP depth of selected samples against the population mean at the chromosome level.
#' Includes significance stars based on paired t-tests per chromosome.
#'
#' @param input Accepts one of:
#'   1. A data.frame with columns: MarkerName, SampleName, Chr, Position, R. Column with z is optional but required if col2test = "z".
#'   2. An object of class 'qploidy_standardization' (from Qploidy standardization functions)
#'   3. A character string: path to a .tsv[.gz] file with standardized Qploidy data
#'   The function will validate and process the input according to its type.
#'
#' @param geno.pos Optional. A data.frame with columns: MarkerName, Chromosome, Position. Required if input data.frame does not contain Chr and Position columns.
#' @param selected_samples Character vector of sample names to compare against the population.
#' @param col2test Column name to use for depth comparison. Options are "R" or "z".
#'
#' @details
#' The function reads standardized SNP depth data from the Qploidy object, normalizes it per sample,
#' calculates per-chromosome averages for the selected sample and the population mean, and performs
#' paired t-tests to assess differences. Significance is annotated using stars (*, **, ***), with
#' Bonferroni-adjusted p-values.
#'
#' The resulting ggplot2 bar plot displays per-chromosome SNP depth for the sample vs population,
#' faceted by sample, with normalized percentages on the y-axis and significance indicated above bars.
#'
#' @return A ggplot2 object representing the chromosome-level SNP depth comparison.
#'
#' @examples
#'\dontrun{
#' chromosome_level_test_plot_qploidy(
#'   input = "standardization.tsv.gz",
#'   selected_samples = "J432",
#'   col2test = "z"
#' )
#' }
#'
#' @importFrom dplyr select distinct bind_rows all_of
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom tibble column_to_rownames
#' @importFrom ggplot2 ggplot aes geom_bar geom_text scale_fill_manual scale_color_manual guides facet_wrap
#' @importFrom stats t.test p.adjust
#' @importFrom scales percent_format
#'
#' @author Josue Chinchilla-Vargas
#'
#' @export
chromosome_level_test_plot_qploidy <- function(
    input,
    geno.pos = NULL,
    selected_samples,
    col2test = c("R","z")
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

  # Marker info & chromosomes
  marker_info <- select(ttest_data, MarkerName, Chr, Position) %>%
    distinct(MarkerName, .keep_all = TRUE)

  chroms <- marker_info$Chr
  chrom_levels <- unique(chroms)

  # Significance helper
  sig_label <- function(pvals) {
    case_when(
      is.na(pvals)  ~ "",
      pvals < 0.001 ~ "***",
      pvals < 0.01  ~ "**",
      pvals < 0.05  ~ "*",
      TRUE          ~ ""
    )
  }

  # Chromosome-level test
  chr_plot_df <- lapply(selected_samples, function(selected_sample) {

    other_samples <- setdiff(all_samples, selected_sample)

    # Population mean
    snp_pop_mean <- rowMeans(snp_prop[, other_samples], na.rm = TRUE)

    # Aggregate per chromosome
    df <- data.frame(
      CHROM      = chrom_levels,
      Sample     = selected_sample,
      SampleProp = vapply(chrom_levels, function(chr) {
        mean(snp_prop[chroms == chr, selected_sample], na.rm = TRUE)
      }, numeric(1)),
      PopAvgProp = vapply(chrom_levels, function(chr) {
        mean(snp_pop_mean[chroms == chr], na.rm = TRUE)
      }, numeric(1))
    )

    # Paired t-test per chromosome
    df$p_raw <- vapply(chrom_levels, function(chr) {
      idx <- chroms == chr
      x_sel <- snp_prop[idx, selected_sample]
      x_pop <- snp_pop_mean[idx]
      ok <- is.finite(x_sel) & is.finite(x_pop)
      if (sum(ok) >= 20)
        t.test(x_sel[ok], x_pop[ok], paired = TRUE)$p.value
      else NA_real_
    }, numeric(1))

    df$p_adj <- p.adjust(df$p_raw, method = "bonferroni")
    df$delta <- df$SampleProp - df$PopAvgProp
    df$Label <- sig_label(df$p_adj)

    df
  }) %>% bind_rows()

  # Reshape for plotting
  chr_plot_long <- pivot_longer(
    chr_plot_df,
    cols = c(SampleProp, PopAvgProp),
    names_to = "ID",
    values_to = "DepthProp"
  )

  # Significance factor (for legend)
  chr_plot_df$SignifCat <- factor(
    ifelse(chr_plot_df$Label == "", "Not significant", chr_plot_df$Label),
    levels = c("Not significant", "*", "**", "***"),
    labels = c("Not significant", "P < 0.05", "P < 0.01", "P < 0.001")
  )

  # Plotting
  ggplot(chr_plot_long, aes(x = CHROM, y = DepthProp, fill = ID)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +

    # Significance stars
    geom_text(
      data = chr_plot_df,
      aes(
        x = CHROM,
        y = pmax(SampleProp, PopAvgProp) * 1.05,
        label = Label
      ),
      inherit.aes = FALSE,
      vjust = 0,
      size = 3.5
    ) +

    # Invisible layer to include significance in legend
    geom_text(
      data = chr_plot_df,
      aes(
        x = CHROM,
        y = 0,
        label = Label,
        color = SignifCat
      ),
      inherit.aes = FALSE,
      show.legend = TRUE,
      size = 0
    ) +

    scale_fill_manual(
      values = c(
        "SampleProp" = "gold",
        "PopAvgProp" = "grey70"
      )
    ) +

    scale_color_manual(
      values = c(
        "Not significant" = "black",
        "P < 0.05"        = "black",
        "P < 0.01"        = "black",
        "P < 0.001"       = "black"
      ),
      drop = FALSE
    ) +

    guides(
      fill = guide_legend(order = 1),
      color = guide_legend(
        order = 2,
        title = "Significance",
        override.aes = list(
          label = c("", "*", "**", "***"),
          size  = 6
        )
      )
    ) +

    facet_wrap(~ Sample, scales = "free_x") +
    scale_y_continuous(labels = percent_format(scale = 1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      x = "Chromosome",
      y = "Normalized SNP depth (%)",
      title = "Chromosome-level SNP depth: sample vs population"
    )
}
