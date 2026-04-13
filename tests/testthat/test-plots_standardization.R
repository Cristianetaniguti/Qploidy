library(testthat)
library(ggplot2)

# ---------------------------------------------------------------------------
# Shared fixtures
# ---------------------------------------------------------------------------

make_std <- function(n_markers = 100, n_samples = 50, ploidy = 2, seed = 42) {
  fake <- simulate_standardization_input(
    n_markers = n_markers,
    n_samples  = n_samples,
    ploidy     = ploidy,
    seed       = seed
  )
  standardize(
    data                    = fake$sample_data,
    genos                   = fake$geno_data,
    geno.pos                = fake$geno_pos,
    threshold.missing.geno  = 1,
    threshold.geno.prob     = 0.5,
    ploidy.standardization  = ploidy,
    threshold.n.clusters    = 2,
    n.cores                 = 1,
    type                    = "counts",
    verbose                 = FALSE
  )
}

make_baf_sample <- function(std, sample_id = "S1") {
  data <- std$data[std$data$SampleName == sample_id, ]
  baf  <- tidyr::pivot_wider(
    data[, c("MarkerName", "SampleName", "Chr", "Position", "baf", "z", "ratio", "R")],
    names_from  = SampleName,
    values_from = baf
  )
  colnames(baf)[ncol(baf)] <- "sample"
  baf
}

# ---------------------------------------------------------------------------
# plot_baf
# ---------------------------------------------------------------------------

test_that("plot_baf returns a ggplot with minimal arguments", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf(bsmp, area_single = 0.75, ploidy = 2)

  expect_s3_class(p, "gg")
})

test_that("plot_baf works with add_expected_peaks and colors", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf(bsmp, area_single = 0.75, ploidy = 2,
                add_expected_peaks = TRUE, colors = TRUE)

  expect_s3_class(p, "gg")
})

test_that("plot_baf works with centromere string input", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf(bsmp, area_single = 0.75, ploidy = 2,
                add_centromeres = TRUE,
                centromeres = c("1" = 500000L))

  expect_s3_class(p, "gg")
})

test_that("plot_baf stops with mismatched ploidy vector length", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  expect_error(
    plot_baf(bsmp, area_single = 0.75, ploidy = c(2, 4), add_expected_peaks = TRUE),
    "Provide a ploidy for each chromosome"
  )
})

# ---------------------------------------------------------------------------
# plot_baf_hist
# ---------------------------------------------------------------------------

test_that("plot_baf_hist returns a ggplot", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf_hist(bsmp, area_single = 0.75, ploidy = 2)

  expect_s3_class(p, "gg")
})

test_that("plot_baf_hist works with BAF_hist_overall = TRUE", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf_hist(bsmp, area_single = 0.75, ploidy = 2,
                     BAF_hist_overall = TRUE)

  expect_s3_class(p, "gg")
})

test_that("plot_baf_hist works with ratio = TRUE", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf_hist(bsmp, area_single = 0.75, ploidy = 2, ratio = TRUE, add_estimated_peaks = FALSE)

  expect_s3_class(p, "gg")
})

test_that("plot_baf_hist works with add_estimated_peaks and add_expected_peaks", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf_hist(bsmp, area_single = 0.75, ploidy = 2,
                     add_estimated_peaks = TRUE, add_expected_peaks = TRUE)

  expect_s3_class(p, "gg")
})

test_that("plot_baf_hist works with rm_homozygous = TRUE", {
  std  <- make_std()
  bsmp <- make_baf_sample(std)

  p <- plot_baf_hist(bsmp, area_single = 0.75, ploidy = 2,
                     rm_homozygous = TRUE)

  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# plot_qploidy_standardization
# ---------------------------------------------------------------------------

test_that("plot_qploidy_standardization stops on non-standardization object", {
  expect_error(
    plot_qploidy_standardization(list(), sample = "S1"),
    "qploidy_standardization"
  )
})

test_that("plot_qploidy_standardization stops when sample is NULL", {
  std <- make_std()
  expect_error(
    plot_qploidy_standardization(std),
    "Define sample ID"
  )
})

test_that("plot_qploidy_standardization stops on unknown chromosome", {
  std <- make_std()
  expect_error(
    plot_qploidy_standardization(std, sample = "S1", chr = "99"),
    "Chromosome not found"
  )
})

test_that("plot_qploidy_standardization returns a figure with type='all'", {
  std <- make_std()
  p   <- plot_qploidy_standardization(std, sample = "S1", type = "all", ploidy = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_qploidy_standardization works for each canonical type", {
  std   <- make_std()
  types <- c("BAF", "zscore", "BAF_hist", "Ratio_hist",
             "BAF_hist_overall", "Ratio_hist_overall", "ratio", "R", "het")

  for (tp in types) {
    p <- plot_qploidy_standardization(std, sample = "S1", type = tp, ploidy = 2)
  }
})

test_that("plot_qploidy_standardization accepts type aliases", {
  std <- make_std()

  aliases <- list(
    het              = c("heterozygosity", "heterozygous"),
    BAF              = c("baf"),
    zscore           = c("z", "z-score"),
    BAF_hist         = c("baf_hist", "baf_histogram"),
    Ratio_hist       = c("ratio_hist", "ratio_histogram", "Ratio_histogram"),
    BAF_hist_overall = c("baf_hist_overall", "baf_hist_sample", "baf_histogram_overall"),
    Ratio_hist_overall = c("ratio_hist_overall", "ratio_histogram_overall")
  )

  for (alias_list in aliases) {
    for (alias in alias_list) {
      p <- plot_qploidy_standardization(std, sample = "S1", type = alias, ploidy = 2)
      expect_s3_class(p, "gg")
    }
  }
})

test_that("plot_qploidy_standardization warns on unknown type and still plots known ones", {
  std <- make_std()
  expect_warning(
    p <- plot_qploidy_standardization(std, sample = "S1",
                                      type = c("BAF", "not_a_type"), ploidy = 2),
    "Unknown plot type"
  )
  expect_s3_class(p, "gg")
})

test_that("plot_qploidy_standardization works when chr is specified by index", {
  std <- make_std()
  p   <- plot_qploidy_standardization(std, sample = "S1", type = "BAF",
                                      chr = 1, ploidy = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_qploidy_standardization works with centromeres", {
  std  <- make_std()
  chrs <- unique(std$data$Chr)
  cent <- setNames(rep(500000L, length(chrs)), chrs)

  p <- plot_qploidy_standardization(
    std, sample = "S1", type = "BAF", ploidy = 2,
    add_centromeres = TRUE, centromeres = cent
  )
  expect_s3_class(p, "gg")
})

# ---------------------------------------------------------------------------
# all_resolutions_plots
# ---------------------------------------------------------------------------

test_that("all_resolutions_plots returns a named list of plots", {
  std  <- make_std()
  chrs <- unique(std$data$Chr)
  cent <- setNames(rep(500000L, length(chrs)), chrs)

  res <- all_resolutions_plots(
    data_standardized = std,
    sample            = "S1",
    ploidy            = 2,
    centromeres       = cent
  )

  expect_type(res, "list")
  expect_named(res, c("chromosome", "chromosome_arm", "sample"))
  expect_s3_class(res$chromosome, "gg")
  expect_s3_class(res$chromosome_arm, "gg")
  expect_s3_class(res$sample, "gg")
})

test_that("all_resolutions_plots works without centromeres", {
  std <- make_std()

  res <- all_resolutions_plots(
    data_standardized = std,
    sample            = "S1",
    ploidy            = 2,
    centromeres       = NULL
  )

  expect_type(res, "list")
  expect_s3_class(res$chromosome, "gg")
  expect_null(res$chromosome_arm)
  expect_s3_class(res$sample, "gg")
})

test_that("all_resolutions_plots stops on non-standardization object", {
  expect_error(
    all_resolutions_plots(list(), sample = "S1", ploidy = 2, centromeres = NULL),
    "qploidy_standardization"
  )
})

test_that("all_resolutions_plots saves PNG files when file_name is provided", {
  std  <- make_std()
  tmp  <- tempfile()

  res <- all_resolutions_plots(
    data_standardized = std,
    sample            = "S1",
    ploidy            = 2,
    centromeres       = NULL,
    file_name         = tmp
  )

  expect_true(file.exists(paste0(tmp, "_res:chromosome.png")))
  expect_true(file.exists(paste0(tmp, "_res:sample.png")))
})

# ---------------------------------------------------------------------------
# plot_xy_with_ploidy_guides
# ---------------------------------------------------------------------------

make_xy_df <- function(n = 50, ploidy = 2, seed = 1) {
  set.seed(seed)
  dosages <- sample(0:ploidy, n, replace = TRUE)
  data.frame(
    X          = pmax(0, rnorm(n, mean = 200 - dosages * 80, sd = 15)),
    Y          = pmax(0, rnorm(n, mean = 20  + dosages * 80, sd = 15)),
    SampleName = sample(c("A", "B"), n, replace = TRUE),
    geno       = dosages
  )
}

test_that("plot_xy_with_ploidy_guides returns a ggplot", {
  df <- make_xy_df()
  p  <- plot_xy_with_ploidy_guides(df, ploidy = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_xy_with_ploidy_guides works with sample='all'", {
  df <- make_xy_df()
  p  <- plot_xy_with_ploidy_guides(df, ploidy = 2, sample = "all")
  expect_s3_class(p, "gg")
})

test_that("plot_xy_with_ploidy_guides works highlighting a single sample", {
  df <- make_xy_df()
  p  <- plot_xy_with_ploidy_guides(df, ploidy = 2, sample = "A")
  expect_s3_class(p, "gg")
})

test_that("plot_xy_with_ploidy_guides works with color_by_geno = TRUE", {
  df <- make_xy_df()
  p  <- plot_xy_with_ploidy_guides(df, ploidy = 2, color_by_geno = TRUE)
  expect_s3_class(p, "gg")
})

test_that("plot_xy_with_ploidy_guides warns when highlighted sample has no points", {
  df <- make_xy_df()
  expect_warning(
    plot_xy_with_ploidy_guides(df, ploidy = 2, sample = "MISSING"),
    "no plotted points"
  )
})

test_that("plot_xy_with_ploidy_guides stops when X/Y columns missing", {
  df <- data.frame(A = 1:3, B = 1:3)
  expect_error(plot_xy_with_ploidy_guides(df, ploidy = 2))
})

test_that("plot_xy_with_ploidy_guides stops when all points are NA", {
  df <- data.frame(X = NA_real_, Y = NA_real_)
  expect_error(plot_xy_with_ploidy_guides(df, ploidy = 2), "No points to plot")
})

# ---------------------------------------------------------------------------
# plot_baf_with_ploidy_guides
# ---------------------------------------------------------------------------

make_baf_df <- function(n = 60, ploidy = 2, seed = 7) {
  set.seed(seed)
  dosages <- sample(0:ploidy, n, replace = TRUE)
  data.frame(
    baf        = pmin(1, pmax(0, rnorm(n, mean = dosages / ploidy, sd = 0.05))),
    R          = rpois(n, lambda = 200),
    SampleName = sample(c("A", "B"), n, replace = TRUE),
    ratio      = pmin(1, pmax(0, rnorm(n, mean = dosages / ploidy, sd = 0.05)))
  )
}

test_that("plot_baf_with_ploidy_guides returns a ggplot", {
  df <- make_baf_df()
  p  <- plot_baf_with_ploidy_guides(df, ploidy = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides works with sample='all'", {
  df <- make_baf_df()
  p  <- plot_baf_with_ploidy_guides(df, ploidy = 2, sample = "all")
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides works highlighting a single sample", {
  df <- make_baf_df()
  p  <- plot_baf_with_ploidy_guides(df, ploidy = 2, sample = "A")
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides works with normalize_depth = FALSE", {
  df <- make_baf_df()
  p  <- plot_baf_with_ploidy_guides(df, ploidy = 2, normalize_depth = FALSE)
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides works with fallback_to_ratio", {
  df      <- make_baf_df()
  df$baf[1:10] <- NA
  p <- plot_baf_with_ploidy_guides(df, ploidy = 2, fallback_to_ratio = TRUE)
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides works with custom radius", {
  df <- make_baf_df()
  p  <- plot_baf_with_ploidy_guides(df, ploidy = 2,
                                    normalize_depth = TRUE, radius = 150)
  expect_s3_class(p, "gg")
})

test_that("plot_baf_with_ploidy_guides stops when baf column is missing", {
  df <- data.frame(R = 1:3, SampleName = "A")
  expect_error(plot_baf_with_ploidy_guides(df, ploidy = 2), "'baf' not found")
})

test_that("plot_baf_with_ploidy_guides stops when R column is missing", {
  df <- data.frame(baf = c(0.1, 0.5), SampleName = "A")
  expect_error(plot_baf_with_ploidy_guides(df, ploidy = 2))
})

# ---------------------------------------------------------------------------
# plot_geno_by_marker
# ---------------------------------------------------------------------------

make_geno_df <- function(ploidy = 2, n_markers = 3, n_per_geno = 10, seed = 99) {
  set.seed(seed)
  markers <- paste0("M", seq_len(n_markers))
  do.call(rbind, lapply(markers, function(mk) {
    dosages <- rep(0:ploidy, each = n_per_geno)
    n       <- length(dosages)
    data.frame(
      MarkerName = mk,
      geno       = dosages,
      ratio      = pmin(1, pmax(0, rnorm(n, mean = dosages / ploidy, sd = 0.05))),
      baf        = pmin(1, pmax(0, rnorm(n, mean = dosages / ploidy, sd = 0.05)))
    )
  }))
}

test_that("plot_geno_by_marker returns a ggplot", {
  df <- make_geno_df()
  p  <- plot_geno_by_marker(df, ploidy = 2)
  expect_s3_class(p, "gg")
})

test_that("plot_geno_by_marker infers ploidy from data when NULL", {
  df <- make_geno_df(ploidy = 4)
  p  <- plot_geno_by_marker(df, ploidy = NULL)
  expect_s3_class(p, "gg")
})

test_that("plot_geno_by_marker handles NA genotypes without error", {
  df          <- make_geno_df()
  df$geno[1:5] <- NA
  p <- plot_geno_by_marker(df, ploidy = 2)
  expect_s3_class(p, "gg")
})
