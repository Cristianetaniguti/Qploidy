# Test suite for hmm_estimate_CN and plot_cn_track
library(testthat)
library(Qploidy2)

test_that("hmm_estimate_CN, plot_cn_track, and other HMM functions work as expected", {
  set.seed(123)
  vcf_file <- tempfile(fileext = ".vcf")
  simulate_vcf(
    seed = 123, file_path = vcf_file,
    n_tetraploid = 5, n_diploid = 2, n_triploid = 2,
    n_markers = 50
  )
  data <- qploidy_read_vcf(vcf_file)
  genos <- qploidy_read_vcf(vcf_file, geno = TRUE)
  genos <- genos[grep("Tetraploid", genos$SampleName), ]
  genos.pos <- qploidy_read_vcf(vcf_file, geno.pos = TRUE)

  stand_file <- tempfile(fileext = ".tsv.gz")
  simu_data_standardized <- standardize(
    data = data,
    genos = genos,
    geno.pos = genos.pos,
    ploidy.standardization = 4,
    threshold.n.clusters = 3,
    n.cores = 1,
    out_filename = stand_file,
    verbose = FALSE
  )
  sample <- unique(simu_data_standardized$data$SampleName)[1]

  # hmm_estimate_CN test
  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = 4
  )
  expect_type(res, "list")
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 6)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = NA
  )

  expect_type(res, "list")
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 6)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = sample,
    chr = 1,
    segment_zscore = FALSE,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    exp_ploidy = NA,
    rm_outliers = FALSE,
  )

  expect_type(res, "list")
  expect_true("by_marker" %in% names(res))
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 5)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  samples <- unique(simu_data_standardized$data$SampleName)
  res <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = samples[1],
    chr = 1,
    cn_grid = c(1:8)
  )

  expect_type(res, "list")
  expect_true("by_marker" %in% names(res))
  expect_true("by_window" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$by_window))
  expect_true(nrow(res$by_window) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$by_window)))
  expect_equal(round(sum(res$by_window$post_CN4),2), 4)
  expect_equal(round(sum(res$by_marker$post_max),2), 50)

  # plot_cn_track test (should return a gg object)
  p <- plot_cn_track(hmm_CN = res,
                     qploidy_standarize_result = simu_data_standardized,
                     sample_id = sample,
                     show_window_lines = TRUE)
  expect_type(p, "list")
  expect_true("gg" %in% class(p$arranged) || "ggarrange" %in% class(p$arranged))

  # Test hmm_estimate_CN_multi
  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = simu_data_standardized,
    sample_ids = unique(simu_data_standardized$data$SampleName),
    n_cores = 1,
    chr = 1,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21
  )

  expect_type(multi_res, "list")
  expect_true(length(multi_res) >= 1)
  expect_true(inherits(multi_res, "hmm_CN"))
  expect_equal(mean(multi_res$by_window$CN_call), 3.1667, tolerance = 1e-3)
  expect_equal(mean(multi_res$by_marker$CN_call), 3.3333, tolerance = 1e-3)

  # Test summarize_cn_mode
  summ <- summarize_cn_mode(df = multi_res, level = "sample")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = simu_data_standardized,
    sample_ids = "all",
    n_cores = 1,
    chr = 1,
    cn_grid = c(1:4)
  )

  expect_type(multi_res, "list")
  expect_true(length(multi_res) >= 1)
  expect_true(inherits(multi_res, "hmm_CN"))
  expect_equal(mean(multi_res$by_window$CN_call), 3.250, tolerance = 1e-3)
  expect_equal(mean(multi_res$by_marker$CN_call), 3.333, tolerance = 1e-3)

  oneCN <- hmm_estimate_CN(
    qploidy_standarize_result = simu_data_standardized,
    sample_id = "Triploid1",
    chr = 1,
    segment_zscore = TRUE,
    snps_per_window =10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    exp_ploidy = 3
  )

  updated_CN <- update_hmm_CN_multi(multi_res, oneCN)

  temp_dir <- tempdir()
  formatted_date <- format(Sys.Date(), "%m_%d_%Y")

  write_hmm_CN(updated_CN, prefix = paste0(temp_dir, "/", formatted_date))

  # Test summarize_cn_mode
  summ <- summarize_cn_mode(df = multi_res, level = "sample")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  summ <- summarize_cn_mode(multi_res, level = "chromosome")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))
  expect_equal(summ$CN_mode, c(2,2,4,4,4,4,4,3,3))

  p <- compare_cn_track(multi_res, samples_to_plot = samples)

  expect_true("gg" %in% class(p))

})
