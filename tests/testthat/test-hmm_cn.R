# Test suite for hmm_estimate_CN and plot_cn_track
library(testthat)
library(Qploidy)

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
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    bw = 0.03,
    exp_ploidy = 4
  )
  expect_type(res, "list")
  expect_true("result" %in% names(res))
  expect_true("params" %in% names(res))
  expect_s3_class(res, "hmm_CN")
  expect_true(is.data.frame(res$result))
  expect_true(nrow(res$result) > 0)
  expect_true(all(c("Sample", "Chr", "WindowID", "CN_call") %in% names(res$result)))

  # plot_cn_track test (should return a gg object)
  p <- plot_cn_track(hmm_CN = res, qploidy_standarize_result = simu_data_standardized, sample_id = sample, show_window_lines = TRUE)
  expect_true("gg" %in% class(p) || "ggarrange" %in% class(p))

  # Test hmm_estimate_CN_multi
  multi_res <- hmm_estimate_CN_multi(
    qploidy_standarize_result = simu_data_standardized,
    sample_ids = unique(simu_data_standardized$data$SampleName),
    n_cores = 1,
    chr = 1,
    snps_per_window = 10,
    min_snps_per_window = 5,
    cn_grid = c(2, 3, 4),
    M = 21,
    bw = 0.03,
    exp_ploidy = 4
  )
  expect_type(multi_res, "list")
  expect_true(length(multi_res) >= 1)
  expect_true(inherits(multi_res, "hmm_CN"))

  # Test summarize_cn_mode
  summ <- summarize_cn_mode(res, level = "sample")
  expect_true(is.data.frame(summ))
  expect_true("CN_mode" %in% names(summ))

  # Test merge_cn_summary_with_estimates
  dummy_area <- area_estimate_ploidy(simu_data_standardized, samples = "all", level = "sample", ploidies = c(2,3,4))

  merged <- merge_cn_summary_with_estimates(summ, dummy_area, level = "sample")
  expect_true(is.data.frame(merged))
  expect_true("CN_HMM" %in% names(merged))
})
