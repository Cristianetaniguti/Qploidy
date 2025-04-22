test_that("area_estimate_ploidy handles missing data correctly", {
  fake_input <- simulate_standardization_input(n_markers = 100, n_samples = 50, ploidy = 4, seed = 2025)

  qploidy_standardization <- standardize(
    data = fake_input$sample_data,
    genos = fake_input$geno_data,
    geno.pos = fake_input$geno_pos,
    threshold.missing.geno = 1,
    threshold.geno.prob = 0.5,
    ploidy.standardization = 4,
    threshold.n.clusters = 4,
    n.cores = 1,
    type = "intensities",
    verbose = FALSE
  )

  result <- area_estimate_ploidy(
    qploidy_standardization = qploidy_standardization,
    samples = "all",
    level = "chromosome",
    ploidies = c(2, 3, 4),
    area = 0.75
  )

  expect_s3_class(result, "qploidy_area_ploidy_estimation")
  expect_true(!is.null(result$ploidy))

  result <- area_estimate_ploidy(
    qploidy_standardization = qploidy_standardization,
    samples = "all",
    level = "chromosome-arm",
    centromeres = c("1" = 500000, "2" = 500000, "3" = 500000),
    ploidies = c(2, 4, 5),
    area = 0.75
  )

  expect_s3_class(result, "qploidy_area_ploidy_estimation")
  expect_true(!is.null(result$ploidy))

  result_merged <- merge_arms_format(result, filter_diff = 0.001)
  result_merged$ploidy

  expect_s3_class(result_merged, "qploidy_area_ploidy_estimation")
  expect_true(!is.null(result_merged$ploidy))
  expect_equal(ncol(result_merged$ploidy), 3)
  expect_equal(nrow(result_merged$ploidy), 50)
  expect_equal(result_merged$ploidy[1, 1], "4")
  expect_equal(result_merged$ploidy[2, 1], "4")

  aneu <- get_aneuploids(result_merged$ploidy)
  expect_equal(all(aneu), FALSE)

  # Capture the printed output
  output <- capture.output(print(result_merged))

  # Check that the output contains expected summary information
  expect_true(any(grepl("Object of class qploidy_area_ploidy_estimation", output)))
  expect_true(any(grepl("Number of samples:", output)))
  expect_true(any(grepl("Chromosomes:", output)))
  expect_true(any(grepl("Tested ploidies:", output)))
  expect_true(any(grepl("Number of euploid samples:", output)))
  expect_true(any(grepl("Number of potential aneuploid samples:", output)))
  expect_true(any(grepl("Number of highly inbred samples:", output)))


  temp2 <- tempfile()
  p_list <- all_resolutions_plots(
    data_standardized = qploidy_standardization,
    sample = "S1",
    types_chromosome = c("BAF", "BAF_hist", "zscore", "het", "ratio"),
    ploidy = 4,
    centromeres = c("1" = 500000, "2" = 500000, "3" = 500000),
    file_name = temp2
  )

  expect_true(all(c("chromosome", "chromosome_arm") %in% names(p_list)))
  expect_true(file.exists(paste0(temp2, "_res:chromosome_arm.png")))
  expect_true(file.exists(paste0(temp2, "_res:chromosome.png")))
  expect_true(file.exists(paste0(temp2, "_res:sample.png")))
  
})
