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

  result1 <- area_estimate_ploidy(
    qploidy_standardization = qploidy_standardization,
    samples = "all",
    level = "chromosome",
    ploidies = c(2, 3, 4),
    area = 0.75
  )

  expect_s3_class(result1, "qploidy_area_ploidy_estimation")
  expect_true(!is.null(result1$ploidy))

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


  multi_esti <- hmm_estimate_CN_multi(qploidy_standarize_result = qploidy_standardization,
                                      sample_ids = "all",
                                      n_cores = 1 ,
                                      chr = NULL,
                                      snps_per_window = 20,
                                      cn_grid =  c(2,3,4) ,
                                      M = 100 ,
                                      bw = 0.03,
                                      exp_ploidy = 2 ,
                                      het_lims = c(0,1),
                                      het_quantile = 0.8,
                                      baf_weight = 1 ,
                                      z_range =  0.2 ,
                                      transition_jump =  0.995  ,
                                      max_iter = 60,
                                      z_only = FALSE)


  summ_hmm <- summarize_cn_mode(multi_esti, level = "chromosome")
  merged_table <- merge_cn_summary_with_estimates(summ_hmm, result1, level="chromosome")
  ploidies <- list(area = result1, hmm = merged_table)

  all_hmm <- TRUE
  if(all_hmm){
    tab <- ploidies$hmm
    num_cols <- which(apply(tab, 2, is.numeric))
    if(length(num_cols) != 0)
      tab[, num_cols] <- apply(tab[, num_cols], 2, function(x) round(x, 4))
    res <- datatable(tab,
                     selection = "single",
                     options = list(pageLength = 10, scrollX = TRUE))
  } else {
    ploidy <- ploidies$area$ploidy
    ploidy[which(ploidies$area$diff_first_second < 0.01)] <- NA

    if(ncol(ploidy) > 1){
      num_cols <- which(apply(ploidy, 2, is.numeric))
      if(length(num_cols) != 0)
        ploidy[, num_cols] <- apply(ploidy[, num_cols], 2, function(x) round(x, 4))
      colnames(ploidy) <- colnames(ploidies$area$ploidy)

    } else {
      ploidy <- as.data.frame(round(ploidy, 4))
      rownames(ploidy) <- rownames(ploidies$area$ploidy)
      colnames(ploidy) <- colnames(ploidies$area$ploidy)
    }

    res <- datatable(ploidy,
                     selection = "single",
                     options = list(pageLength = 10, scrollX = TRUE))
  }

  expect_equal(mean(tab[,3]), 4)
  expect_equal(mean(tab[,6]), 0.1818, tolerance = 0.001)

})
