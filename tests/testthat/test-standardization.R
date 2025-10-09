library(dplyr)
library(testthat)

test_that("get_centers returns expected result", {
  fake_input <- simulate_standardization_input(n_markers = 10, n_samples = 5, ploidy = 2, seed = 2025)

  # Example that marker has all 3 possible dosages
  m <- fake_input$standardization_input %>% filter(MarkerName == "m2")
  result <- get_centers(m, ploidy = 2)

  expect_type(result, "list")
  expect_named(result, c("rm", "centers_theta", "MarkerName", "n.clusters"))
  expect_true(length(result$centers_theta) == 3)
  expect_equal(result$centers_theta, c(0.03685, 0.47120, 0.97875), tolerance = 1e-4)

  # Example that marker has only 2 of the 3 possible dosages - 1 imputed
  m1 <- fake_input$standardization_input %>% filter(MarkerName == "m1")
  result <- get_centers(m1, ploidy = 2) # marker get flag 2 - will be removed
  expect_true(result$rm == 2)
  result <- get_centers(m1, ploidy = 2, n.clusters.thr = 2)


  expect_type(result, "list")
  expect_named(result, c("rm", "centers_theta", "MarkerName", "n.clusters"))
  expect_true(length(result$centers_theta) == 3)
  expect_equal(result$centers_theta, c(0.086750, 0.520875, 0.955000), tolerance = 1e-4)
})

test_that("get_baf works as expected", {
  centers_theta <- c(0.1, 0.5, 0.9)
  theta <- c(0.05, 0.3, 0.6, 0.95)

  bafs <- get_baf(theta, centers_theta, ploidy = 2)

  expect_equal(length(bafs), length(theta))
  expect_true(all(bafs >= 0 & bafs <= 1))
  expect_equal(sum(bafs), 1.875, tolerance = 1e-4)
})

test_that("get_baf_par returns correct BAFs for 5 markers and 4 samples", {
  set.seed(123)

  n_markers <- 5
  n_samples <- 4
  marker_ids <- paste0("m", 1:n_markers)
  sample_ids <- paste0("S", 1:n_samples)

  # Simulate theta values between 0 and 1 (e.g., sample allelic ratios)
  theta_matrix <- matrix(runif(n_markers * n_samples, min = 0.05, max = 0.95), nrow = n_markers)
  theta_df <- data.frame(ID = marker_ids, theta_matrix)
  colnames(theta_df)[-1] <- sample_ids

  # Simulate cluster centers (3 centers for diploid: 0, 0.5, 1)
  # We'll make them a bit variable to simulate real centroids
  centers_matrix <- t(replicate(n_markers, sort(c(
    runif(1, 0.05, 0.15),
    runif(1, 0.45, 0.55),
    runif(1, 0.85, 0.95)
  ))))
  centers_df <- data.frame(ID = marker_ids, centers_matrix)
  colnames(centers_df)[-1] <- 1:3

  # Prepare input list
  par_input <- list(theta_df, centers_df)

  # Call the function
  result <- get_baf_par(par_input, ploidy = 2)

  # Assertions
  expect_type(result, "list")
  expect_length(result, n_markers)
  expect_equal(sum(result[[1]]), 2.1509, tolerance = 0.001)
})

test_that("get_zscore returns z-score table", {
  fake_input <- simulate_standardization_input(n_markers = 10, n_samples = 5, ploidy = 2, seed = 2025)

  result <- get_zscore(data = fake_input$sample_data, geno.pos = fake_input$geno_pos)

  expect_s3_class(result, "data.frame")
  expect_true("z" %in% colnames(result))
  expect_equal(nrow(result), nrow(fake_input$sample_data))
  expect_equal(var(result$z), 0.8163, tolerance = 0.001)
})

test_that("standardize runs and returns expected object", {
  fake_input <- simulate_standardization_input(n_markers = 10, n_samples = 5, ploidy = 2, seed = 2025)

  result <- standardize(
    data = fake_input$sample_data,
    genos = fake_input$geno_data,
    geno.pos = fake_input$geno_pos,
    threshold.missing.geno = 1,
    threshold.geno.prob = 0.5,
    ploidy.standardization = 2,
    threshold.n.clusters = 2,
    n.cores = 1,
    type = "intensities",
    verbose = FALSE
  )

  expect_s3_class(result, "qploidy_standardization")
  expect_true("data" %in% names(result))
  expect_equal(sum(result$data$baf, na.rm = TRUE), 27.3644 ,tolerance = 0.001)


  temp <- tempfile(fileext = ".tsv")
  result <- standardize(
    data = fake_input$sample_data,
    genos = fake_input$geno_data,
    geno.pos = fake_input$geno_pos,
    threshold.missing.geno = 1,
    threshold.geno.prob = 0.5,
    ploidy.standardization = 2,
    threshold.n.clusters = 2,
    n.cores = 1,
    type = "counts",
    verbose = TRUE,
    out_filename = temp
  )

  expect_s3_class(result, "qploidy_standardization")
  expect_true("data" %in% names(result))
  expect_equal(sum(result$data$baf, na.rm = TRUE), 27.3644 ,tolerance = 0.001)

  # Check print.qploidy_standardization
  txt <- print(result)
  expect_equal(dim(txt), c(6,3))
  expect_equal(as.character(txt[4,3]), "(0 %)  ")

  # Check read.qploidy_standardization
  result <- read_qploidy_standardization(temp)

  expect_s3_class(result, "qploidy_standardization")
  expect_true("data" %in% names(result))
  expect_equal(sum(result$data$baf, na.rm = TRUE), 27.3644 ,tolerance = 0.001)

  # Check print.qploidy_standardization
  txt <- print(result)
  expect_equal(dim(txt), c(6,3))
  expect_equal(as.character(txt[4,3]), "(0 %)  ")
})

test_that("rm_outlier removes extreme values", {
  # Create a dataset with 19 'normal' values and 1 extreme outlier
  normal_theta <- rnorm(19, mean = 0.2, sd = 0.01)
  outlier_theta <- 3.0
  test_theta <- c(normal_theta, outlier_theta)

  test_df <- data.frame(theta = test_theta)
  cleaned <- rm_outlier(test_df)

  expect_true(length(cleaned) < nrow(test_df))
})


test_that("updog_centers returns correct output structure", {
  skip_if_not_installed("updog")

  library(dplyr)
  library(updog)

  set.seed(123)
  n_markers <- 5
  n_individuals <- 10
  ploidy <- 2

  # Simulate varying total read depth (size) for each marker × individual
  sizemat <- matrix(rpois(n_markers * n_individuals, lambda = 100), nrow = n_markers, ncol = n_individuals)
  sizemat[sizemat == 0] <- 1  # avoid 0 depth

  # Simulate reference reads with moderate allele bias
  refmat <- matrix(rbinom(n_markers * n_individuals, size = c(sizemat), prob = 0.7),
                   nrow = n_markers, ncol = n_individuals)
  altmat <- sizemat - refmat

  # Name dimensions
  marker_ids <- paste0("m", 1:n_markers)
  sample_ids <- paste0("ind", 1:n_individuals)

  rownames(refmat) <- rownames(altmat) <- rownames(sizemat) <- marker_ids
  colnames(refmat) <- colnames(altmat) <- colnames(sizemat) <- sample_ids

  # Create multidog object
  md <- multidog(
    refmat = refmat,
    sizemat = sizemat,
    ploidy = ploidy,
    model = "norm"
  )

  # No markers removed
  rm.mks <- character(0)

  result <- updog_centers(md, threshold.n.clusters = 2, rm.mks = rm.mks)

  # Basic checks
  expect_type(result, "list")
  expect_length(result, 5)
  expect_equal(sum(result$m1$centers_theta), 1.502, tolerance = 0.001)
})

