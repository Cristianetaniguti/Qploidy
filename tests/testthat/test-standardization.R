context("standardization")

library(dplyr)

# Setup
set.seed(2025)

n_markers <- 10
n_samples <- 5
marker_ids <- paste0("m", 1:n_markers)
sample_ids <- paste0("S", 1:n_samples)

# Expand marker x sample
data_grid <- expand.grid(MarkerName = marker_ids, SampleName = sample_ids)

# Simulate genotype and values
final_dataset <- data_grid %>%
  rowwise() %>%
  mutate(
    geno = sample(0:2, 1),
    X = case_when(
      geno == 0 ~ round(rnorm(1, mean = 280, sd = 10)),
      geno == 1 ~ round(rnorm(1, mean = 150, sd = 10)),
      geno == 2 ~ round(rnorm(1, mean = 20,  sd = 5))
    ),
    Y = case_when(
      geno == 0 ~ round(rnorm(1, mean = 20,  sd = 5)),
      geno == 1 ~ round(rnorm(1, mean = 150, sd = 10)),
      geno == 2 ~ round(rnorm(1, mean = 280, sd = 10))
    ),
    R = X + Y,
    ratio = round(Y / R, 4),
    prob = round(runif(1, 0.85, 1.00), 2)
  ) %>%
  ungroup()

# Create datasets
sample_data <- final_dataset %>%
  select(MarkerName, SampleName, X, Y, R, ratio)

geno_data <- final_dataset %>%
  select(MarkerName, SampleName, geno, prob)

geno_pos <- tibble::tibble(
  MarkerName = marker_ids,
  Chromosome = rep(c("1", "1", "1", "2", "2"), length.out = n_markers),
  Position = sample(1e5:1e6, n_markers)
)

standardization_input <- inner_join(
  sample_data %>% select(MarkerName, SampleName, ratio) %>% rename(theta = ratio),
  geno_data %>% select(MarkerName, SampleName, geno),
  by = c("MarkerName", "SampleName")
)

test_that("get_centers returns expected result", {
  # Example that marker has all 3 possible dosages
  m3 <- standardization_input %>% filter(MarkerName == "m3")
  result <- get_centers(m3, ploidy = 2)

  expect_type(result, "list")
  expect_named(result, c("rm", "centers_theta", "MarkerName", "n.clusters"))
  expect_true(length(result$centers_theta) == 3)
  expect_equal(result$centers_theta, c(0.06485, 0.49305, 0.92860), tolerance = 1e-4)

  # Example that marker has only 2 of the 3 possible dosages - 1 imputed
  m1 <- standardization_input %>% filter(MarkerName == "m1")
  result <- get_centers(m1, ploidy = 2) # marker get flag 2 - will be removed
  expect_true(result$rm == 2)
  result <- get_centers(m1, ploidy = 2, n.clusters.thr = 2)


  expect_type(result, "list")
  expect_named(result, c("rm", "centers_theta", "MarkerName", "n.clusters"))
  expect_true(length(result$centers_theta) == 3)
  expect_equal(result$centers_theta, c(0.0831, 0.5327, 0.9823), tolerance = 1e-4)
})

test_that("get_baf works as expected", {
  centers_theta <- c(0.1, 0.5, 0.9)
  theta <- c(0.05, 0.3, 0.6, 0.95)

  bafs <- get_baf(theta, centers_theta, ploidy = 2)

  expect_equal(length(bafs), length(theta))
  expect_true(all(bafs >= 0 & bafs <= 1))
  expect_equal(sum(bafs), 1.875, tolerance = 1e-4)
})

test_that("get_baf_par returns a list of correct length", {
  # test-get_baf_par.R

  library(testthat)
  library(dplyr)

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
})

test_that("get_zscore returns z-score table", {
  result <- get_zscore(data = sample_data, geno.pos = geno_pos)

  expect_s3_class(result, "data.frame")
  expect_true("z" %in% colnames(result))
  expect_equal(nrow(result), nrow(sample_data))
  expect_equal(var(result$z), 0.8163, tolerance = 0.001)
})


test_that("standardize runs and returns expected object", {
  result <- standardize(
    data = sample_data,
    genos = geno_data,
    geno.pos = geno_pos,
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
  expect_equal(sum(result$data$baf, na.rm = TRUE), 23.40818 ,tolerance = 0.001)
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
