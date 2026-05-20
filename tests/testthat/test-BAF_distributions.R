# Test suite for BAF_distributions.R functions
library(testthat)
library(Qploidy2)

test_that("generate_baf_template returns valid template for all distributions", {
  cn <- 4
  M <- 21
  set.seed(123)
  dists <- c("gaussian", "beta", "beta_binomial", "negative_binomial")
  for (dist in dists) {
    templ <- generate_baf_template(cn, M = M, bw = 0.03, dist = dist)
    expect_type(templ, "double")
    expect_length(templ, M)
    expect_true(all(templ > 0))
    expect_equal(sum(templ), 1, tolerance = 1e-6)
  }
})

test_that("baf_log_likelihood returns correct value for known input", {
  templ <- rep(1/10, 10)
  counts <- rep(10, 10)
  ll <- baf_log_likelihood(counts, templ)
  expect_equal(ll, -230.2585, tolerance = 1e-4)
  expect_equal(ll, sum(counts * log(templ)), tolerance = 1e-6)
  expect_equal(baf_log_likelihood(rep(0, 10), templ), 0)
})

test_that("compute_baf_likelihoods returns expected structure and probabilities sum to 1", {
  set.seed(1)
  baf_vec <- rbeta(200, 2, 2)
  cn_grid <- 2:4
  res <- compute_baf_likelihoods(baf_vec, cn_grid, M = 21, bw = 0.03, dist = "gaussian", plot = FALSE)
  expect_type(res, "list")
  expect_true(all(c("ll_vec", "prob_vec") %in% names(res)))
  expect_length(res$ll_vec, length(cn_grid))
  expect_length(res$prob_vec, length(cn_grid))
  expect_equal(sum(res$prob_vec), 1, tolerance = 1e-6)
  expect_equal(res$ll_vec[3], -947.7981, tolerance = 1e-4)
})

test_that("select_best_baf_model returns best model and grid results", {
  set.seed(2)
  baf_vec <- rbeta(200, 2, 2)
  cn_grid <- 2:4
  res <- select_best_baf_model(baf_vec = baf_vec,cn_grid =  cn_grid,
                               dists = c("gaussian", "beta"), bw_grid = c(0.02, 0.03),
                               add_uniform_grid = c(FALSE, TRUE), uniform_weight_grid = c(0.01, 0.05),
                               M = 21, plot = FALSE)
  expect_type(res, "list")
  expect_true(all(c("n_obs", "best", "grid_results") %in% names(res)))
  expect_true(is.data.frame(res$grid_results))
  expect_true(!is.null(res$best))
  expect_true(res$n_obs > 0)
  expect_true(res$best$best_cn %in% cn_grid)
  expect_true(res$best$dist %in% c("beta"))
  expect_true(res$best$bw %in% c(0.03))
  expect_true(res$best$add_uniform %in% c(FALSE, TRUE))
  expect_true(res$best$uniform_weight %in% c(0))
  # If plot=TRUE, returns a ggplot object
  res_plot <- select_best_baf_model(baf_vec = baf_vec, cn_grid = cn_grid,
                                    dists = c("gaussian"), bw_grid = 0.03,
                                    add_uniform_grid = FALSE, uniform_weight_grid = 0,
                                    M = 21, plot = TRUE)
  expect_true("plot" %in% names(res_plot))
  expect_true(any(class(res_plot$plot) %in% c("gg", "ggplot")))
})
