#' EM loop for HMM ploidy estimation
#'
#' This function runs the EM algorithm for the HMM ploidy estimation, given all required parameters.
#' Returns a list with updated parameters and results.
#'
#' @param cn_grid Integer vector of copy-number states to consider
#' @param mu Initial means for each state
#' @param K Number of states
#' @param state_ids Character vector of state IDs
#' @param sig Initial shared standard deviation
#' @param z Numeric vector of z-scores
#' @param z_only Logical, use z emission only
#' @param ll_baf_matrix BAF log-likelihood matrix
#' @param n_baf Numeric vector, BAF normalization
#' @param w_baf Numeric, BAF weight
#' @param correct_scale Logical, correct BAF scale
#' @param A Initial transition matrix
#' @param pi0 Initial state probabilities
#' @param W Number of windows
#' @param max_iter Maximum EM iterations
#' @param verbose Logical, print progress
#'
#' @return List with updated parameters: mu, cn_grid, K, state_ids, sig, gamma, ll_em, pi0, A, ll_hist
em_hmm_cn <- function(cn_grid, mu, K, state_ids, sig, z, z_only, ll_baf_matrix, n_baf, w_baf, correct_scale, A, pi0, W, max_iter, verbose) {
  ll_hist <- numeric(max_iter)
  for (iter in 1:max_iter) {
    # Emissions
    ll_em <- matrix(NA_real_, nrow=W, ncol=K, dimnames=list(NULL, state_ids))
    for (k in seq_len(K)) {
      c <- cn_grid[k]
      llz <- dnorm(z, mean=mu[as.character(c)], sd=sig, log=TRUE)
      if(any(is.nan(llz))) llz[which(is.nan(llz))] <- 0
      if (z_only) {
        ll_em[,k] <- llz
      } else {
        llb <- ll_baf_matrix[,k]
        if(correct_scale) {
          llb <- llb/ n_baf
        }
        ll_em[, k] <- (1-w_baf) * llz + w_baf * llb
      }
    }
    if (!all(is.finite(ll_em))) {
      bad_w <- which(!is.finite(rowSums(ll_em)))[1]
      bad_k <- which(!is.finite(ll_em[bad_w, ]))
      stop(sprintf("Non-finite emission at window %d, states: %s.", bad_w, paste(colnames(ll_em)[bad_k], collapse=", ")))
    }
    logA <- log(A); logpi0 <- log(pi0)
    log_alpha <- matrix(-Inf, W, K); log_beta <- matrix(0, W, K)
    # forward
    log_alpha[1, ] <- logpi0 + ll_em[1, ]
    for (i in 2:W) {
      for (k in 1:K) {
        log_alpha[i,k] <- ll_em[i,k] + logsumexp(log_alpha[i-1, ] + logA[,k])
      }
    }
    # backward
    for (i in (W-1):1) {
      for (k in 1:K) {
        log_beta[i,k] <- logsumexp(logA[k, ] + ll_em[i+1, ] + log_beta[i+1, ])
      }
    }
    loglik <- logsumexp(log_alpha[W, ])
    ll_hist[iter] <- loglik
    # E-step: posteriors
    log_gamma <- log_alpha + log_beta
    log_gamma <- sweep(log_gamma, 1, apply(log_gamma, 1, logsumexp), "-")
    gamma <- exp(log_gamma)
    # pairwise
    xi_sum <- matrix(0, K, K)
    for (i in 1:(W-1)) {
      M_ij <- outer(log_alpha[i, ], log_beta[i+1, ], "+") +
        logA + matrix(ll_em[i+1, ], K, K, byrow=TRUE)
      M_ij <- M_ij - logsumexp(as.vector(M_ij))
      xi_sum <- xi_sum + exp(M_ij)
    }
    # M-step: update parameters
    pi0 <- gamma[1, ] / sum(gamma[1, ])
    A <- xi_sum / pmax(rowSums(xi_sum), 1e-12)
    A[!is.finite(A)] <- 0
    A <- sweep(A, 1, pmax(rowSums(A), 1e-12), "/")
    A <- pmax(A, 1e-12); A <- sweep(A, 1, rowSums(A), "/")
    # update mu and sigma
    mu <- numeric(K)
    for (k in 1:K) {
      w <- gamma[,k]
      mu[k] <- sum(w * z) / pmax(sum(w), 1e-12)
    }
    mu <- setNames(mu, as.character(cn_grid))
    mu <- mu[state_ids]
    # update shared sigma
    sig <- sqrt(sum(gamma * (matrix(z, W, K) - rep(mu, each=W))^2) /
                  pmax(sum(gamma), 1e-12))
    sig <- max(sig, 1e-3)
    # Convergence check
    if (iter > 4 && is.finite(ll_hist[iter]) && is.finite(ll_hist[iter-1]) &&
        abs(ll_hist[iter] - ll_hist[iter-1]) < 1e-4) break
  }
  vmsg("Testing CN grid %s:", verbose = verbose, level = 2, type = ">>", paste(cn_grid, collapse=", "))
  vmsg("EM converged in %d iterations", verbose = verbose, level = 2, type = ">>", iter)
  vmsg("Final log-likelihood: %.2f", verbose = verbose, level = 2, type = ">>", ll_hist[iter])

  return(list(mu=mu, cn_grid=cn_grid, K=K, state_ids=state_ids, sig=sig, gamma=gamma, ll_em=ll_em, pi0=pi0, A=A, ll_hist=ll_hist))
}
