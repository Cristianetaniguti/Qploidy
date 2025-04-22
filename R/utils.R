globalVariables(c("rnom", "X", "Y", "runif", "prob"))

#' Pascal triangle to define the expected peaks for each ploidy
#'
#' @param h ploidy value
#'
#' @export
pascalTriangle <- function(h) {
  lapply(0:h, function(i) choose(i, 0:i))
}


#' Calculate the Statistical Mode
#'
#' This function returns the most frequent (modal) value in a vector.
#' If there are multiple values with the same highest frequency,
#' it returns the first one encountered.
#'
#' @param x A vector of numeric, character, or factor values.
#'
#' @return A single value representing the mode of the input vector.
#'
#' @examples
#' mode(c(1, 2, 2, 3, 4))   # Returns 2
#' mode(c("apple", "banana", "apple"))  # Returns "apple"
#'
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#' Adjust allele dosage using sequencing error and allele bias - function adapted from updog
#'
#' Adjusts allele dosage `p` by the sequencing error rate `eps`
#' and the allele bias `h`. This reproduces the behavior of the
#' `xi_fun()` from the updog package.
#'
#' @param p A numeric vector of allele dosages (e.g., 0 to 1).
#' @param eps A numeric vector (or scalar) of sequencing error rates.
#' @param h A numeric vector (or scalar) of allele bias values.
#'
#' @return A numeric vector of adjusted probabilities of the reference read.
#'
#' @examples
#' xi_fun(p = c(0, 0.5, 1), eps = 0.01, h = 0.8)
#'
#' @author David Gerard
#' @keywords internal
#' @noRd
xi_fun <- function(p, eps, h) {
  n <- length(p)

  if (!(length(eps) %in% c(1, n))) {
    stop("xi_fun: eps must either have length 1 or the same length as p.")
  }

  if (!(length(h) %in% c(1, n))) {
    stop("xi_fun: h must either have length 1 or the same length as p.")
  }

  # Vector recycling
  eps <- rep(eps, length.out = n)
  h <- rep(h, length.out = n)

  # Compute adjusted probabilities
  xi <- p * (1 - h) + (1 - p) * h * (1 - eps)

  return(xi)
}


