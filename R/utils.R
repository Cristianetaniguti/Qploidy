globalVariables(c("rnom", "X", "Y", "runif", "prob"))

#' Pascal Triangle for Expected Peaks Calculation
#'
#' This function generates the Pascal triangle for a given ploidy value. The Pascal triangle is used to define the expected peaks for each ploidy level, which can be useful in various genetic analyses.
#'
#' @param h An integer representing the ploidy value.
#'
#' @return A list where each element corresponds to a row of the Pascal triangle, up to the specified ploidy value.
#'
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


#' Check if a File is Compressed
#'
#' This function checks whether a given file is compressed by inspecting its
#' magic number (first few bytes). It detects common compression formats such as
#' gzip (`.gz`), bzip2 (`.bz2`), and xz (`.xz`).
#'
#' @param file_path A character string giving the path to the file to be checked.
#'
#' @return A character string indicating the compression type (`"gzip (.gz)"`,
#' `"bzip2 (.bz2)"`, or `"xz (.xz)"`) if the file is compressed, or `FALSE` if
#' the file is not recognized as compressed.
#'
is_compressed_file <- function(file_path) {
  con <- file(file_path, "rb")
  magic <- readBin(con, "raw", n = 4)
  close(con)

  # Known magic numbers for compressed formats
  if (all(magic[1:2] == as.raw(c(0x1f, 0x8b)))) {
    return("gzip (.gz)")
  } else if (all(magic[1:3] == as.raw(c(0x42, 0x5a, 0x68)))) {
    return("bzip2 (.bz2)")
  } else if (all(magic[1:6] == as.raw(c(0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00)))) {
    return("xz (.xz)")
  } else {
    return(FALSE)  # Not a known compressed format
  }
}
