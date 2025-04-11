#' Pascal triangle to define the expected peaks for each ploidy
#'
#' @param h ploidy value
#'
#' @export
pascalTriangle <- function(h) {
  lapply(0:h, function(i) choose(i, 0:i))
}

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


#' Simulate Genotyping Data with Flexible Ploidy
#'
#' Generates synthetic genotyping and signal intensity data for a given ploidy level.
#' Returns a structured list containing input data suitable for standardization analysis.
#'
#' @param n_markers Integer. Number of markers to simulate (default: 10).
#' @param n_samples Integer. Number of individuals/samples to simulate (default: 5).
#' @param ploidy Integer. Ploidy level of the organism (e.g., 2 for diploid, 4 for tetraploid).
#' @param seed Integer. Random seed for reproducibility (default: 2025).
#'
#' @return A named list with:
#' \describe{
#'   \item{sample_data}{Allelic signal intensities (X, Y, R, ratio).}
#'   \item{geno_data}{Genotype dosage and probability data.}
#'   \item{geno_pos}{Genomic coordinates for each marker.}
#'   \item{standardization_input}{Merged input data with theta and genotype.}
#' }
#'
#' @export
simulate_standardization_input <- function(n_markers = 10,
                                           n_samples = 5,
                                           ploidy = 2,
                                           seed = 2025) {
  set.seed(seed)

  marker_ids <- paste0("m", 1:n_markers)
  sample_ids <- paste0("S", 1:n_samples)

  data_grid <- expand.grid(MarkerName = marker_ids, SampleName = sample_ids)

  final_dataset <- data_grid %>%
    rowwise() %>%
    mutate(
      geno = sample(0:ploidy, 1),
      # Simulate allele signal as interpolation between "pure" profiles
      X = round(rnorm(1, mean = 280 - geno * ((280 - 20) / ploidy), sd = 10)),
      Y = round(rnorm(1, mean = 20 + geno * ((280 - 20) / ploidy),  sd = 10)),
      R = X + Y,
      ratio = round(Y / R, 4),
      prob = round(runif(1, 0.85, 1.00), 2)
    ) %>%
    ungroup()

  sample_data <- final_dataset %>%
    select(MarkerName, SampleName, X, Y, R, ratio)

  geno_data <- final_dataset %>%
    select(MarkerName, SampleName, geno, prob)

  geno_pos <- tibble::tibble(
    MarkerName = marker_ids,
    Chromosome = rep(c("1", "2", "3"), length.out = n_markers),
    Position = sort(sample(1e5:1e6, n_markers))
  )

  standardization_input <- inner_join(
    sample_data %>% select(MarkerName, SampleName, ratio) %>% rename(theta = ratio),
    geno_data %>% select(MarkerName, SampleName, geno),
    by = c("MarkerName", "SampleName")
  )

  list(
    sample_data = sample_data,
    geno_data = geno_data,
    geno_pos = geno_pos,
    standardization_input = standardization_input
  )
}
