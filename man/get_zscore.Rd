% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardization.R
\name{get_zscore}
\alias{get_zscore}
\title{Calculate Z-Scores for Allele Intensities or Counts}
\usage{
get_zscore(data = NULL, geno.pos = NULL)
}
\arguments{
\item{data}{A data.frame containing signal intensity and ratio values with the following columns:
\describe{
  \item{MarkerName}{Marker identifiers.}
  \item{SampleName}{Sample identifiers.}
  \item{X}{Reference allele intensity or count.}
  \item{Y}{Alternative allele intensity or count.}
  \item{R}{Total signal or depth (i.e., \code{X + Y}).}
  \item{ratio}{Allelic ratio, typically \code{Y / (X + Y)}.}
}}

\item{geno.pos}{A data.frame with marker genomic positions, containing the following columns:
\describe{
  \item{MarkerName}{Marker identifiers.}
  \item{Chromosome}{Chromosome identifier where the marker is located.}
  \item{Position}{Genomic position (base-pair coordinate) of the marker.}
}}
}
\value{
A data.frame containing the following columns:
\describe{
  \item{MarkerName}{Marker ID.}
  \item{Chr}{Chromosome corresponding to the marker.}
  \item{Position}{Genomic position (bp).}
  \item{SampleName}{Sample ID.}
  \item{z}{Z-score computed per marker across all samples.}
}
Markers with missing chromosome or position information are excluded from the final output.
}
\description{
This function computes per-marker Z-scores based on the total signal intensity (R),
which typically represents the sum of reference (X) and alternative (Y) allele signals.
The Z-score measures how much each sample deviates from the mean intensity of that marker.
}
\details{
The function also merges positional metadata from the `geno.pos` input, adding chromosome
and physical position for each marker.
}
\examples{
data <- data.frame(
  MarkerName = rep("m1", 5),
  SampleName = paste0("S", 1:5),
  X = c(100, 110, 90, 95, 85),
  Y = c(200, 190, 210, 205, 215),
  R = c(300, 300, 300, 300, 300),
  ratio = c(0.67, 0.63, 0.70, 0.68, 0.72)
)
geno.pos <- data.frame(MarkerName = "m1", Chromosome = "1", Position = 123456)
get_zscore(data, geno.pos)

}
