% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardization.R
\name{read_qploidy_standardization}
\alias{read_qploidy_standardization}
\title{Read Qploidy Standardization File}
\usage{
read_qploidy_standardization(qploidy_standardization_file)
}
\arguments{
\item{qploidy_standardization_file}{A string specifying the path to the file generated
by the `standardize` function. The file should be in CSV format and include metadata,
filters, and data sections.}
}
\value{
An object of class `qploidy_standardization`, which is a list containing:
  - `info`: A named vector of standardization parameters.
  - `filters`: A named vector summarizing the number of markers removed at each filtering step.
  - `data`: A data frame containing the standardized dataset with BAF, Z-scores, and genotype information.
}
\description{
This function reads a file generated by the `standardize` function and reconstructs
a `qploidy_standardization` object. The file contains metadata, filtering information,
and the standardized dataset.
}
\details{
The function uses the `vroom` package to efficiently read the file in chunks.
The first row contains metadata (`info`), the second row contains filtering information (`filters`),
and the remaining rows contain the standardized dataset (`data`).
}
