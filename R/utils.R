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


#' print qploidy_area_ploidy_estimation object
#'
#' @param x qploidy_area_ploidy_estimation object
#' @param ... print parameters
#'
#' @method print qploidy_area_ploidy_estimation
#'
#' @return No return value, called for side effects.
#'
#' @export
print.qploidy_area_ploidy_estimation <- function(x, ...){

  count_aneu <- get_aneuploids(x$ploidy)

  df <- data.frame(c1 = c("Number of samples:",
                          "Chromosomes:",
                          "Tested ploidies:",
                          "Number of euploid samples:",
                          "Number of potential aneuploid samples:",
                          "Number of highly inbred samples:"),
                   c2 = c(dim(x$ploidy)[1],
                          {if(!is.null(colnames(x$ploidy)))
                            paste0(colnames(x$ploidy), collapse = ",") else
                              paste0(x$chr, collapse = ",")},
                          paste0(x$tested, collapse = ","),
                          sum(!count_aneu, na.rm = TRUE),
                          sum(count_aneu, na.rm = TRUE),
                          x$n.inbred
                   ))

  colnames(df) <- rownames(df) <- NULL

  cat("Object of class qploidy_area_ploidy_estimation")
  print(format(df, justify = "left", digits = 2))
}


#' indexes for aneuploids
#'
#' @param ploidy_df ploidy table (chromosome in columns and individuals in rows)
#'
#' @return A logical vector where each element corresponds to an individual in the
#'         input ploidy table. The value is `TRUE` if the individual is identified
#'         as potentially aneuploid, and `FALSE` otherwise.
#'
#' @export
get_aneuploids  <- function(ploidy_df){

  if(any(grepl("/NA", ploidy_df) | grepl("NA/",ploidy_df))){
    ploidy_df <- gsub("/NA", "", ploidy_df)
    ploidy_df <- gsub("NA/", "", ploidy_df)
  }

  count_aneu <- !apply(ploidy_df, 1, function(y) {
    if(any(is.na(y))) {
      temp <- length(unique(y[-which(is.na(y))]))
      if(temp == 1) TRUE else if(temp == 0) NA else FALSE
    } else length(unique(y)) == 1
  })
  return(count_aneu)
}

##' Verbose Message Utility
##'
##' Prints a formatted verbose message with timestamp, indentation, and type label, if verbose is TRUE.
##'
##' @param text Character string, the message to print (supports sprintf formatting).
##' @param verbose Logical. If TRUE, prints the message; if FALSE, suppresses output.
##' @param level Integer, indentation level (0=header, 1=main step, 2=detail, 3=sub-detail).
##' @param type Character string, message type (e.g., "INFO", "WARN", "ERROR"). Only shown for level 0.
##' @param ... Additional arguments passed to sprintf for formatting.
##'
##' @details Use the verbose argument to control message output. Typically, pass the function's verbose parameter to vmsg.
##'
##' @return No return value, called for side effects.
##' @export
vmsg <- function(text, verbose = FALSE, level = 1, type = ">>", ...) {
  if (!verbose) return(invisible())
  # Format timestamp
  timestamp <- format(Sys.time(), "[%H:%M:%S]")

  # Create indentation based on level
  indent <- switch(as.character(level),
    "0" = "",           # Section headers
    "1" = "  * ",       # Main steps (medium bullet)
    "2" = "    - ",     # Details
    "3" = "      > ",   # Sub-details
    paste0(paste(rep("  ", level), collapse = ""), "* ")  # Fallback for level > 3
  )

  # Format type label (only show for level 0)
  type_label <- if (level == 0) sprintf("%-1s ", type) else ""

  # Format message text
  dots <- list(...)
  if(length(dots) == 0) {
    msg_text <- text
  } else {
    msg_text <- sprintf(text, ...)
  }
  # Combine everything
  formatted_msg <- sprintf("%s %s%s%s", timestamp, type_label, indent, msg_text)
  message(formatted_msg)
}