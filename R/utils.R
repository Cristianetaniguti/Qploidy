pascalTriangle <- function(h) {
  lapply(0:h, function(i) choose(i, 0:i))
}

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
