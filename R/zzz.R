.onAttach <- function(libname, pkgname){
  msg <- paste(
    "Migration Notice: Support for this package on CRAN ends 12/30/2026.",
    "Please migrate to **Qploidy2**: https://github.com/Breeding-Insight/Qploidy2"
  )

  packageStartupMessage(msg)
}
