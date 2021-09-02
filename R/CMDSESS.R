#' Calculates ESS for a single summary variable obtained from classical multidimensional scaling.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @keywords internal
CMDSESS <- function(dmat,...) {
  dmat <- dmat^2
  # cmdscale does not work well if all samples are the same
  # So, we pre-empt this and declare the ESS to be 1 (there's 1 distinct value)
  if ( all(dmat == 0) ) {
    return(1.0)
  }
  x <- stats::cmdscale(dmat,k=1)[,1]
  return(coda::effectiveSize(x))
}