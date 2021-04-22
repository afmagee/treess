#' Calculates ESS for a single summary variable obtained from classical multidimensional scaling.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @keywords internal
CMDSESS <- function(dmat,...) {
  dmat <- dmat^2
  x <- stats::cmdscale(dmat,k=1)[,1]
  return(coda::effectiveSize(x))
}