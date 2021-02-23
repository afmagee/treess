#' Calculates ESS for a single summary variable obtained from classical multidimensional scaling.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples For compatibility with eval and call, not used here.
#' @keywords internal
CMDSESS <- function(dmat,trees=NA,nsim=NA,alpha=NA,min.nsamples=NA) {
  dmat <- dmat^2
  x <- stats::cmdscale(dmat,k=1)[,1]
  return(coda::effectiveSize(x))
}