#' Calculates the ESS of a transformed variable y where y[i] = sum_j(d(x[i],x[j])).
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples For compatibility with eval and call, not used here.
#' @keywords internal
totalDistanceESS <- function(dmat,...) {
  # recover()
  
  x <- rowSums(dmat)
  
  return(coda::effectiveSize(x))
}