#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' Instead of calculating the pseudo-ESS 100 times and reporting summaries, returns the minimum over all possible reference trees.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples For compatibility with eval and call, not used here.
#' @keywords internal
minPseudoESS <- function(dmat,trees=NA,nsim=NA,alpha=NA,min.nsamples=NA) {
  # Shortcut for pseudoESS(summary="min")
  pseudoESS(dmat,summary="min")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' Instead of calculating the pseudo-ESS 100 times and reporting summaries, returns the median over all possible reference trees.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples For compatibility with eval and call, not used here.
#' @keywords internal
medianPseudoESS <- function(dmat,trees=NA,nsim=NA,alpha=NA,min.nsamples=NA) {
  # Shortcut for pseudoESS(summary="median")
  pseudoESS(dmat,summary="median")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' Instead of calculating the pseudo-ESS 100 times and reporting summaries, returns the minimum or median over all possible reference trees at user's discretion.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param summary Should we return the minimum or the median?
#' @keywords internal
pseudoESS <- function(dmat,summary=c("min","median")) {
  # recover()
  
  all_ess <- apply(dmat,1,coda::effectiveSize)
  
  ess <- 0
  if ( summary == "min" ) {
    ess <- min(all_ess)
  } else if ( summary == "median" ) {
    ess <- median(all_ess)
  } else {
    stop("Unrecognized summary option for pseudo-ESS")
  }
  
  return(ess)
}
