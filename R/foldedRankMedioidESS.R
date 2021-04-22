#' Calculates an ESS in the spirit of the folded rank-transformed ESS of Vehtari et al. (2019), replacing the median with the medioid.
#' Where there is no unique medioid, calculates the ESS for all and returns the minimum.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @keywords internal
foldedRankMedioidESS <- function(dmat,...) {
  
  # there may not be one unique medioid
  the_medioids <- which(rowSums(dmat) == min(rowSums(dmat)))
  
  all_ess <- sapply(the_medioids,function(the_medioid) {
    xr <- rank(dmat[,the_medioid])
    X <- qnorm((xr-0.375)/(max(xr)-0.25))
    coda::effectiveSize(X)
  })
  
  return(min(all_ess))
}
