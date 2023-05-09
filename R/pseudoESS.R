#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
#' Instead of calculating the pseudo-ESS 100 times and reporting summaries, returns the minimum over all possible reference trees.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @keywords internal
minPseudoESS <- function(dmat,...) {
  # Shortcut for pseudoESS(summary="min")
  pseudoESS(dmat,summary="min")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
#' Instead of calculating the pseudo-ESS 100 times and reporting summaries, returns the median over all possible reference trees.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @keywords internal
medianPseudoESS <- function(dmat,...) {
  # Shortcut for pseudoESS(summary="median")
  pseudoESS(dmat,summary="median")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
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

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
#' Unlike calling the \link{treess} function with method="medianPseudoESS" or "minPseudoESS", this function only uses a subsample.
#' The default number of subsamples is 100, as in Lanfear et al. (2015).
#'
#' @param trees The phylogenies, pre-processed for fast distance computation.
#' @param dist.fn A function for computing distances on the already pre-processed trees
#' @param n.subsample.references The number of reference posterior samples to use. Default 100.
#' @keywords internal
subsampledMinPseudoESS <- function(processed.trees,dist.fn,n.subsample.references=100) {
  subsampledPseudoESS(processed.trees,dist.fn,n.subsample.references,"min")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
#' Unlike calling the \link{treess} function with method="medianPseudoESS" or "minPseudoESS", this function only uses a subsample.
#' The default number of subsamples is 100, as in Lanfear et al. (2015).
#'
#' @param trees The phylogenies, pre-processed for fast distance computation.
#' @param dist.fn A function for computing distances on the already pre-processed trees
#' @param n.subsample.references The number of reference posterior samples to use. Default 100.
#' @keywords internal
subsampledMedianPseudoESS <- function(processed.trees,dist.fn,n.subsample.references=100) {
  subsampledPseudoESS(processed.trees,dist.fn,n.subsample.references,"median")
}

#' Calculates the pseudo-ESS of Lanfear et al. (2015).
#' 
#' Unlike calling the \link{treess} function with method="medianPseudoESS" or "minPseudoESS", this function only uses a subsample.
#' The default number of subsamples is 100, as in Lanfear et al. (2015).
#'
#' @param processed.trees The phylogenies, pre-processed for fast distance computation.
#' @param dist.fn A function for computing distances on the already pre-processed trees
#' @param n.subsample.references The number of reference posterior samples to use. Default 100.
#' @param summary Should we return the minimum or the median?
#' @keywords internal
subsampledPseudoESS <- function(processed.trees,dist.fn,n.subsample.references=100,summary=c("min","median")) {
  # recover()
  n <- length(processed.trees)
  
  idx <- round(seq(1,length(processed.trees),length.out=n.subsample.references))
  idx <- do.call(rbind,lapply(idx,function(i){
    cbind(rep(i,n),1:n)
  }))
  all_dists <- getSubsampledDistanceMatrix(processed.trees,dist.fn,idx)
  
  all_ess <- sapply(1:n.subsample.references,function(i){
    coda::effectiveSize(all_dists[((i-1)*n+1):(i*n)])
  })
  
  ess <- 0
  if ( "min" %in% summary && "median" %in% summary ) {
    ess <- c(min(all_ess),median(all_ess))
    names(ess) <- c("min","median")
  } else if ( summary == "min" ) {
    ess <- min(all_ess)
    names(ess) <- "min"
  } else if ( summary == "median" ) {
    ess <- median(all_ess)
    names(ess) <- "median"
  } else {
    stop("Unrecognized summary option for pseudo-ESS")
  }
  
  return(ess)
}
