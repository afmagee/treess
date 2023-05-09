#' Calculates the approximate ESS of Lanfear et al. (2015).
#'
#' @param dmat Sample-to-sample distance matrix for chains, NOT SQUARED.
#' @param max.approximateESS.timelag The maximum time lag for computing jump distances
#' @param alpha The cutoff proportion of the asymptote used to estimate the average pairwise squared distance.
#' @keywords internal
approximateESS <- function(dmat,max.approximateESS.timelag=100,alpha=0.05,...) {
  # This function essentially strings together the following RWTY functions
  # topological.autocorr, topological.approx.ess, approx.ess.multi, approx.ess.single, tree.autocorr

  dmat <- dmat^2
  
  N <- n <- dim(dmat)[1]
  
  # make sure we don't over shoot
  max_off_diag <- max.approximateESS.timelag
  if ( max_off_diag > dim(dmat)[1]-1 ) {
    max_off_diag <- dim(dmat)[1]-1
  }
  
  # get empirical curve of lag time t vs distance
  # see also topological.autocorr
  t <- 1:max_off_diag
  d_t <- sapply(t,function(t_){
    mean(dmat[row(dmat) == col(dmat) + t_])
  })
  
  m <- .computeApproximateESSThreshold(d_t, t, alpha)
  
  return(.computeApproxESSFromComponents(N, m, d_t))
}

#' Helper function for approximate ESS
#'
#' @param d_t Average (squared) distances at time lags
#' @param t The time lags
#' @param alpha 1 - threshold
#' @keywords internal
.computeApproximateESSThreshold <- function(d_t, t, alpha) {
  # alpha to threshold
  thresh = 1 - alpha
  
  # Optimize the exponential curve fit
  fn <- function(par,t,d_t) {
    E_d <- par[1] * (1 - exp(-t/par[2]))
    return( sum((E_d - d_t)^2) )
  }
  par <- optim(c(1,1),fn,t=t,d_t=d_t)
  
  # For finding the cutoff lag
  thresh <- thresh * par$par[1]
  
  m <- ifelse( max(d_t) < thresh, length(d_t)+1, min(which(d_t >= thresh)))
  
  return(m)
}

#' Helper function for approximate ESS
#'
#' @param N Number of samples in the chain
#' @param m The time lag at which samples become independent
#' @param d_t Average (squared) distances at time lags
#' @keywords internal
.computeApproxESSFromComponents <- function(N, m, d_t) {
  # The below code is from approx.ess.single, with d_t in place of df$topo.distance 
  D <- max(d_t)
  
  S = 0
  
  if(m>1){
    for(k in 1:(m - 1)){
      f = d_t[k]
      S = S + ((N - k) * f)
    }
  }
  
  S = S + (N - m + 1) * (N - m) * D / 2
  S = S / 2 / N^2
  ESS = 1 / (1 - 4 * S / D)
  
  return(ESS)
}

#' Calculates the approximate ESS of Lanfear et al. (2015).
#'
#' Uses subsampling to avoid computing the full distance matrix, there are several options.
#' See details
#'
#' @param processed.trees The phylogenies, pre-processed for fast distance computation.
#' @param dist.fn A function for computing distances on the already pre-processed trees
#' @param n.per.diag The number of samples taken per diagonal band.
#' @param max.approximateESS.timelag The maximum time lag for computing jump distances.
#' @param alpha The cutoff proportion of the asymptote used to estimate the average pairwise squared distance.
#' @details Depending on how n.per.diag and max.approximateESS.timelag are set, there are four possible behaviors.
#' Note that n.per.diag >= length(processed.trees) is a shortcut for "use every sample on the diagonal", and max.approximateESS.timelag >= length(processed.trees) - 1 means there is no truncation of the computation.
#' 0. n.per.diag = Inf, max.approximateESS.timelag = Inf: No subsampling, whole distance matrix is used.
#' 1. n.per.diag = Inf, max.approximateESS.timelag = 100: No subsampling, whole distance matrix is used.
#' 2. n.per.diag = 100 , max.approximateESS.timelag = Inf: No subsampling, whole distance matrix is used.
#' 3. n.per.diag = 100, max.approximateESS.timelag = 100: No subsampling, whole distance matrix is used.
#' Option 1 recapitulates the behavior of treess::approximateESS, while option 3 is more akin to the defaults in RWTY.
#' Experiments suggest that option 2 produces notably more pessimistic estimates than 0, 1, or 3.
#' @keywords internal
subsampledApproximateESS <- function(processed.trees,dist.fn,n.per.diag,max.approximateESS.timelag=100,alpha=0.05,...) {
  i_j_dist <- computeDiagonallySubsampledDistances(processed.trees,dist.fn,n.per.diag,same.n.per.diag=FALSE,max.diag=max.approximateESS.timelag)
  d_t <- unlist(lapply(i_j_dist,function(xx){mean(xx[,3]^2)}))
  
  N <- n <- length(processed.trees)
  
  # make sure we don't over shoot
  max_off_diag <- max.approximateESS.timelag
  if ( max_off_diag > n-1 ) {
    max_off_diag <- n-1
  }
  
  t <- 1:max_off_diag
  
  m <- .computeApproximateESSThreshold(d_t, t, alpha)
  
  return(.computeApproxESSFromComponents(N, m, d_t))
  
}
