#' Computes ESS by finding the time lag at which the "jump distance" curve becomes indistinguishable from the jump distances of an uncorrelated set of samples from the same target distribution.
#' Uses permutation resampling to compute the null distribution.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim Number of simulations (bootstraps/permutations) to compute the null distribution.
#' @param alpha The cutoff quantile for the null distribution, above which samples are considered independent.
#' @param min.nsamples We only consider distances between x_t and x_{t+T} if there are at least this many samples.
#' @param interpolate: Whether to use linear interpolation (default) on the (smoothed) jump distance curve (otherwise it is piecewise constant).
#' @param use.median: Whether to summarize the jump distance distribution at each time lag by the median (if TRUE), otherwise the mean is used.
#' @keywords internal
jumpDistancePermutationESS <- function(dmat,trees=NA,min.nsamples=5,nsim=1000,alpha=0.05) {
  # A shortcut for jumpDistanceEquivalenceTestESS using permutations
  jumpDistanceEquivalenceTestESS(dmat=dmat,min.nsamples=min.nsamples,nsim=nsim,alpha=alpha,bootstrap=FALSE,interpolate=TRUE)
}

#' Computes ESS by finding the time lag at which the "jump distance" curve becomes indistinguishable from the jump distances of an uncorrelated set of samples from the same target distribution.
#' Uses (nonparametric) bootstrapping to compute the null distribution.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim Number of simulations (bootstraps/permutations) to compute the null distribution.
#' @param alpha The cutoff quantile for the null distribution, above which samples are considered independent.
#' @param min.nsamples We only consider distances between x_t and x_{t+T} if there are at least this many samples.
#' @param interpolate: Whether to use linear interpolation (default) on the (smoothed) jump distance curve (otherwise it is piecewise constant).
#' @param use.median: Whether to summarize the jump distance distribution at each time lag by the median (if TRUE), otherwise the mean is used.
#' @keywords internal
jumpDistanceBootstrapESS <- function(dmat,trees=NA,min.nsamples=5,nsim=1000,alpha=0.05) {
  # A shortcut for jumpDistanceEquivalenceTestESS using bootstrapping
  jumpDistanceEquivalenceTestESS(dmat=dmat,min.nsamples=min.nsamples,nsim=nsim,alpha=alpha,bootstrap=TRUE,interpolate=TRUE)
}

#' Computes ESS by finding the time lag at which the "jump distance" curve becomes indistinguishable from the jump distances of an uncorrelated set of samples from the same target distribution.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim Number of simulations (nonparametric bootstraps/permutations) to compute the null distribution.
#' @param alpha The cutoff quantile for the null distribution, above which samples are considered independent.
#' @param min.nsamples We only consider distances between x_t and x_{t+T} if there are at least this many samples.
#' @param bootstrap: Whether to use nonparametric bootstrapping (if TRUE) for resampling or permutation resampling (if FALSE).
#' @param interpolate: Whether to use linear interpolation (default) on the (smoothed) jump distance curve (otherwise it is piecewise constant).
#' @param use.median: Whether to summarize the jump distance distribution at each time lag by the median (if TRUE), otherwise the mean is used.
#' @keywords internal
jumpDistanceEquivalenceTestESS <- function(dmat,min.nsamples,nsim,alpha,bootstrap,interpolate,use.median=TRUE) {
  # recover()
  
  central_tendency <- mean
  if ( use.median ) {
    central_tendency <- median
  }
  
  n <- dim(dmat)[1]
  
  # permute matrices to get appropriate autocorrelation
  off_diag <- row(dmat) == col(dmat)+1
  null_dist <- sapply(1:nsim,function(i){
    permute <- sample.int(n,replace=bootstrap)
    central_tendency(sapply(1:(n-1),function(i){dmat[permute[i],permute[i+1]]}))
  })
  threshold <- quantile(null_dist,probs=alpha)
  
  # get distances at time lags until we hit the threshold
  lagged <- numeric(n-min.nsamples+1)
  for (i in 1:(n-min.nsamples+1)) {
    # distance at this time lag
    this_lag <- central_tendency(dmat[row(dmat) == col(dmat)+i])
    # ensure distance by lag curve is monotonic
    lagged[i] <- ifelse(i == 1,this_lag,max(this_lag,lagged[i-1]))
    # lagged[i] <- ifelse(i == 1,this_lag,max(this_lag,lagged[i-1]))
    # early terimination to avoid unneeded computation
    if (lagged[i] > threshold) {
      break
    }
  }
  
  first_larger <- sapply(threshold,function(thresh){
    min(which(lagged > thresh))
  })
  
  thin <- numeric(length(alpha))
  names(thin) <- paste0(alpha*100,"%")
  for (i in 1:length(alpha)) {
    if (interpolate && first_larger[i] > 1) {
      x2 <- first_larger[i]
      x1 <- x2 - 1
      y2 <- lagged[x2]
      y1 <- lagged[x1]
      slope <- (y2 - y1)/(x2 - x1)
      thin[i] <- x1 + (threshold[i] - y1)/slope
    } else {
      thin[i] <- first_larger[i]
    }
  }
  
  return(n/thin)
  
}
