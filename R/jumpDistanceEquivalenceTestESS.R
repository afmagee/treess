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
  
  # guarantees if we don't find a better thinning value we set ESS = 1
  thin <- n 
  # get distances at time lags until we hit the threshold
  G_s <- rep(NA,n-min.nsamples+1)
  G_s[1] <- central_tendency(dmat[row(dmat) == col(dmat)+1])
  for (i in 2:(n-min.nsamples+1)) {
    # distance at this time lag
    g_s <- central_tendency(dmat[row(dmat) == col(dmat)+i])
    G_s[i] <- max(g_s,G_s[i-1])
    # early terimination to avoid unneeded computation
    if (G_s[i] > threshold) {
      thin <- i
      break
    }
  }
  
  if ( interpolate && thin != n && thin != 1 ) {
    # interpolation assumes we have first entry at lag 0
    G_s <- c(0,G_s[!is.na(G_s)])
    thin <- finds0Smoothed(x=0:(length(G_s)-1),y=G_s,threshold)
  }

  return(n/thin)
  
}
