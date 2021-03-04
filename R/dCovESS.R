#' # Computes ESS by finding the time lag at which the distance covariance is not distinguishable from 0
#' # Uses the energy package implementation which is written in C for speed
#' #' Computes ESS by finding the time lag at which the "jump distance" curve becomes indistinguishable from the jump distances of an uncorrelated set of samples from the same target distribution.
#' #' Uses (nonparametric) bootstrapping to compute the null distribution.
#' #'
#' #' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' #' @param trees For compatibility with eval and call, not used here.
#' #' @param nsim Number of simulations (bootstraps/permutations) to compute the null distribution.
#' #' @param alpha The cutoff quantile for the null distribution, above which samples are considered independent.
#' #' @param min.nsamples We only consider distances between x_t and x_{t+T} if there are at least this many samples.
#' #' @param interpolate: Whether to use linear interpolation (default) on the (smoothed) jump distance curve (otherwise it is piecewise constant).
#' #' @param use.median: Whether to summarize the jump distance distribution at each time lag by the median (if TRUE), otherwise the mean is used.
#' #' @keywords internal
#' dCovESS <- function(dmat,trees=NA,min.nsamples=5,nsim=1000,alpha=0.05,bootstrap=TRUE,interpolate=TRUE) {
#'   # recover()
#'   
#'   n <- dim(dmat)[1]
#'   
#'   # permute matrices to get appropriate autocorrelation
#'   null_dist <- sapply(1:nsim,function(i){
#'     permute <- permuteDistanceMatrix(dmat)
#'     energy::dcov(permute[1:(n-1),1:(n-1)],permute[2:n,2:n],index=1.0)
#'   })
#'   threshold <- quantile(null_dist,probs=alpha)
#'   
#'   # guarantees if we don't find a better thinning value we set ESS = 1
#'   thin <- n 
#'   dc <- rep(NA,n-min.nsamples+1)
#'   # Use -dCov() such that curve is increasing (nondecreasing) instead of decreasing (nonincreasing)
#'   dc[1] <- -energy::dcov(dmat[1:(n-1),1:(n-1)],dmat[2:n,2:n],index=1.0)
#'   for (i in 2:(n-min.nsamples+1)) {
#'     dc[i] <- max(dc[i-1],-energy::dcov(dmat[-c(1:i),-c(1:i)],dmat[-c((n-i+1):n),-c((n-i+1):n)],index=1.0))
#'     if ( dc[i] > threshold ) {
#'       thin <- i
#'       break
#'     }
#'   }
#'   
#'   if ( interpolate && thin != n && thin != 1 ) {
#'     # interpolation assumes we have first entry at lag 0
#'     dc <- c(energy::dcov(dmat,dmat),dc[!is.na(dc)])
#'     thin <- finds0Smoothed(x=0:(length(dc)-1),y=dc,threshold)
#'   }
#'   
#'   return(n/thin)
#'   
#' }
#' 
#' # sumOfdCorESS <- function(dmat) {
#' #   n <- dim(dmat)[1]
#' #   cors <- sapply(1:(n-5+1),function(i){
#' #     energy::bcdcor(dmat[-c(1:i),-c(1:i)],dmat[-c((n-i+1):n),-c((n-i+1):n)])
#' #   })
#' #   cors <- c(1,cors)
#' #   for (i in 2:length(cors)) {
#' #     cors[i] <- min(cors[i-1],cors[i])
#' #   }
#' #   cmax <- length(cors)
#' #   if (any(cors < 0)) {
#' #     cmax <- min(which(cors < 0)) - 1
#' #   }
#' #   return(n/sum(cors[1:cmax]))
#' # }