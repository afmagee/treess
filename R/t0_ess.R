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
    if (lagged[i] >= threshold) {
      break
    }
  }
  
  first_larger <- sapply(threshold,function(thresh){
    min(which(lagged >= thresh))
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

# Computes ESS by finding the time lag at which the distance covariance is not distinguishable from 0
# Uses the energy package implementation which is written in C for speed
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
dCovESS <- function(dmat,trees=NA,min.nsamples=5,nsim=1000,alpha=0.05,bootstrap=TRUE,interpolate=TRUE) {
  # recover()
  
  n <- dim(dmat)[1]
  
  # permute matrices to get appropriate autocorrelation
  null_dist <- sapply(1:nsim,function(i){
    permute <- permuteDistanceMatrix(dmat)
    energy::dcov(permute[1:(n-1),1:(n-1)],permute[2:n,2:n],index=1.0)
  })
  threshold <- quantile(null_dist,probs=alpha)

  thin <- n #guarantees if we don't find a better thinning value we set ESS = 1
  dc <- rep(NA,n-min.nsamples)
  dc[1] <- energy::dcov(dmat[1:(n-1),1:(n-1)],dmat[2:n,2:n],index=1.0)
  for (i in 2:(n-min.nsamples)) {
    dc[i] <- min(dc[i-1],energy::dcov(dmat[-c(1:i),-c(1:i)],dmat[-c((n-i+1):n),-c((n-i+1):n)],index=1.0))
    if ( dc[i] <= threshold ) {
      thin <- i
      break
    }
  }
  
  if ( interpolate && thin != n ) {
    dc <- c(energy::dcov(dmat,dmat),dc[!is.na(dc)])
    thin <- findt0Smoothed(x=0:(length(dc)-1),y=dc,y.crit=threshold,above.or.below="below")
  }

  return(n/thin)
  
}

#' Finds t0 for a piecewise constant curve by linear interpolation.
#'
#' @param x The time lags at which measurements are taken (MUST INCLUDE 0!)
#' @param y The measurement associate with each x that we are using for to compute the ESS.
#' @param y.crit The critical value above (or below) which independence is achieved.
#' @param above.or.below Is t0 the time at which y >= y.crit ("above") or y <= y.crit ("below")?
#' @keywords internal
findt0Smoothed <- function(x,y,y.crit,above.or.below) {
  if ( (!x[1] == 0 && length(x) == length(y)) ) {
    stop("findt0Smoothed requires x[1] = 0 and length(x) = length(y)")
  }
  if ( above.or.below == "below" ) {
    y <- -y
    y.crit <- -y.crit
  } else if ( above.or.below != "above" ) {
    stop("Unrecognized option to findt0Smoothed.")
  }
  
  # de-duplicate entries
  n <- length(y)
  if ( (any(y[-1] == y[-n])) ) {
    first <- match(unique(y),y)
    x <- x[first]
    y <- y[first]
  }
  
  first_larger <- min(which(y >= y.crit))

  if (first_larger == 1) {
    stop("t0 cannot be negative")
  }
  
  x2 <- x[first_larger]
  x1 <- x[first_larger-1]
  y2 <- y[first_larger]
  y1 <- y[first_larger-1]
  slope <- (y2 - y1)/(x2 - x1)
  thin <- x1 + (y.crit - y1)/slope

  return(thin)
}

# # Computes ESS by finding the time lag at which the distance covariance is not distinguishable from 0
# # Uses the energy package implementation which is written in C for speed
# dCovTestESS <- function(dmat,trees=NA,min.nsamples=5,nsim=200,alpha=0.05) {
#   n <- dim(dmat)[1]
#   
#   thin <- n
#   for (i in 1:(n-min.nsamples)) {
#     p <- energy::dcov.test(dmat[-c(1:i),-c(1:i)],dmat[-c((n-i+1):n),-c((n-i+1):n)],index=1.0,R=nsim)$p.value
#     if ( p >= alpha ) {
#       thin <- i
#       break
#     }
#   }
#   return(n/thin)
# }
# 
# # Computes the MI between trees in terms of splits (trees as split coordinates)
# # splits is a matrix of all splits, trees in rows, unique splits in columns, as presence/absence
# # x_indices and y_indices refer to which rows in splits are x and which are y, may be overlapping
# treeCoordMI <- function(x_indices,y_indices,splits,summary.fn=max) {
#   # recover()
#   
#   n <- length(y_indices)
#   nsplits <- dim(splits)[2]
#   
#   p_x_1 <- colSums(splits[x_indices,])/n
#   p_y_1 <- colSums(splits[y_indices,])/n
#   
#   is_00 <- matrix(0,nrow=n,ncol=nsplits)
#   is_01 <- matrix(0,nrow=n,ncol=nsplits)
#   is_10 <- matrix(0,nrow=n,ncol=nsplits)
#   is_11 <- matrix(0,nrow=n,ncol=nsplits)
#   
#   for (i in 1:n) {
#     is_00[i,] <- as.numeric(splits[x_indices[i],] == 0 & splits[y_indices[i],] == 0)
#     is_01[i,] <- as.numeric(splits[x_indices[i],] == 0 & splits[y_indices[i],] == 1)
#     is_10[i,] <- as.numeric(splits[x_indices[i],] == 1 & splits[y_indices[i],] == 0)
#     is_11[i,] <- as.numeric(splits[x_indices[i],] == 1 & splits[y_indices[i],] == 1)
#   }
#   
#   is_00 <- colSums(is_00)
#   is_01 <- colSums(is_01)
#   is_10 <- colSums(is_10)
#   is_11 <- colSums(is_11)
#   
#   all_mi <- sapply(1:nsplits,function(i){
#     f_xy <- c(is_00[i],is_01[i],is_10[i],is_11[i])
#     f_xy <- f_xy/sum(f_xy)
#     
#     fx_fy <- c((1-p_x_1[i])*(1-p_y_1[i]),(1-p_x_1[i])*(p_y_1[i]),(p_x_1[i])*(1-p_y_1[i]),p_x_1[i]*p_y_1[i])
#     
#     # terms to sum
#     summand <- f_xy*log(f_xy/fx_fy)
#     # For MI, we may have 0 * log(0) = 0, but R produces NaN, so strip out the 0s
#     summand <- summand[!f_xy < .Machine$double.eps]
#     
#     return(sum(summand))
#   })
#   
#   return(summary.fn(all_mi))
# }

# # Computes ESS by finding the time at which MI(X_t,X_t+tau) becomes negligible
# # Interprets trees as a binary vector of split presence/absence, and computes MI between columns of this matrix
# treeAutoMIESS <- function(trees,min.nsamples=5,nsim=1000,alpha=0.05,bootstrap=TRUE,interpolate=TRUE) {
#   # recover()
#   
#   n <- length(trees)
#   
#   splits <- trees2Coords(trees)
#   
#   # Get null distribution
#   null_dist <- sapply(1:nsim,function(i){
#     permute <- sample.int(n,replace=bootstrap)
#     meanTreeCoordMI(permute[-n],permute[-1],splits)
#   })  
#   threshold <- quantile(null_dist,probs=1-alpha)
#   # cat("threshold = ",threshold,"\n")
#   
#   # get MI at time lags until we hit the threshold(s)
#   min_thresh <- min(threshold)
#   lagged <- numeric(n-min.nsamples+1)
#   for (i in 1:(n-min.nsamples+1)) {
#     # MI at this time lag
#     this_lag <- meanTreeCoordMI(1:(n-i),(i+1):n,splits)
#     # ensure distance by lag curve is monotonic
#     lagged[i] <- ifelse(i == 1,this_lag,min(this_lag,lagged[i-1]))
#     # early terimination to avoid unneeded computation
#     if (lagged[i] <= min_thresh) {
#       break
#     }
#   }
#   
#   first_smaller <- sapply(threshold,function(thresh){
#     min(which(lagged <= thresh))
#   })
#   
#   thin <- numeric(length(alpha))
#   names(thin) <- paste0(alpha*100,"%")
#   for (i in 1:length(alpha)) {
#     if (interpolate && first_smaller[i] > 1) {
#       x2 <- first_smaller[i]
#       x1 <- x2 - 1
#       y2 <- lagged[x2]
#       y1 <- lagged[x1]
#       slope <- (y2 - y1)/(x2 - x1)
#       thin[i] <- x1 + (threshold[i] - y1)/slope
#     } else {
#       thin[i] <- first_smaller[i]
#     }
#   }
#   
#   return(n/thin)
#   
# }

# # Computes ESS by finding the time lag at which the Chatterjee correlation measure (xi) is not discernably different from 0
# # The correlation is computed for the variable theta, with theta[i] = d(x[i],medioid(x))
# # If there is no unique medioid, calculates ESS for all medioids and returns the minimum
# # Uses bootstrapping to compute the null distribution of correlation coefficients
# # VERY SLOW, to be widely usable needs a C/C++ implementation of the correlation test
# chatterjeeCorrelationESS <- function(dmat,min.nsamples=5,alpha=0.05) {
#   # recover()
#   
#   n <- dim(dmat)[1]
#   
#   # there may not be one unique medioid
#   the_medioids <- which(rowSums(dmat) == min(rowSums(dmat)))
#   
#   thin <- sapply(the_medioids,function(the_medioid) {
#     lag <- n
#     for (i in 1:(n-min.nsamples)) {
#       sig <- chatterjeeCor.isSig(dmat[the_medioid,-c(1:i)],dmat[the_medioid,-c((n-i+1):n)],conf.level=1-alpha)
#       if ( !sig ) {
#         lag <- i
#         break
#       }
#     }
#     return(lag)
#   })
#   
#   
#   return(n/max(thin))
# }

# # Computes ESS by finding the time lag at which there is no discernable correlation between the adjacencies in sample iteration and the adjacencies in the sample space
# # Has much coarser resolution than others (can't go below ESS ~= min.nsamples)
# spaceTimeRankESS <- function(dmat,min.nsamples=5,alpha=0.05) {
#   
#   if ( method == "tau" ) {
#     method <- "kendall"
#   } else if ( method == "rho" ) {
#     method <- "spearman"
#   } else {
#     stop("Invalid 'method' specified")
#   }
#   
#   n <- ifelse(is.null(dim(x)),length(x),dim(x)[1])
#   
#   space_dists <- as.matrix(dist.fn(x))
#   time_dists <- as.matrix(dist(1:n))
#   
#   # To have at least min.nsamples in our testing, we need there to be at least this many samples in the distance comparison at that time lag
#   max_computable <- max(which(choose(floor(n/(1:n)),2) >= min.nsamples))
#   
#   lag <- Inf
#   for (i in 1:max_computable) {
#     
#     # get the submatrices for the distances at this time lag
#     lagged_space_dists <- lapply(1:i,function(j){
#       idx <- seq(j,n,i)
#       if ( length(idx) == floor(n/i) ) {
#         space_dists[idx,idx]
#       } else {
#         space_dists[idx[1:floor(n/i)],idx[1:floor(n/i)]]
#       }
#       
#     })
#     # we summarize the distances at each time lag by the average
#     # since we're about to take the rank, mean and sum are equivalent
#     lagged_space_dists <- Reduce("+",lagged_space_dists)
#     lagged_space_dists <- lagged_space_dists[upper.tri(lagged_space_dists)]
#     
#     lagged_time_dists <- lapply(1:i,function(j){
#       idx <- seq(j,n,i)
#       if ( length(idx) == floor(n/i) ) {
#         time_dists[idx,idx]
#       } else {
#         time_dists[idx[1:floor(n/i)],idx[1:floor(n/i)]]
#       }
#       
#     })
#     lagged_time_dists <- Reduce("+",lagged_time_dists)
#     lagged_time_dists <- lagged_time_dists[upper.tri(lagged_time_dists)]
#     
#     # there will be ties and that will make cor.test complain
#     p <- suppressWarnings(cor.test(lagged_space_dists,lagged_time_dists,method="spearman")$p.value)
#     
#     if (p >= alpha) {
#       lag <- i
#       break
#     }
#   }
#   
#   if (is.finite(lag)) {
#     return(n/lag)
#   } else {
#     return(paste0("<",min.nsamples))
#   }
# }
