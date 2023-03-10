#' Calculates the approximate ESS of Lanfear et al. (2015).
#'
#' @param dmat Sample-to-sample distance matrix for chains, NOT SQUARED.
#' @param max.approximateESS.timelag The maximum time lag for computing jump distances
#' @param alpha The cutoff proportion of the asymptoted used to estimate the average pairwise squared distance.
#' @keywords internal
approximateESS <- function(dmat,max.approximateESS.timelag=100,alpha=0.05,...) {
  # This function essentially strings together the following RWTY functions
  # topological.autocorr, topological.approx.ess, approx.ess.multi, approx.ess.single, tree.autocorr

  dmat <- dmat^2
  
  N <- n <- dim(dmat)[1]
  
  # alpha to threshold
  thresh = 1 - alpha
  
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
  
  # Optimize the exponential curve fit
  fn <- function(par,t,d_t) {
    E_d <- par[1] * (1 - exp(-t/par[2]))
    return( sum((E_d - d_t)^2) )
  }
  par <- optim(c(1,1),fn,t=t,d_t=d_t)
    
  # For finding the cutoff lag
  thresh <- thresh * par$par[1]
  
  m <- ifelse( max(d_t) < thresh, length(d_t)+1, min(which(d_t >= thresh)))
  
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
