#' Calculates the approximate ESS of Lanfear et al. (2015).
#'
#' @param dmat Sample-to-sample distance matrix for chains, NOT SQUARED.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha The cutoff proportion of the asymptoted used to estimate the average pairwise squared distance.
#' @param min.nsamples We only compute the average squared distance between x_t and x_{t+T} if there are at least this many samples.
#' @keywords internal
approximateESS <- function(dmat,trees=NA,nsim=NA,min.nsamples=5,alpha=0.05) {
  # recover()
  
  dmat <- dmat^2
  
  n <- dim(dmat)[1]
  
  # alpha to threshold
  thresh = 1 - alpha
  
  # turn minimum number of samples into the maximum off-diagonal
  max_off_diag <- n - min.nsamples
  
  # get empirical curve of lag time t vs distance
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

  E_delta_squared_iid <- 0
  
  if ( max(d_t) < thresh ) {
    E_delta_squared_iid <- max(dmat)
  } else {
    thin <- min(which(d_t >= thresh))
    E_delta_squared_iid <- sum(unlist(lapply(thin:n,function(i){
      dmat[row(dmat) == col(dmat) + i]
    })))/choose(n-thin,2)
  }
  
  # Solve for ESS 
  rhs <- sum(dmat)/(n^2)
  
  ess <- 1/(1 - rhs/E_delta_squared_iid)
  
  return(ess)
}
