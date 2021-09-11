#' Calculates ESS using a Frechet-like generalization of a univariate ESS based on the autocorrelation of the chain at varying timesteps.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples We only compute the correlation between x_t and x_{t+T} if there are at least this many samples.
#' @param lower.bound Should we use the lower-bound ESS or attempt to estimated E(x_t) - E(x_{t+T})^2? Lower bound is faster.
#' @keywords internal
frechetCorrelationESS <- function(dmat,min.nsamples=5,lower.bound=TRUE,...) {
  # recover()
  
  # If there's only 1 tree sampled, the ESS is 1
  if ( sum(dmat) == 0 ) {
    return(1)
  }
  
  n <- dim(dmat)[1]
  
  dmat <- dmat^2
  
  # Compute covariances only as far as we need to
  cors <- numeric(n-min.nsamples-1)
  P <- rep(NA,floor((n-min.nsamples)/2))
  # # TODO: the following should be a more efficient way to get var1 and var2, at the cost of only working for a lower bound
  # ssq_front_back <- diag(apply(apply(dmat, 2, cumsum), 1, cumsum))
  # ssq_back_front <- diag(apply(apply(dmat[n:1,n:1], 2, cumsum), 1, cumsum))
  # for (i in 1:(n-min.nsamples-1)) {
  #   var1 <- ssq_front_back[n-i]/(2*(n-i)*(n-i-1))
  #   var2 <- ssq_back_front[n-i]/(2*(n-i)*(n-i-1))
  for (i in 1:(n-min.nsamples-1)) {
    rs1 <- rowSums(dmat[-c(1:i),-c(1:i)])
    rs2 <- rowSums(dmat[-c((n-i+1):n),-c((n-i+1):n)])
    var1 <- sum(rs1)/(2*(n-i)*(n-i-1)) # extra factor of 2 because we're summing over the whole square and not an upper/lower triangular portion
    var2 <- sum(rs2)/(2*(n-i)*(n-i-1))
    
    d12 <- mean(dmat[row(dmat) == (col(dmat)+i)])
    
    # lower bound on 2 x covariance
    covar <- var1 + var2 - d12
    
    if ( !lower.bound ) {
      stop("This code is untested/unproven")
      # Find index of Karcher mean tree(s), may not be unique
      mu1 <- i + which(rs1 == min(rs1))
      mu2 <- which(rs2 == min(rs2))
      
      # Estimate (mu1 - mu2)^2
      mu1_mu2_sq <- mean(unlist(lapply(mu1,function(m1){
        lapply(mu2,function(m2){
          dmat[m1,m2]
        })
      })))
      
      # non-lower-bound estimate of 2 x covar
      covar <- covar + mu1_mu2_sq
    }
    
    covar <- covar/2
    
    # If either of the variances are 0, then one set of trees is all the same and the correlation is technically undefined
    # Practically, this means there is very little variability in these trees, the sampler is mixing poorly, and we can report this with a high correlation
    if ( var1 == 0 || var2 == 0 ) {
      cors[i] <- 1
    } else {
      cors[i] <- covar/sqrt(var1*var2) 
    }
    
    # We want to combine (0,1), (2,3), ... but we're indexing starting at 1
    if ( i %% 2 == 1 ) {
      if ( i > 1) {
        P[(i+1)/2] <- cors[i] + cors[i-1]
      } else {
        # cors[0] is 1.0
        P[(i+1)/2] <- cors[i] + 1.0
      }
      if ( P[(i+1)/2] < 0 ) {
        break
        # We only sum over P > 0 so we don't need further terms
      }
    } 
  }
  
  # Smoothed P, aka P' in Vehtari et al.
  # Remove any P we did not compute and thus do not need
  P <- P[!is.na(P)]
  for (i in 2:length(P)) {
    P[i] <- min(P[i],P[i-1])
  }
  
  # Unless we summed over all time lags, we stopped when P[length(P)] < 0
  k <- length(P) - 1
  if ( P[length(P)] > 0 ) {
    k <- length(P)
  }
  tau_hat <- -1 + 2 * sum(P[1:k])

  # Paranoid exception handling
  if (tau_hat < 0) {
    tau_hat <- 1 
  }
  
  return(n/tau_hat)
  
}
