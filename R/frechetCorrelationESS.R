#' Calculates ESS using a Frechet-like generalization of a univariate ESS based on the autocorrelation of the chain at varying timesteps.
#'
#' @param dmat Sample-to-sample distance matrix for the MCMC chain.
#' @param trees For compatibility with eval and call, not used here.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples We only compute the correlation between x_t and x_{t+T} if there are at least this many samples.
#' @param lower.bound Should we use the lower-bound ESS or attempt to estimated E(x_t) - E(x_{t+T})^2? Lower bound is faster.
#' @keywords internal
frechetCorrelationESS <- function(dmat,trees=NA,nsim=NA,alpha=NA,min.nsamples=5,lower.bound=TRUE) {
  # recover()
  
  # If there's only 1 tree sampled, the ESS is 1
  if ( sum(dmat) == 0 ) {
    return(1)
  }
  
  n <- dim(dmat)[1]
  
  dmat <- dmat^2
  
  cors <- numeric(n-min.nsamples-1)
  for (i in 1:(n-min.nsamples-1)) {
    rs1 <- rowSums(dmat[-c(1:i),-c(1:i)])
    rs2 <- rowSums(dmat[-c((n-i+1):n),-c((n-i+1):n)])
    var1 <- sum(rs1)/choose(n-i,2)/4
    var2 <- sum(rs2)/choose(n-i,2)/4
    
    d12 <- mean(dmat[row(dmat) == (col(dmat)+i)])
    
    # lower bound on 2 x covariance
    covar <- var1 + var2 - d12
    
    if ( !lower.bound ) {
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
  }
  
  P <- numeric(floor((n-min.nsamples)/2))
  P[1] <- cors[1] + cors[2]
  for (i in 2:floor((n-min.nsamples)/2)) {
    P[i] <- min(P[i-1],cors[2*i-1] + cors[2*i])
  }
  
  tau_hat <- 1
  if ( P[1] > 0 ) {
    k <- max(which(P > 0))
    tau_hat <- 1 + 2 * sum(cors[1:(2*k+1)])
  }
  
  return(n/tau_hat)
  
}
