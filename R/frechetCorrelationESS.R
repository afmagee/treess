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

# Very slow, but seems to be working (along with all necessary functions it rests upon)
# Cannot do anything but the lower bound here
subsampledFrechetCorrelationESS <- function(x,dist.fn,n.per.diag,min.nsamples=5,regularize.average.weight=1,lower.bound=TRUE,...) {
  # recover()
  
  n <- length(x)

  n_per_diag <- n.per.diag # this should be a function parameter and will need to be experimented with
  wt <- regularize.average.weight # weight of samples outside this box, this should get fed to the appropriate internal function
  
  # Strategy: subsample
  #  1. think in bands along the diagonals
  #  1.1 question of whether to spread points uniformly or to concentrate towards the topleft/bottomright where most "variance boxes" live
  #  2. store matrix with i,j and distance (squared)
  #  2.2 access kth diagonal with i == j + k
  #  3. between E[d^2] bit is easy, pre-compute because will need for variances
  #  4. variances within sub-chunks is a bit trickier
  #  4.1 requires knowing which off-diagonals a box uses
  #  4.2 once we know that, we can access the correct diagonal as (i == j + k) & i >= thresh
  #  4.3 some diagonals will be empty or sparse, probably want a weighted average of points within the box and the whole diagonal
  #  4.4 weighting scheme could treat the "global" average for the diagonal as c cells (try c = 1 to start)
  
  i_j_dist <- computeDiagonallySubsampledDistances(x,dist.fn,n.per.diag)
  i_j_dist <- i_j_dist[,-4]
  i_j_dist[,4] <- i_j_dist[,4]^2
  
  
  d12 <- sapply(1:(n-1),function(i){
    mean(i_j_dist[i_j_dist[,3] == i + 1,4])
  })
  
  # recover()
  
  # Compute covariances only as far as we need to
  cors <- numeric(n-min.nsamples-1)
  P <- rep(NA,floor((n-min.nsamples)/2))
  for (i in 1:(n-min.nsamples-1)) {
    n_on_diagonal <- n - i
    
    var1 <- estimateSumFromDiagonalSubsamples(i_j_dist,
                                              diag.means=d12,
                                              first=i+1,
                                              last=n,
                                              weight.per.point=1,
                                              weight.unsampled.average=wt)/(2*(n-i)*(n-i-1)) # extra factor of 2 because we're summing over the whole square and not an upper/lower triangular portion
    var2 <- estimateSumFromDiagonalSubsamples(i_j_dist,
                                              diag.means=d12,
                                              first=1,
                                              last=n-i,
                                              weight.per.point=1,
                                              weight.unsampled.average=wt)/(2*(n-i)*(n-i-1))

    # lower bound on 2 x covariance
    covar <- var1 + var2 - d12[i]
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

computeDiagonallySubsampledDistances <- function(x,dist.fn,n.per.diag) {
  n <- length(x)
  
  # recover()
  
  # A matrix containing all the pairs for which we will compute distances as (d,k) coordinates along the diagonal
  dk <- NULL
  
  if (n.per.diag < length(x) - 1) {
    
    # # number of distances to compute: number from diagonals we subsample + number from those we fully sample
    # n_compute <- n.per.diag * n
    # 
    # # ensure at least one sample in every diagonal band but weight sampling towards longer bands
    # diag_samples <- sample(2:(n-1),size=n_compute - (n - 1),replace=TRUE,prob=((n-2):1))
    # n_samples <- 1 + sapply(2:n,function(d){
    #   sum(diag_samples == d)
    # })
    # 
    # dk <- matrix(nrow=n_compute,ncol=2)
    # idx <- 0
    # for (d in 2:n) {
    #   n_on_diagonal <- n - d + 1
    #   subseq <- round(seq(1,n_on_diagonal,length.out=n_samples[d-1]))
    #   indices <- idx + (1:n_on_diagonal)
    #   dk[indices,1] <- rep(d,n_on_diagonal)
    #   dk[indices,2] <- subseq
    #   idx <- idx + n_on_diagonal
    # }
    # 
    # idx <- n.per.diag * (last_subsampled_diagonal - 1)
    # for (d in (last_subsampled_diagonal+1):n) {
    #   n_on_diagonal <- n - d + 1
    #   dk[idx+(1:n_on_diagonal),1] <- rep(d,n_on_diagonal)
    #   dk[idx+(1:n_on_diagonal),2] <- 1:n_on_diagonal
    #   idx <- idx + n_on_diagonal
    # }
    
    # Uniformly over each diagonal
    # Find the diagonal beyond which we compute all distances, not a subsampled set
    last_subsampled_diagonal <- n - n.per.diag

    # number of distances to compute: number from diagonals we subsample + number from those we fully sample
    n_compute <- n.per.diag * (last_subsampled_diagonal - 1) + sum(n.per.diag:1)

    dk <- matrix(nrow=n_compute,ncol=2)
    for (d in 2:last_subsampled_diagonal) {
      n_on_diagonal <- n - d + 1
      subseq <- round(seq(1,n_on_diagonal,length.out=n.per.diag))
      dk[((d-2)*n.per.diag+1):((d-1)*n.per.diag),1] <- rep(d,n.per.diag)
      dk[((d-2)*n.per.diag+1):((d-1)*n.per.diag),2] <- subseq
    }

    idx <- n.per.diag * (last_subsampled_diagonal - 1)
    for (d in (last_subsampled_diagonal+1):n) {
      n_on_diagonal <- n - d + 1
      dk[idx+(1:n_on_diagonal),1] <- rep(d,n_on_diagonal)
      dk[idx+(1:n_on_diagonal),2] <- 1:n_on_diagonal
      idx <- idx + n_on_diagonal
    }
  } else {
    n_compute <- choose(n,2)
    dk <- matrix(nrow=n_compute,ncol=2)
    idx <- 0
    for (d in 2:n) {
      n_on_diagonal <- n - d + 1
      dk[idx+(1:n_on_diagonal),1] <- rep(d,n_on_diagonal)
      dk[idx+(1:n_on_diagonal),2] <- 1:n_on_diagonal
      idx <- idx + n_on_diagonal
    }
  }
  
  # (i,j) pairs for which we will compute distances
  ij <- dk2ij(dk[,1],dk[,2])
  
  dists <- getSubsampledDistanceMatrix(x,dist.fn,ij)
  
  return(cbind(ij,dk,dists))
  
}

# Finds the index of the kth element on the dth diagonal
# The first diagonal is the actual diagonal, d such that j = i
# Diagonals increase into the upper right corner, the second is d such that j = i + 1, the last is (1,n)
dk2ij <- function(d,k) {
    # c(k,k+d-1)
  # apply(cbind(d,k),1,function(dk){
  #   c(dk[2], dk[2] + dk[1] - 1)
  # })
  unname(cbind(k,k+d-1))
}

# The reverse mapping of dk2ij
ij2dk <- function(i,j) {
  # c(j-i+1,i)
  # apply(cbind(i,j),1,function(ij){
  #   c(ij[2] - ij[1] + 1, ij[1])
  # })
  unname(cbind(j-i+1,i))
}

estimateSumFromDiagonalSubsamples <- function(x,
                                              diag.means,
                                              first,
                                              last,
                                              weight.per.point=1,
                                              weight.unsampled.average=1) {
  i_j_dist <- x

  # Size of this submatrix
  m <- last - first + 1
  
  # recover()
  is_in_box <- i_j_dist[,1] >= first & i_j_dist[,1] <= last & i_j_dist[,2] >= first & i_j_dist[,2] <= last
  i_j_dist <- i_j_dist[is_in_box,]
  summand <- sapply(2:m,function(d){
    n_on_diag <- m - (d - 1)
    to_use <- i_j_dist[,3] == d
    n_in_box <- sum(to_use)
    if (n_in_box > 0) {
      sample_mean <- mean(i_j_dist[to_use,4])
    } else {
      sample_mean <- 0
    }
    w_in <- n_in_box * weight.per.point
    w_out <- weight.unsampled.average
    w_tot <- w_in + w_out
    n_on_diag * (sample_mean * w_in/w_tot + diag.means[d-1] * w_out/w_tot)
  })
  return(2 * sum(summand))
}

# # Cannot do anything but the lower bound here
# subsampledFrechetCorrelationESS <- function(x,dist.fn,n.per.diag,min.nsamples=5,lower.bound=TRUE,rooted=FALSE,...) {
#   recover()
#   
#   n <- length(x)
#   if ( class(x) == "matrix" ) {
#     n <- dim(x)[1]
#   }
#   
#   n_per_diag <- n.per.diag # this should be a function parameter and will need to be experimented with
#   wt <- 1 # weight of samples outside this box
#   
#   # Strategy: subsample
#   #  1. think in bands along the diagonals
#   #  1.1 question of whether to spread points uniformly or to concentrate towards the topleft/bottomright where most "variance boxes" live
#   #  2. store matrix with i,j and distance (squared)
#   #  2.2 access kth diagonal with i == j + k
#   #  3. between E[d^2] bit is easy, pre-compute because will need for variances
#   #  4. variances within sub-chunks is a bit trickier
#   #  4.1 requires knowing which off-diagonals a box uses
#   #  4.2 once we know that, we can access the correct diagonal as (i == j + k) & i >= thresh
#   #  4.3 some diagonals will be empty or sparse, probably want a weighted average of points within the box and the whole diagonal
#   #  4.4 weighting scheme could treat the "global" average for the diagonal as c cells (try c = 1 to start)
#   
#   # Compute a subsample of distances
#   n_compute <- n_per_diag*(n-n_per_diag-1) + sum(n_per_diag:1) # number of diagonals we subsample + number we fully sample
#   i_j_dist <- matrix(nrow=n_compute,ncol=4)
#   i_j_dist_contributes_to1 <- vector("list",n_compute)
#   i_j_dist_contributes_to2 <- vector("list",n_compute)
#   
#   # Loop over diagonals
#   # The overhead per function call is staggering, lots of unneeded error checking and maybe excess RCPP overhead from all the calls
#   # We need a way around calling RF.dist
#   idx <- 1
#   for (k in 1:(n-1)) {
#     n_on_diag <- n - k
#     index_on_diag <- 1:n_on_diag
#     if (n_on_diag > n_per_diag) {
#       index_on_diag <- round(seq(1,n_on_diag,length.out=n_per_diag))
#     }
#     for (l in index_on_diag) {
#       i <- l
#       j <- k + l
#       i_j_dist[idx,1] <- i
#       i_j_dist[idx,2] <- j
#       i_j_dist[idx,3] <- k
#       # Not working yet exactly, definitely should only remove the <=0 indices when we know they happen
#       i_j_dist_contributes_to1[[idx]] <- (i+1):n
#       i_j_dist_contributes_to1[[idx]] <- i_j_dist_contributes_to1[[idx]][i_j_dist_contributes_to1[[idx]] > 0]
#       i_j_dist_contributes_to2[[idx]] <- 1:(n-j-1)
#       i_j_dist_contributes_to2[[idx]] <- i_j_dist_contributes_to2[[idx]][i_j_dist_contributes_to2[[idx]] > 0]
#       idx <- idx + 1
#     }
#   }
#   
#   i_j_dist[,4] <- select.RF.dist(trees,i_j_dist[,1:2],FALSE,rooted=rooted)^2
#   
#   # To avoid repeated computation, we loop over the off-diagonals and assign the values where they're needed
#   vars1 <- numeric(n-1)
#   nv1 <- numeric(n-1)
#   vars2 <- numeric(n-1)
#   nv2 <- numeric(n-1)
#   d_sq <- numeric(n-1)
#   
#   # We should be looping over the cells of i_j_dist directly and adding the squared values into d_sq, vars1/2
#   # Then, we can loop over and divide by the number of contributions
#   # So, loop 1: compute d_sq and add to appropriate vars 1/2
#   # Loop 2, use i_j_dist_contributes_to1/2 to know which vars 1/2 to add this value to
#   # Loop 3, divide
#   #### This isn't quite right, we need to take the weighted average for the values on that diagonal, then add in that times the number on that diagonal!
#   for (i in 1:(n-2)) {
#     idx_diag <- which(i_j_dist[,3] == i)
#     d_sq[i] <- mean(i_j_dist[idx_diag,4])
#     vars1[1:(n-1-i)] <- vars1[1:(n-1-i)] + d_sq[i]
#     vars2[1:(n-1-i)] <- vars2[1:(n-1-i)] + d_sq[i]
#     nv1[1:(n-1-i)] <- nv1[1:(n-1-i)] + wt
#     nv2[1:(n-1-i)] <- nv2[1:(n-1-i)] + wt
#   }
#   
#   for (i in 1:dim(i_j_dist)[1]) {
#     vars1[i_j_dist_contributes_to1[[i]]] <- vars1[i_j_dist_contributes_to1[[i]]] + i_j_dist[i,4]
#     vars2[i_j_dist_contributes_to2[[i]]] <- vars2[i_j_dist_contributes_to2[[i]]] + i_j_dist[i,4]
#     nv1[i_j_dist_contributes_to1[[i]]] <- nv1[i_j_dist_contributes_to1[[i]]] + 1
#     nv2[i_j_dist_contributes_to2[[i]]] <- nv2[i_j_dist_contributes_to2[[i]]] + 1
#   }
#   
#   for (i in 1:(n-1)) {
#     vars1[i] <- vars1[i]/nv1[i]
#     vars2[i] <- vars2[i]/nv2[i]
#   }
#   
#   vars1 <- vars1 * 2
#   vars2 <- vars2 * 2
#   
#   cors <- (vars1 + vars2 - d_sq)/sqrt(vars1*vars2)
#   
#   
#   cors <- numeric(n-1)
#   P <- rep(NA,floor((n-1)/2))
#   
#   # Smoothed P, aka P' in Vehtari et al.
#   # Remove any P we did not compute and thus do not need
#   P <- P[!is.na(P)]
#   for (i in 2:length(P)) {
#     P[i] <- min(P[i],P[i-1])
#   }
#   
#   # Unless we summed over all time lags, we stopped when P[length(P)] < 0
#   k <- length(P) - 1
#   if ( P[length(P)] > 0 ) {
#     k <- length(P)
#   }
#   tau_hat <- -1 + 2 * sum(P[1:k])
#   
#   # Paranoid exception handling
#   if (tau_hat < 0) {
#     tau_hat <- 1 
#   }
#   
#   return(n/tau_hat)
#   
# }

