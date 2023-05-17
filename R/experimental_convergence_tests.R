#' Lists all available methods for computing test of convergence of multiple MCMC chains of trees/multivariate objects.
#' 
#' @keywords internal
#' @examples
#' getConvergenceTestMethods()
getConvergenceTestMethods <- function() {
  all_methods <- c("pseudoPSRF",
                   "totalDistancePSRF",
                   "foldedRankMedioidPSRF",
                   "CMDSPSRF",
                   "medianDistanceWilcoxTest",
                   "ASDSFPermutationTest")
  return(sort(all_methods))
}

#' Test for global convergence of multiple MCMC chains.
#'
#' @param x List of all chains of the tree-valued (or multivariate or non-Euclidean) variable(s) for which to compute convergence.
#' @param dist.fn A function suitable to calculate distances between all values in x.
#' @param method The method for computing the convergence diagnostic (see details).
#' @param alpha For methods using  (only matters when simplify.output == TRUE).
#' @param threshold A maximum threshold for PSRF to de
#' @param identify.outliers If TRUE, and global convergence is not found (for any method), identifies which chain(s) are different from the others. Only works if there are at least 3 chains
#' @param min.split.freq For ASDSF-based tests only. Minimum split frequency to be included in ASDSF.
#' @details Allowable methods are "pseudoPSRF", "totalDistancePSRF", "foldedRankMedioidPSRF", and "medianDistanceWilcoxTest"
#' @return A list, where $converged is a logical vector indicating whether global convergence was found for each method used, $statistic is the computed PSRF or test statistic for each method, and (optionally) $pairwise is a list of matrices with the results of all pairs of tests for all methods.
#' @keywords internal
#' @examples
#' testConvergence(list(rnorm(100),rnorm(100)))
#' testConvergence(list(rnorm(100,0,1),list(rnorm(100,0,1),list(rnorm(100,0,1),rnorm(100,1,1)))
testConvergence <- function(x,dist.fn,methods=getConvergenceTestMethods(),alpha=0.05,threshold=1.01,identify.outliers=TRUE,min.split.freq=0.01) {
  # recover()
  
  # Check for valid inputs
  if ( !(any(c("list","mcmc.list") %in% class(x))) ) {
    stop("Argument 'x' must be a list.")
  }
  
  is_tree <- "multiPhylo" %in% class(x[[1]])
  is_matrix <- "maxtrix" %in% class(x[[1]])
  
  nchains <- length(x)
  
  # ensure all chains are of same dimension(s)
  lens <- unlist(lapply(x,length))
    
  if ( length(unique(lens)) != 1 ) {
    stop("All elements of `x` must be of same length/dimension")
  }
  
  if ( is_matrix ) {
    ngen <- dim(x[[1]])[1]
  } else {
    ngen <- lens[1]
  }
  
  # Make sure the chosen method(s) are valid
  if ( any(!(methods %in% getConvergenceTestMethods())) ) {
    stop("Invalid option to argument 'methods'.")
  }
  
  dmat <- matrix()
  if ( any(c("pseudoPSRF","totalDistancePSRF","foldedRankMedioidPSRF","CMDSPSRF","medianDistanceWilcoxTest") %in% methods) ) {
    # Compute distance matrix, convergence diagnostics
    if ( is_tree ) {
      x <- do.call(c,x)
      class(x) <- "multiPhylo"
    } else if ( is_matrix ) {
      x <- do.call(rbind,x)
    } else {
      x <- do.call(c,x)
    }
    
    dists <- dist.fn(x)
    
    if ( !("dist" %in% class(dists)) ) {
      stop("dist.fn must return object of class dist")
    }
    
    dmat <- as.matrix(dists)
  }
  
  coords <- matrix()
  if ( "ASDSFPermutationTest" %in% methods ) {
    if ( !("multiPhylo" %in% class(x[[1]])) ) {
      stop("Cannot compute ASDSFPermutationTest for non-tree-valued MCMC chains.")
    }
    coords <- trees2Coords(do.call(c,x))
  }
  
  
  convergence <- sapply(methods,function(psrf_method){
    res <- eval(call(psrf_method,dmat=dmat,ngen=ngen,nchains=nchains,coords=coords,min.split.freq=min.split.freq))
    if ( grepl("test",tolower(psrf_method)) ) {
      return(c(res$p.value >= alpha,res$statistic))
    } else {
      return(c(res$psrf[1,1] < threshold,res$psrf[1,1]))
    }
  })
  
  converged <- as.logical(convergence[1,])
  names(converged) <- methods
  
  statistic <- convergence[2,]
  
  res <- list(converged=converged,statistic=statistic)
  
  if ( identify.outliers && any(!converged) ) {
    pairwise <- lapply(methods,function(psrf_method){
      
      all_comparisons <- matrix(TRUE,nchains,nchains)
      
      for (i in 1:(nchains-1)) {
        for (j in (i+1):nchains) {
          # construct the submatrix containing only these 2 chains
          two_by_two <- matrix(NA,2*ngen,2*ngen)
          two_by_two[1:ngen,1:ngen] <- dmat[((i-1)*ngen+1):(i*ngen),((i-1)*ngen+1):(i*ngen)]
          two_by_two[(ngen+1):(2*ngen),(ngen+1):(2*ngen)] <- dmat[((j-1)*ngen+1):(j*ngen),((j-1)*ngen+1):(j*ngen)]
          two_by_two[1:ngen,(ngen+1):(2*ngen)] <- dmat[((i-1)*ngen+1):(i*ngen),((j-1)*ngen+1):(j*ngen)]
          two_by_two[(ngen+1):(2*ngen),1:ngen] <- dmat[((j-1)*ngen+1):(j*ngen),((i-1)*ngen+1):(i*ngen)]
          
          # test convergence and record
          res <- eval(call(psrf_method,dmat=two_by_two,ngen=ngen,nchains=2,coords=coords,min.split.freq=min.split.freq))
          passes <- ifelse(grepl("test",tolower(psrf_method)),res$p.value >= alpha,res$psrf < threshold)
          all_comparisons[i,j] <- passes
          all_comparisons[j,i] <- all_comparisons[i,j]
        }
      }
      return(all_comparisons)
    })
    names(pairwise) <- methods
    res$pairwise <- pairwise
  }
  
  return(res)
}

#' Tests the null hypothesis that a set of chains come from the same distributions against the alternative that they are different.
#'
#' @param coords A 0/1 matrix of the splits in each tree (RF coordinate matrix).
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @param nrep Number of permutation replicates.
#' @keywords internal
ASDSFPermutationTest <- function(coords,ngen,nchains,min.split.freq,nrep=500,...) {
  # recover()
  
  # observed ASDSF
  chain_labels <- sort(rep(1:nchains,ngen))
  observed_probs <- lapply(1:nchains,function(i){
    colMeans(coords[chain_labels == i,])
  })
  observed_asdsf <- ASDSF(observed_probs,min.freq=min.split.freq)

  # permutation test
  null_dist <- sapply(1:nrep,function(i){
    permuted_labels <- sample(chain_labels)
    probs <- lapply(1:nchains,function(i){
      colMeans(coords[permuted_labels == i,])
    })
    asdsf <- ASDSF(probs,min.freq=min.split.freq)
    return(asdsf)
  })

  p <- (sum(null_dist >= observed_asdsf)+0.5)/(nrep+1)

  return(list(statistic=observed_asdsf,p.value=p))

}

#' Calculates the PSRF in the spirit of the pseudo-ESS of Lanfear et al. (2015)
#' Returns the maximum (over all reference trees in the chain) of all PSRF values.
#'
#' @param dmat Sample-to-sample distance matrix for ALL chains.
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @keywords internal
pseudoPSRF <- function(dmat,ngen,nchains,...) {
  # recover()
  
  all_psrf <- sapply(1:dim(dmat)[1],function(idx) {
    X <- dmat[,idx]
    X_mcmc <- lapply(1:nchains,function(i){
      coda::as.mcmc(X[((i-1)*ngen+1):(i*ngen)])
    })
    as.numeric(coda::gelman.diag(X_mcmc)$psrf[1,1])
  })
  
  max_index <- which.max(all_psrf)
  
  return(all_psrf[[max_index]])
}

#' Calculates the PSRF in the spirit of the totalDistanceESS
#'
#' @param dmat Sample-to-sample distance matrix for ALL chains.
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @keywords internal
totalDistancePSRF <- function(dmat,ngen,nchains,...) {
  # recover()
  
  X <- rowSums(dmat)
  
  X_mcmc <- lapply(1:nchains,function(i){
    coda::as.mcmc(X[((i-1)*ngen+1):(i*ngen)])
  })
  
  return(coda::gelman.diag(X_mcmc))
}

#' Calculates an PSRF in the spirit of the folded rank-transformed PSRF of Vehtari et al. (2019), replacing the median with the medioid.
#' Where there is no unique medioid, calculates the PSRF for all and returns the maximum.
#'
#' @param dmat Sample-to-sample distance matrix for ALL chains.
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @keywords internal
foldedRankMedioidPSRF <- function(dmat,ngen,nchains,...) {
  
  # there may not be one unique medioid
  the_medioids <- which(rowSums(dmat) == min(rowSums(dmat)))
  
  all_psrf <- lapply(the_medioids,function(the_medioid) {
    xr <- rank(dmat[,the_medioid])
    X <- qnorm((xr+0.5)/(max(xr)+1))
    X_mcmc <- lapply(1:nchains,function(i){
      coda::as.mcmc(X[((i-1)*ngen+1):(i*ngen)])
    })
    coda::gelman.diag(X_mcmc)
  })
  
  # find maximum of PSRF point estimates, return that PSRF object
  point_ests <- unlist(lapply(all_psrf,function(psrf){psrf$psrf[1,1]}))
  max_index <- which.max(point_ests)
  
  return(all_psrf[[max_index]])
}

#' Calculates PSRF for a single summary variable obtained from classical multidimensional scaling
#'
#' @param dmat Sample-to-sample distance matrix for ALL chains.
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @keywords internal
CMDSPSRF <- function(dmat,ngen,nchains,...) {
  X <- stats::cmdscale(dmat,k=1)[,1]

  X_mcmc <- lapply(1:nchains,function(i){
    coda::as.mcmc(X[((i-1)*ngen+1):(i*ngen)])
  })
  
  return(coda::gelman.diag(X_mcmc))
}

#' Computes a Wilcoxon signed rank test on a pair of transformed variables B and W
#' The original sample is an array X and x_ij is the ith sample of the variable from chain j (x_ij can be a vector, a tree, or any object for which we can compute distances)
#' Then W_i = median(d(x_ij,x_.j)) and B_i = median(d(x_ij,x_i.))
#' That is, W is the vector of within-chain median distances (median distance from each point to every other point in the same chain), and B is the vector of between chain median distances
#' If all chains are sampling the same distribution, then W and B should come from the same distribution
#' Calculates PSRF for a single summary variable obtained from classical multidimensional scaling
#'
#' @param dmat Sample-to-sample distance matrix for ALL chains.
#' @param ngen Number of generations in each chain (cannot mix chain sizes).
#' @param nchains Number of chains that are in the distance matrix.
#' @keywords internal
medianDistanceWilcoxTest <- function(dmat,ngen,nchains,...) {
  # recover()
  
  B <- numeric(ngen*nchains)
  W <- numeric(ngen*nchains)
  idx <- 0
  for (i in 1:nchains) {
    within_indices <- c((1+(i-1)*ngen):(i*ngen))
    for (j in 1:ngen) {
      idx <- idx + 1
      # within-chain comparison 
      W[idx] <- median(dmat[idx,within_indices])
      # within-chain comparison 
      B[idx] <- median(dmat[idx,-within_indices])
    }
  }
  
  return(wilcox.test(x=B,y=W,paired=TRUE))
}
