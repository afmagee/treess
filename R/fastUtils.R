#' Get all distances for pairs (x[[i]],x[[j]]) of objects specified by the (k x 2) matrix ij
#' 
#' @param x The values
#' @param dist.fn A function which can be called dist.fn(x[[i]],x[[j]])
#' @param ij A 2-column matrix specifying (i,j) pairs for which the distances will be computed, or a list of such matrices
#' @return Distances in order specified by matrix ij, either as a vector (class(ij) == "matrix") or a list (class(ij) == "list")
getSubsampledDistanceMatrix <- function(x,dist.fn,ij) {
  res <- NULL
  if ( "list" %in% class(ij) ) {
    res <- lapply(ij,function(ij_){
      sapply(1:dim(ij_)[1],function(idx){
        dist.fn(x[[ij_[idx,1]]],x[[ij_[idx,2]]])
      })
    })
  } else if ( "matrix" %in% class(ij) ) {
    res <- sapply(1:dim(ij)[1],function(idx){
      dist.fn(x[[ij[idx,1]]],x[[ij[idx,2]]])
    })
  } else {
    stop("Invalid input format for argument \"ij\"")
  }
  return(res)
}

# Get a tree ready for RF distance calculations, or for splitFrequencyESS
preprocessForRF <- function(trees,rooted) {
  treeSplits(trees,rooted=rooted)
}

# Get a tree ready for path distance calculation
preprocessForPD <- function(trees,rooted,use.edge.length=FALSE) {
  ntaxa <- length(trees[[1]]$tip.label)
  trees <- checkRootedOption(trees,rooted)
  trees <- ape::.compressTipLabel(trees)

  # This would be faster with phangorn:::coph, but for now playing it safe using only exported functions
  ut <- upper.tri(matrix(NA,ntaxa,ntaxa))
  tip_to_tip_dists <- lapply(trees,function(phy){
    if ( !use.edge.length ) {
      phy$edge.length <- rep(1.0,length(phy$edge.length))
    }
    dm <- cophenetic.phylo(phy)
    return(dm[ut])
  })
  return(tip_to_tip_dists)
}

# Compute the RF distance from a set of splits
# fastmatch::fmatch is faster than the builtin matching and speeds this up
rfDistFromSplits <- function(splits1,splits2) {
  ntot <- length(splits1) + length(splits2)
  nmatch <- sum(!is.na(fastmatch::fmatch(splits1,splits2)))
  return(ntot - 2 * nmatch)
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

#' Compute subsample of distances from a distance matrix.
#' 
#' @param x The values
#' @param dist.fn A function which can be called dist.fn(x[[i]],x[[j]])
#' @param n.per.diag The number of samples to take on each diagonal band
#' @param uniform.over.matrix If FALSE, each diagonal is sampled equally, if TRUE samples are dispersed over the matrix more uniformly
#' @return List of matrices, for the 2nd through nth diagonal band (named accordingly), each a matrix. Matrix columns are row, column, distance.
#' @keywords internal
computeDiagonallySubsampledDistances <- function(x,dist.fn,n.per.diag,uniform.over.matrix=FALSE) {
  n <- length(x)
  
  # recover()
  
  n_compute_per_diag <- NULL
  if (n.per.diag < length(x) - 1) {
    if ( !uniform.over.matrix ) {
      # Uniformly over each diagonal
      # Find the diagonal beyond which we compute all distances, not a subsampled set
      last_subsampled_diagonal <- n - n.per.diag
      
      # number of distances to compute on each diagonal
      n_compute_per_diag <- c(rep(n.per.diag,(last_subsampled_diagonal-1)),n.per.diag:1)
    } else {
      n_compute <- n.per.diag * n
      
      w <- (n-2):0
      w <- w/sum(w)
      w <- round(n.per.diag * (n - 1) * w)
      n_compute_per_diag <- 1 + w
    }
  } else {
    n_compute_per_diag <- (n-1):1
  }
  
  res <- lapply(1:(n-1),function(i){
    n_on_diagonal <- n - i
    d <- i + 1
    idx <- round(seq(1,n_on_diagonal,length.out=n_compute_per_diag[i]))
    ij <- dk2ij(rep(d,n_compute_per_diag[i]),idx)
    dists <- getSubsampledDistanceMatrix(x,dist.fn,ij)
    return(cbind(ij,dists))
  })
  
  names(res) <- 2:n

  return(res)
}

























estimateSumFromDiagonalSubsamplesSlow <- function(x,
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

computeDiagonallySubsampledDistancesSlow <- function(x,dist.fn,n.per.diag) {
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