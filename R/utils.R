#' Finds s0 for a piecewise constant curve by linear interpolation.
#' 
#' Curve is assumed to be nondecreasing.
#'
#' @param x The time lags at which measurements are taken (MUST INCLUDE 0!)
#' @param y The measurement associate with each x that we are using for to compute the ESS.
#' @param y.crit The critical value above which independence is achieved.
#' @keywords internal
finds0Smoothed <- function(x,y,y.crit) {
  n <- length(y)
  
  if ( (!x[1] == 0 && length(x) == n) ) {
    stop("findt0Smoothed requires x[1] = 0 and length(x) = length(y)")
  }
  if ( any(y[-1] < y[-n]) ) {
    stop("findt0Smoothed requires nondecreasing y")
  }
  
  # Find change points and remove replicate entries, track x coordinates
  # x=(1,2,3,4,5,6),y=(1,1,1,2,2,3) -> x=(1,4,6),y=(1,2,3)
  if ( (any(y[-1] == y[-n])) ) {
    first <- match(unique(y),y)
    x <- x[first]
    y <- y[first]
  }
  
  first_larger <- min(which(y > y.crit))
  
  if (x[first_larger] == 0) {
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

#' KL divergence
#' 
#' Calculates Kullbeck-Leibler divergence KL(p||q) for discrete p and q.
#'
#' @param p probability mass function
#' @param q probability mass function
#' @return KL divergence.
#' @export
KL <- function(p,q) {
  summand <- p * log(q/p)
  # 0 * log(0) is 0 for KL purposes
  summand[p < .Machine$double.eps] <- 0
  return(-sum(summand))
}

expandDistanceMatrix <- function(dmat,indices) {
  # recover()
  
  # When computing the variance of the posterior for different ESS, NAs get introduced into the indices
  # This is because ESS <= n in our construction
  indices <- indices[!is.na(indices)]
  
  ntrees <- dim(dmat)[1]
  ngen <- length(indices)
  
  # If all indices are the same index, we know the new distance matrix immediately
  if ( length(unique(indices)) == 1 ) {
    return(matrix(0,ngen,ngen))
  }
  
  # Label the matrix so we know which tree is represented in each row/column
  newdmat <- matrix(NA,ngen,ngen)
  diag(newdmat) <- 0
  for (i in 1:(ngen-1)) {
    for (j in (i+1):ngen) {
      newdmat[i,j] <- newdmat[j,i] <- dmat[indices[i],indices[j]]
    }
  }
  
  colnames(newdmat) <- indices
  row.names(newdmat) <- indices
  
  return(newdmat)
}


# expandDistanceMatrix <- function(dmat,indices) {
#   # recover()
#   
#   # When computing the variance of the posterior for different ESS, NAs get introduced into the indices
#   # This is because ESS <= n in our construction
#   indices <- indices[!is.na(indices)]
#   
#   ntrees <- dim(dmat)[1]
#   ngen <- length(indices)
#   
#   # If all indices are the same index, we know the new distance matrix immediately
#   if ( length(unique(indices)) == 1 ) {
#     return(matrix(0,ngen,ngen))
#   }
#   
#   # Label the matrix so we know which tree is represented in each row/column
#   row.names(dmat) <- 1:ntrees
#   colnames(dmat) <- 1:ntrees
#   
#   # Delete all unneeded rows/columns
#   to_remove <- sort(which(sapply(1:ntrees,function(i){sum(indices == i) == 0})),decreasing=TRUE)
#   for (i in to_remove) {
#     dmat <- dmat[,-i]
#     dmat <- dmat[-i,]
#   }
#   
#   n_used_trees <- dim(dmat)[1]
#   unique_indices <- as.numeric(colnames(dmat))
#   
#   # Duplicate all the rows/columns in order
#   for (i in 1:n_used_trees) {
#     count <- sum(indices == unique_indices[i])
#     # If count == 1, the matrix is unchanged for now. Otherwise we must expand the entries here
#     if (count > 1) {
#       current_dim <- dim(dmat)[1]
#       new_dim <- current_dim + count - 1
#       
#       tmp <- matrix(nrow=new_dim,ncol=new_dim)
#       
#       # Expanding the matrix breaks it up into three chunks (some may have size 0), left, center, and right (or upper, middle, and lower)
#       current_indices <- as.numeric(colnames(dmat))
#       
#       idx <- unique_indices[i]
#       
#       left_current <- which(current_indices < idx)
#       center_current <- which(current_indices == idx)
#       right_current <- which(current_indices > idx)
#       
#       new_names <- c(current_indices[left_current],rep(current_indices[center_current],count),current_indices[right_current])
#       
#       left_new <- left_current
#       center_new <- center_current:(center_current+count-1)
#       right_new <- c()
#       if ( max(center_new) < new_dim ) {
#         right_new <- (max(center_new)+1):(max(center_new)+length(right_current))
#       }
#       
#       # Fill upper-left portion (already expanded) and get first chunk of central portion (the portion to expand)
#       first <- c()
#       if ( length(left_current) > 0 ) {
#         tmp[left_new,left_new] <- dmat[left_current,left_current]
#         first <- dmat[left_new,center_current]
#       }
#       
#       # Fill bottom-right portion (not yet expanded) and get last chunk of central portion (the portion to expand)
#       last <- c()
#       if ( length(right_current) > 0 ) {
#         tmp[right_new,right_new] <- dmat[right_current,right_current]
#         last <- dmat[right_current,center_current]
#       }
#       
#       # Fill upper-right and bottom-left portions
#       if ( (length(left_current) > 0) && (length(right_current) > 0) ) {
#         # upper-right
#         tmp[left_new,right_new] <- dmat[left_current,right_current]
#         # bottom-left
#         tmp[right_new,left_new] <- dmat[right_current,left_current]
#       }
#       
#       # Fill newly-expanded portion
#       new_row <- c(first,rep(0,count),last)
#       tmp[,center_new] <- new_row
#       tmp <- t(tmp)      
#       tmp[,center_new] <- new_row
#       
#       # Done
#       dmat <- tmp
#       colnames(dmat) <- new_names
#       row.names(dmat) <- new_names
#     }
#   }
#   # Re-order the matrix to match the indices
#   blocked_indices <- as.numeric(colnames(dmat))
#   blocked_order <- rank(blocked_indices,ties.method="first")
#   target_order <- rank(indices,ties.method="first")
#   key <- match(target_order,blocked_order)
#   dmat <- dmat[,key]
#   dmat <- dmat[key,]
#   
#   return(dmat)
# }

# Finds all unique trees (ignoring branch lengths) in the input set
uniqueTopologies <- function(trees) {
  
  remaining_trees <- trees
  
  topos <- vector("list")
  counts <- numeric()
  
  iter <- 0
  while( length(remaining_trees) > 0 ){
    iter <- iter + 1
    # Strip branch lengths
    clado <- remaining_trees[[1]]
    clado$edge.length <- NULL
    
    # Add this topology
    topos[[iter]] <- clado
    
    # count matches
    dists <- suppressWarnings(RF.dist(remaining_trees,clado))
    counts <- c(counts,sum(dists == 0))
    
    # Remove matches
    remaining_trees <- remaining_trees[dists > 0]
  }
  
  class(topos) <- "multiPhylo"
  
  return(list(topologies=topos,counts=counts))
}

# Finds all unique trees (with branch lengths) in the input set
uniqueTrees <- function(trees) {
  
  remaining_trees <- trees
  
  utrees <- vector("list")
  counts <- numeric()
  
  iter <- 0
  while( length(remaining_trees) > 0 ){
    iter <- iter + 1
    
    # Add this topology
    utrees[[iter]] <- remaining_trees[[1]]
    
    # count matches
    dists <- suppressWarnings(KF.dist(remaining_trees,remaining_trees[[1]]))
    counts <- c(counts,sum(dists == 0))
    
    # Remove matches
    remaining_trees <- remaining_trees[!(dists == 0)]
  }
  
  class(utrees) <- "multiPhylo"
  
  return(list(trees=utrees,counts=counts))
}
