# For block permutation tests
block.permute <- function(x,b) {
  if (b == 1) {
    return(sample(x))
  } else {
    n <- length(x)
    n.blocks <- floor(n/b)
    first.cut <- sample.int(b,1)
    cuts <- c(0,first.cut,first.cut + b * (1:(n.blocks-1)))
    if (cuts[length(cuts)] != n) {
      cuts <- c(cuts,n)
    }
    blocks <- lapply(1:(length(cuts)-1),function(i){
      (cuts[i]+1):(cuts[i+1])
    })
    blocks <- blocks[sample.int(length(blocks))]
    return(x[unlist(blocks)])
  }
}

# Checks compatibility of rooting of trees with whether user has specified they should be treated as rooted or not
# Unroots trees if needed, throws error if needed, throws warning if needed
checkRootedOption <- function(trees,rooted) {
  if ( any(ape::is.rooted.multiPhylo(trees)) ) {
    if ( !rooted ) {
      warning("Counting splits on rooted trees, instead of clades.")
      trees <- ape::unroot.multiPhylo(trees)
    } 
  } else if ( rooted ) {
    stop("Trees are unrooted but rooted=TRUE")
  }
  return(trees)
}


#' Gets a sparse representation of a spanning tree from a distance matrix.
#' 
#' @param x A distance matrix
#' @param R Either a number of MSTs to collect (to account for non-uniqueness) or "automatic" to chose the number based on the stabilization of the Holmes test statistic
#' @param labels The labels, as in holmesTest
#' @param se.cutoff For eventual use in auto-determination of number of MSTs needed.
#' @return Vector of all edges in the spanning tree (contains both i->j and j->i for any connected values i and j)
#' @keywords internal
getSparseSpanningTreeListHolmes <- function(x,R,labels,se.cutoff=0.25) {
  sparse_spanning_trees <- list()
  if ( is.numeric(R) ) {
    sparse_spanning_trees <- lapply(1:R,function(i){
      getSparseSpanningTree(x,shuffle=TRUE)
    })
  } else if ( grepl("auto",R) ) {
    stop("Auto-terminating is not ready for use.")
    sparse_spanning_trees <- vector("list",100)
    stats <- numeric(100)
    
    # Start with a large enough sample to get a sense of the true average
    nstart <- 20
    sparse_spanning_trees[1:nstart] <- lapply(1:nstart,function(i){getSparseSpanningTree(x,shuffle=TRUE)})
    stats[1:nstart] <- sapply(1:nstart,function(i){holmesTestStat(sparse_spanning_trees[[i]],labels)})
    
    idx <- nstart
    se <- sd(stats[1:idx])/sqrt(idx)
    while (se > se.cutoff) {
      idx <- idx + 1
      if ( idx > length(sparse_spanning_trees) ) {
        sparse_spanning_trees <- c(sparse_spanning_trees,vector("list",100))
        stats <- c(stats,numeric(100))
      }
      sparse_spanning_trees[[idx]] <- getSparseSpanningTree(x,shuffle=TRUE)
      stats[idx] <- holmesTestStat(sparse_spanning_trees[[idx]],labels)
      se <- sd(stats[1:idx])/sqrt(idx)
      # cat(idx,": ",se,"\n")
    }
    sparse_spanning_trees <- sparse_spanning_trees[1:idx]
  } else {
    stop("Invalid input for argument \"R\"")
  }
  return(sparse_spanning_trees)
}


#' Gets a sparse representation of a spanning tree from a distance matrix.
#' 
#' @param x A distance matrix
#' @param shuffle If true, the distance matrix is permuted before the MST is obtained, then re-sorted to match the original order. Potentially useful if MST is not unique.
#' @return Vector of all edges in the spanning tree (contains both i->j and j->i for any connected values i and j)
#' @keywords internal
getSparseSpanningTree <- function(x,shuffle=FALSE) {
  spanning_tree <- NULL
  
  # recover()
  
  # Permuting the distance matrix allows us to get different MSTs when the MST is not unique
  permute <- NULL
  key <- NULL
  if (shuffle) {
    n <- dim(x)[1]
    permute <- sample.int(n)
    x <- x[permute,permute]
  }
  
  if ( requireNamespace("igraph",quietly=TRUE) ) {
    igraph_x <- igraph::graph.adjacency(x,weighted=TRUE,mode="undirected")
    igraph_spanning_tree <- igraph::mst(igraph_x)
    spanning_tree <- igraph::as_adjacency_matrix(igraph_spanning_tree,sparse=FALSE)
  } else {
    spanning_tree <- ape::mst(x)
  }
  
  # recover()
  
  # This will make things much faster for computing test statistics
  sparse_spanning_tree <- lapply(1:dim(spanning_tree)[1],function(i){
    j <- which(spanning_tree[i,] == 1)
    return(cbind(i,j))
  })
  sparse_spanning_tree <- do.call(rbind,sparse_spanning_tree)
  
  # Sort the spanning tree to match the original ordering, if needed
  if ( shuffle ) {
    # sparse_spanning_tree[i,j] gives us the index in permute, k for the tree at hand
    # we want the index, l, in the original ordering of trees
    # because the original ordering is just 1:n, permute[k] = l
    sparse_spanning_tree[,1] <- permute[sparse_spanning_tree[,1]]
    sparse_spanning_tree[,2] <- permute[sparse_spanning_tree[,2]]
  }
  row.names(sparse_spanning_tree) <- NULL
  
  return(sparse_spanning_tree)
}

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
  
  # TODO can we not apply the same dmat[indices,indices] approach here as in permuteDistanceMatrix?
  
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
