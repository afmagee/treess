#' Reads a MrBayes .trprobs file
#'
#' @param filepath Filepath to MrBayes .trprobs file
#' @param CI.width All trees within the posterior CI.width*100 percentile will be returned (1.0 for all trees).
#' @return A list. $trees contains all unique topologies, and $probs their posterior probabilities
#' @export
readTreeProbs <- function(filepath,CI.width=1.0) {
  # recover()
  trees <- ape::read.nexus(filepath)
  txt <- scan(filepath,what=character(),sep="\n")
  prob_txt <- txt[grepl("[&W",txt,fixed=TRUE)]
  prob_txt <- do.call(rbind,strsplit(prob_txt,"[&W",fixed=TRUE))[,2]
  prob_txt <- do.call(rbind,strsplit(prob_txt,"]",fixed=TRUE))[,1]
  probs <- as.numeric(prob_txt)
  if ( CI.width < sum(probs) ) {
    cumprobs <- cumsum(probs)
    last_index <- min(which(cumprobs >= CI.width))
    trees <- trees[1:last_index]
    probs <- probs[1:last_index]
  }
  probs <- probs/sum(probs)
  return(list(trees=trees,probs=probs))
}

#' Trees as vectors of splits.
#'
#' Transforms sample of trees into their coordinates in (reduced) RF space. Reduction is to the set of coordinates (splits/bipartitions) seen in at least one tree.
#'
#' @param trees list of trees or multiPhylo object
#' @return Output: a list. $coords is a matrix of trees as coordinates (see details), and $taxa is a list of taxa in each split in the matrix (ordered 1:n).
#' @details In the output$coords, every row is a tree, every column the presence/absence of a split.
#' @export
as.RFcoords <- function(trees) {
  splits <- trees2Coords(trees,namesplits=TRUE)
  split_taxa <- colnames(splits)
  split_taxa <- lapply(1:dim(splits)[2],function(i){
    strsplit(split_taxa[i],";")[[1]]
  })
  colnames(splits) <- NULL
  return(list(coords=splits,taxa=split_taxa))
}

#' Tree distance
#'
#' Calculates unrooted RF distance between a set of trees.
#'
#' @param trees list of trees or multiPhylo object.
#' @return Distances as an object of class dist.
#' @details Unlike phangorn::RF.dist, does not offer rooted RF distances of warn when unrooting trees.
#' @export
fastRF <- function(trees) {
  coords <- trees2Coords(trees)
  return(dist(coords,method="manhattan"))
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
  indices <- indices[!is.na(indices)]
  
  ntrees <- dim(dmat)[1]
  ngen <- length(indices)
  
  # If all indices are the same index, we know the new distance matrix immediately
  if ( length(unique(indices)) == 1 ) {
    return(matrix(0,ngen,ngen))
  }
  
  # Label the matrix so we know which tree is represented in each row/column
  row.names(dmat) <- 1:ntrees
  colnames(dmat) <- 1:ntrees
  
  # Delete all unneeded rows/columns
  to_remove <- sort(which(sapply(1:ntrees,function(i){sum(indices == i) == 0})),decreasing=TRUE)
  for (i in to_remove) {
    dmat <- dmat[,-i]
    dmat <- dmat[-i,]
  }
  
  n_used_trees <- dim(dmat)[1]
  unique_indices <- as.numeric(colnames(dmat))
  
  # Duplicate all the rows/columns in order
  for (i in 1:n_used_trees) {
    count <- sum(indices == unique_indices[i])
    # If count == 1, the matrix is unchanged for now. Otherwise we must expand the entries here
    if (count > 1) {
      current_dim <- dim(dmat)[1]
      new_dim <- current_dim + count - 1
      
      tmp <- matrix(nrow=new_dim,ncol=new_dim)
      
      # Expanding the matrix breaks it up into three chunks (some may have size 0), left, center, and right (or upper, middle, and lower)
      current_indices <- as.numeric(colnames(dmat))
      
      idx <- unique_indices[i]
      
      left_current <- which(current_indices < idx)
      center_current <- which(current_indices == idx)
      right_current <- which(current_indices > idx)
      
      new_names <- c(current_indices[left_current],rep(current_indices[center_current],count),current_indices[right_current])
      
      left_new <- left_current
      center_new <- center_current:(center_current+count-1)
      right_new <- c()
      if ( max(center_new) < new_dim ) {
        right_new <- (max(center_new)+1):(max(center_new)+length(right_current))
      }
      
      # Fill upper-left portion (already expanded) and get first chunk of central portion (the portion to expand)
      first <- c()
      if ( length(left_current) > 0 ) {
        tmp[left_new,left_new] <- dmat[left_current,left_current]
        first <- dmat[left_new,center_current]
      }
      
      # Fill bottom-right portion (not yet expanded) and get last chunk of central portion (the portion to expand)
      last <- c()
      if ( length(right_current) > 0 ) {
        tmp[right_new,right_new] <- dmat[right_current,right_current]
        last <- dmat[right_current,center_current]
      }
      
      # Fill upper-right and bottom-left portions
      if ( (length(left_current) > 0) && (length(right_current) > 0) ) {
        # upper-right
        tmp[left_new,right_new] <- dmat[left_current,right_current]
        # bottom-left
        tmp[right_new,left_new] <- dmat[right_current,left_current]
      }
      
      # Fill newly-expanded portion
      new_row <- c(first,rep(0,count),last)
      tmp[,center_new] <- new_row
      tmp <- t(tmp)      
      tmp[,center_new] <- new_row
      
      # Done
      dmat <- tmp
      colnames(dmat) <- new_names
      row.names(dmat) <- new_names
    }
  }
  # Re-order the matrix to match the indices
  blocked_indices <- as.numeric(colnames(dmat))
  blocked_order <- rank(blocked_indices,ties.method="first")
  target_order <- rank(indices,ties.method="first")
  key <- match(target_order,blocked_order)
  dmat <- dmat[,key]
  dmat <- dmat[key,]
  
  return(dmat)
}

permuteDistanceMatrix <- function(dmat,bootstrap=FALSE) {
  # recover()
  
  n <- dim(dmat)[1]
  
  new_indices <- sample.int(n,n,replace=bootstrap)
  
  tmp <- dmat
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp[new_indices[i],new_indices[j]] <- dmat[i,j]
    }
    
  }
  newmat <- t(tmp)
  newmat[lower.tri(newmat)] <- tmp[lower.tri(tmp)]
  
  return(newmat)
}


# Transforms sample of trees into their coordinates in (reduced) KF space.
# Reduction is to the set of coordinates/bipartitions seen in at least one tree
# Arguments: 
#   trees: list of trees or multiPhylo object
#   namesplits: should we name the splits in the matrix?
# Output: trees as coordinates
trees2Coords <- function(trees,namesplits=FALSE) {
  # recover()
  
  taxa <- trees[[1]]$tip.label
  
  ntax <- length(taxa)
  
  # Get master list of all splits
  all_splits <- as.matrix(phangorn::as.splits(trees))

  # Order alphabetically
  all_splits <- all_splits[,order(colnames(all_splits))]
  
  # Remove trivial splits
  trivial <- rowSums(all_splits) == 1 | rowSums(all_splits) == ntax  | rowSums(all_splits) == (ntax - 1)
  
  all_splits <- all_splits[!trivial,]
  
  # Polarize, our rule here is that all splits should include the first taxon
  to_polarize <- all_splits[,1] == 0
  
  all_splits[to_polarize,] <- -1 * (all_splits[to_polarize,] - 1)
  
  # Collapse to strings
  split_names <- c()
  if ( namesplits ) {
    taxa <- colnames(all_splits)
    split_names <- apply(all_splits,1,function(x){
      paste0(taxa[as.logical(x)],collapse=";")
    })
  }
  all_splits <- apply(all_splits,1,paste0,collapse="")
  
  coords <- matrix(0,nrow=length(trees),ncol=length(all_splits))
  
  for (i in 1:length(trees)) {
    these_splits <- as.matrix(phangorn::as.splits(trees[[i]]))
    
    # alphabetize
    these_splits <- these_splits[,order(colnames(these_splits))]
    
    # remove trivial splits (only one taxon or all taxa)
    trivial <- rowSums(these_splits) == 1 | rowSums(these_splits) == ntax  | rowSums(all_splits) == (ntax - 1)
    
    these_splits <- these_splits[!trivial,]
    
    # Polarize, our rule here is that splits should be <= 50% 1s, and if 50% the first element should be a 0
    to_polarize <- these_splits[,1] == 0
    
    these_splits[to_polarize,] <- -1 * (these_splits[to_polarize,] - 1)
    
    these_splits <- apply(these_splits,1,paste0,collapse="")
    
    seen <- all_splits %in% these_splits
    
    coords[i,seen] <- 1
    
  }
  
  if ( namesplits ) {
    colnames(coords) <- split_names
  }
  
  return(coords)
  
}


#######
# everything under here needs to be cleaned up if it is to be included, and whether or not to export needs to be decided
#######

spaceTimeRankAutocorrelationTest <- function(x,dist.fn,method=c("tau|kendall","rho|spearman")) {
  # recover()
  n <- ifelse(is.null(dim(x)),length(x),dim(x)[1])
  
  space_dists <- as.matrix(dist(x))
  time_dists <- as.matrix(dist(1:n))
  
  if ( method == "tau" ) {
    method <- "kendall"
  } else if ( method == "rho" ) {
    method <- "spearman"
  } else {
    stop("Invalid 'method' specified")
  }
  
  cor.test(space_dists,time_dists,method=method)
  
}

chatterjeeCor <- function(x,y) {
  if ( !((is.numeric(x) || is.logical(x)) && (is.numeric(y) || is.logical(y))) ) {
    stop("x and y must be numeric")
  }
  
  if ( length(x) != length(y) ) {
    stop("x and y must be of same length")
  }
  
  n <- length(x)
  
  key <- order(x)
  x <- x[key]
  y <- y[key]
  r <- rank(y,ties.method="max")
  l <- rank(-y,ties.method="max")
  
  num <- n * sum(abs(r[2:n] - r[1:(n-1)]))
  denom <- 2 * sum(l * (n - l))
  
  return(1 - num/denom)
}

chatterjeeCor.isSig <- function(x,y,conf.level=0.95,nboot=200) {
  # recover()
  
  if ( conf.level < 0 || conf.level > 1 ) {
    stop("'conf.level' must be a single number between 0 and 1")
  }
  alpha <- 1 - conf.level
  
  obs <- chatterjeeCor(x,y)
  null <- sapply(1:nboot,function(i){
    x_rep <- sample(x,replace=TRUE)
    y_rep <- sample(y,replace=TRUE)
    chatterjeeCor(x_rep,y_rep)
  })
  upper <- quantile(null,1-alpha/2,na.rm=TRUE)
  lower <- quantile(null,alpha/2,na.rm=TRUE)
  sig <- obs < lower | obs > upper
  names(sig) <- NULL
  return( sig )
}

l1dist <- function(x) {
  dist(x,method="manhattan")
}

s.sq.dist.cov <- function(x,y,dist.fn) {
  n <- length(x)
  
  xdists <- dist.fn(x)^2
  if ( class(xdists) != "dist" ) {
    stop("dist.fn must return object of class dist")
  }
  ydists <- dist.fn(y)^2
  
  xvar <- sum(xdists)/choose(n,2)/2
  yvar <- sum(ydists)/choose(n,2)/2
  
  dxy <- mean(sapply(1:n,function(i){dist.fn(c(x[[i]],y[[i]]))^2}))
  
  return((xvar + yvar - dxy)/2)
}

s.sq.dist.cor <- function(x,y,dist.fn) {
  n <- length(x)
  
  xdists <- dist.fn(x)^2
  if ( class(xdists) != "dist" ) {
    stop("dist.fn must return object of class dist")
  }
  ydists <- dist.fn(y)^2
  
  xvar <- sum(xdists)/choose(n,2)/2
  yvar <- sum(ydists)/choose(n,2)/2
  
  dxy <- mean(sapply(1:n,function(i){dist.fn(c(x[[i]],y[[i]]))^2}))
  
  covar <- (xvar + yvar - dxy)/2
  
  return(covar/sqrt(xvar*yvar))
}

distanceMMPADCorrelation <- function(x,y,dist.fn,use.sum=FALSE) {
  xdists <- dist.fn(x)
  if ( class(xdists) != "dist" ) {
    stop("dist.fn must return object of class dist")
  }
  xdmat <- as.matrix(xdists)
  
  ydists <- dist.fn(y)
  ydmat <- as.matrix(ydists)
  
  # recover()
  
  if ( use.sum ) {
    xr <- rank(rowSums(xdmat))
    yr <- rank(rowSums(ydmat))
  } else {
    xr <- rank(apply(xdmat,1,median))
    yr <- rank(apply(ydmat,1,median))
  }
  
  X <- qnorm((xr+0.5)/(max(xr)+1))
  Y <- qnorm((yr+0.5)/(max(yr)+1))
  
  return(cor(X,Y))
  
}

distanceBinomialCorrelationTest <- function(x,y,dist.fn) {
  # recover()
  n <- length(x)
  neighbors <- sapply(1:n,function(i){dist.fn(c(x[i],y[i]))})
  all <- sapply(1:n,function(i){
    x_to_y <- dist.fn(c(x[i],y))[1:n]
    y_to_x <- dist.fn(c(y[i],x))[1:n]
    return(median(c(x_to_y,y_to_x)))
  })
  binom.test(sum(neighbors < all),n=n,p=0.5)$p.value
}

distanceBinomialCorrelationTestESS <- function(x,dist.fn) {
  # recover()
  dmat <- as.matrix(dist.fn(x))
  n <- dim(dmat)[1]
  thin <- n
  
  for (i in 1:(n-5)) {
    neighbors <- dmat[row(dmat) == col(dmat) + i]
    all <- sapply(1:(n-i),function(j){
      median(c(dmat[i,-c(1:i)],dmat[(n-i+1),1:(n-i+1)]))
    })
    p <- binom.test(sum(neighbors < all),n=n,p=0.5)$p.value
    if (p >= 0.05) {
      thin <- i
      break
    }
  }
  return(n/thin)
}

# Takes a sample of pairwise distances between sampled trees, estimates the mean and variance of the average pairwise distance
# These distances must be PURELY TOPOLOGICAL (RF/SPR not KF/geodesic)
# Option "squared" controls whether distances (FALSE) or squared distances (TRUE) are used
estimateTreeTopologicalDistanceMoments <- function(dmat,squared=TRUE) {

  if ( squared ) {
    dmat <- dmat^2
  }
  
  counts <- numeric()
  
  # de-duplicate all trees
  iter <- 0
  while ( sum(dmat == 0) > dim(dmat)[1] ) {
    iter <- iter + 1
    
    # Number of trees that are of this topology
    matches <- which(dmat[iter,] == 0)
    
    if ( length(matches) > 1 ) {
      dmat <- dmat[-matches[-1],]
      dmat <- dmat[,-matches[-1]]
    }
    counts <- c(counts,length(matches))
  }
  
  # count last trees
  counts <- c(counts,rep(1,dim(dmat)[1] - length(counts)))
  
  # tree probabilities
  tprobs <- counts/sum(counts)
  
  # probabilities of each kind of tree to tree comparison
  cprobs <- outer(tprobs,tprobs)
  
  # mean and variance
  m <- sum(cprobs * dmat)
  v <- sum(cprobs * (dmat - m)^2)
  
  res <- c(m,v)
  names(res) <- c("mean","variance")
  
  return(res)
  
}

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

# Takes in a set of trees and finds the medioid tree(s).
# Not guaranteed to find a unique tree for all distance functions.
# Arguments
#   trees: a multiPhylo object
#   method: wRF or KF for computing distances between trees
# Returns: the median tree(s)
treeMedioid <- function(trees,method) {
  if (class(trees) != "multiPhylo") {stop("trees must be of class \"multiPhylo\" ")}
  
  # Distance matrix of all trees
  if ( method %in% c("kf","KF") ) {
    all_dists <- as.matrix(KF.dist(trees))
  } else if ( method %in% c("wrf","wRF","WRF") ) {
    all_dists <- as.matrix(wRF.dist(trees))
  } else {
    stop("Unrecognized method for tree distance")
  }
  
  # Per tree sum of distances
  per_tree_total <- rowSums(all_dists)
  
  return(trees[[which(per_tree_total == min(per_tree_total))]])
  
}
