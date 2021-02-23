#' Calculates the average standard deviation of split frequencies (and other related measures).
#'
#' @param x A list. Either a list of multiPhylo objects or a list of split probability vectors.
#' @param min.freq Splits less than this frequency (on average) will not be included.
#' @param summary.fn Function. The summary of the SDSF. For ASDSF, use mean, for max SDSF, use max.
#' @param weights Allows the ASDSF to be weighted. Options are "equal" for standard calculations or "entropy" to weight by split entropy. Only valid for summary.fn=mean
#' @return Output: the ASDSF (or other summary of the SDSF).
#' @export
# Takes a distance matrix on unique trees and makes a distance matrix on all trees in the order sampled in indices
ASDSF <- function(x,min.freq=0.01,summary.fn=mean,weights="equal") {
  # recover()
  
  if ( !("list" %in% class(x)) && length(unique(lengths(x)) == 1) ) {
    stop("Argument \"x\ must be a list and all elements must be the same length.")
  }
  
  if ( "multiPhylo" %in% class(x[[1]]) ) {
    ntrees <- length(x[[1]])
    nchains <- length(x)
    all_coords <- trees2Coords(do.call(c,x))
    x <- lapply(1:nchains,function(i){
      these_coords <- all_coords[((i-1)*ntrees+1):(i*ntrees),]
      return(colMeans(these_coords))
    })
  } else if ( "numeric" %in% class(x[[1]]) ) {
  } else {
    stop("Argument \"x\" must either contain trees or split probabilities.")
  }
  
  # Truncation to minimum split frequency
  splits <- do.call(cbind,x)
  splits <- splits[rowMeans(splits) > min.freq,]
  
  # SDSF
  sdsf <- apply(splits,1,sd)
  
  # Compute summary
  if ( weights == "equal" ) {
    return(summary.fn(sdsf))
  } else if ( weights == "entropy" ) {
    if ( !all.equal(summary.fn,mean) ) {
      stop("Non-equal weighting only allowed for option mean")
    }
    p <- rowMeans(splits)
    q <- 1 - p
    e <- -p*log(p)-q*log(q)
    e[is.nan(e)] <- 0
    w <- e/sum(e)
    return(sum(w*sdsf))
  } else {
    stop("Unrecognized option to argument \"weights\".")
  }
  
}
#' Finds t0 for a piecewise constant curve by linear interpolation.
#'
#' @param x The time lags at which measurements are taken (MUST INCLUDE 0!)
#' @param y The measurement associate with each x that we are using for to compute the ESS.
#' @param y.crit The critical value above (or below) which independence is achieved.
#' @param above.or.below Is t0 the time at which y >= y.crit ("above") or y <= y.crit ("below")?
#' @keywords internal
findt0Smoothed <- function(x,y,y.crit,above.or.below) {
  if ( (!x[1] == 0 && length(x) == length(y)) ) {
    stop("findt0Smoothed requires x[1] = 0 and length(x) = length(y)")
  }
  if ( above.or.below == "below" ) {
    y <- -y
    y.crit <- -y.crit
  } else if ( above.or.below != "above" ) {
    stop("Unrecognized option to findt0Smoothed.")
  }
  
  # de-duplicate entries
  n <- length(y)
  if ( (any(y[-1] == y[-n])) ) {
    first <- match(unique(y),y)
    x <- x[first]
    y <- y[first]
  }
  
  first_larger <- min(which(y >= y.crit))
  
  if (first_larger == 1) {
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
#' @param coords matrix of trees as coordinates, as.RFcoords(trees)$coords.
#' @return Distances as an object of class dist.
#' @export
RF.from.coords <- function(coords) {
  if ( (length(dim(coords)) !=2) || any(!(as.numeric(coords) %in% c(0,1))) ) {
    stop("Input must be coordinate matrix")
  }
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
  # all_splits <- apply(all_splits,1,paste0,collapse="")
  # This is faster. I don't know why, but it is
  all_splits <- apply(all_splits,1,function(x){paste0(as.raw(x),collapse="")})
  
  coords <- matrix(0,nrow=length(trees),ncol=length(all_splits))
  
  for (i in 1:length(trees)) {
    these_splits <- as.matrix(phangorn::as.splits(trees[[i]]))
    
    # alphabetize
    these_splits <- these_splits[,order(colnames(these_splits))]
    
    # remove trivial splits (only one taxon or all taxa)
    trivial <- rowSums(these_splits) == 1 | rowSums(these_splits) == ntax  | rowSums(these_splits) == (ntax - 1)
    
    these_splits <- these_splits[!trivial,]
    
    # Polarize, our rule here is that splits should be <= 50% 1s, and if 50% the first element should be a 0
    to_polarize <- these_splits[,1] == 0
    
    these_splits[to_polarize,] <- -1 * (these_splits[to_polarize,] - 1)
    
    # these_splits <- apply(these_splits,1,paste0,collapse="")
    # This is faster. I don't know why, but it is
    these_splits <- apply(these_splits,1,function(x){paste0(as.raw(x),collapse="")})
    
    # seen <- all_splits %in% these_splits
    
    coords[i,all_splits %in% these_splits] <- 1
    
  }
  
  if ( namesplits ) {
    colnames(coords) <- split_names
  }
  
  return(coords)
  
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
