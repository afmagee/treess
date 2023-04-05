#' Get all distances for pairs (x[[i]],x[[j]]) of objects specified by the (k x 2) matrix ij
#' 
#' @param x The values
#' @param dist.fn A function which can be called dist.fn(x[[i]],x[[j]])
#' @param ij A 2-column matrix specifying (i,j) pairs for which the distances will be computed
#' @return Distances in order specified by matrix ij
getSubsampledDistanceMatrix <- function(x,dist.fn,ij) {
  sapply(1:dim(ij)[1],function(idx){
    dist.fn(x[[ij[idx,1]]],x[[ij[idx,2]]])
  })
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
