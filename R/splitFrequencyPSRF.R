#' Split-based multi-chain convergence diagnostic.
#'
#' Generalizes 
#' Works like a version of ape::prop.part which returns per-tree (and which considers rooting).
#'
#' @param trees A list of chains containing trees, or list of outputs of \link{perTreeSplits}.
#' @param rooted For rooted trees, if rooted=TRUE, clades are considered rather than splits. Otherwise trees are treated as unrooted.
#' @return The PSRF.
#' @export
#' @seealso \link{splitProbs}, \link{perTreeSplits}
splitFrequencyPSRF <- function(trees,rooted=FALSE) {
  # Check for valid inputs
  if ( !("list" %in% class(trees)) ) {
    stop("Argument 'trees' must be a list.")
  }
  
  nchains <- length(trees)
  
  # ensure all chains are of same dimension(s)
  lens <- unlist(lapply(trees,length))
  
  if ( length(unique(lens)) != 1 ) {
    stop("All elements of 'trees' must be of same length/dimension")
  }
  chainlength <- lens[1]
  
  splits <- NULL
  if ( class(trees[[1]]) == "multiPhylo" || (class(trees[[1]]) == "list" && class(trees[[1]][[1]]) == "phylo") ) {
    all_trees <- do.call(c,trees)
    all_splits <- perTreeSplits(trees,rooted)
    splits <- lapply(1:nchains,function(i){
      all_splits[((i-1)*chainlength+1):(i*chainlength)]
    })
  } else if ( class(trees[[1]]) == "perTreeSplits" ) {
    splits <- trees
  } else {
    stop("Invalid input.")
  }
  
  # recover()
  means <- lapply(splits,splitFrechetMean)
  
  W <- 1/nchains * sum(sapply(1:nchains,function(i){
    splitFrechetVarianceBinomial(splits[[i]],means[[i]])
  }))
  
  grand_mean <- splitFrechetGrandMean(means,TRUE)
  
  B_n <- 1/(nchains - 1) * grand_mean$ssq
  
  V <- (1 - 1/chainlength) * W + B_n
  
  psrf <- sqrt(V/W)
  
  return(psrf)
  
}