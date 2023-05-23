#' Implicitly generate permutations of distance matrices.
#'
#' For a distance matrix D = d_ij, generates new distance matrix resulting from permuting or bootstrapping the samples.
#' Avoids the need to recompute distances when generating many permuted/bootstrapped distance matrices.
#' 
#' @param dmat matrix of distances
#' @param bootstrap TRUE/FALSE, if TRUE sampling is done with replacement, if FALSE sampling is done without replacement (permutation only)
#' @return The permuted or bootstrapped distance matrix.
#' @export
permuteDistanceMatrix <- function(dmat,bootstrap=FALSE) {
  # recover()
  
  n <- dim(dmat)[1]
  
  new_indices <- sample.int(n,n,replace=bootstrap)
  
  return(dmat[new_indices,new_indices])
}
