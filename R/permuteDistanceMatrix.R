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
  
  # Can we not just replace this with dmat[new_indices,new_indices]?
  # That would also allow this to be more generally used as permuteMatrix
  tmp <- matrix(0,n,n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp[j,i] <- tmp[i,j] <- dmat[new_indices[i],new_indices[j]]
    }
    
  }

  return(tmp)
}
