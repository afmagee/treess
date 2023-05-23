#' Calculates the average standard deviation of split frequencies (and other related measures).
#'
#' @param x A list. Either a list of multiPhylo objects, a list of perTreeSplit objects, or a list of split probability vectors.
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
  
  classes <- sapply(x,class)
  if ( length(unique(classes)) != 1 ) {
    stop("All elements in argument \"x\ must be the same length.")
  }
  
  mean_probs <- NULL
  sdsf <- NULL
  if ( classes[1] == "multiPhylo" || classes[1] == "perTreeSplits" ) {
    if ( classes[1] == "multiPhylo" ) {
      x <- lapply(x,perTreeSplits)
    }
    means <- lapply(x,splitFrechetMean)
    tmp <- perTreeSplitsConvergenceHelper(means,return.sdsf=TRUE)
    mean_probs <- tmp$mean
    sdsf <- tmp$sdsf
  } else if ( classes[1] == "numeric" ) {
    splits <- do.call(cbind,x)
    mean_probs <- rowMeans(x)
    sdsf <- apply(splits,1,sd)
  } else {
    stop("Argument \"x\" must either contain trees or split probabilities.")
  }
  
  # Truncation to minimum split frequency
  sdsf <- sdsf[mean_probs > min.freq]

  # Compute summary
  if ( weights == "equal" ) {
    return(summary.fn(sdsf))
  } else if ( weights == "entropy" ) {
    if ( !all.equal(summary.fn,mean) ) {
      stop("Non-equal weighting only allowed for option mean")
    }
    p <- mean_probs[mean_probs > min.freq]
    q <- 1 - p
    e <- -p*log(p)-q*log(q)
    e[is.nan(e)] <- 0
    w <- e/sum(e)
    return(sum(w*sdsf))
  } else {
    stop("Unrecognized option to argument \"weights\".")
  }
  
}