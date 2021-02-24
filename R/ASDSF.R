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