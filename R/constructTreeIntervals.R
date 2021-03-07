#' #' Calculate a single effective sample size measure for a tree, or other multivariate/non-Euclidean object.
#' #'
#' #' @param trees The tree-valued (or multivariate or non-Euclidean) variable for which to compute ESS.
#' #' @param dist.fn A function suitable to calculate distances between all values in x.
#' #' @param method The method for computing ESS (see details).
#' #' @param alpha Type I error rate for methods using CIs/hypothesis tests, the proportion of the asymptote used in the approximateESS.
#' #' @param nsim For methods using permutation/bootstrap resampling, the number of resampling iterations.
#' #' @param min.nsamples The minimum number of samples do be used in calculating summaries (median distance, correlation, etc.). Essentially the maximum time lag considered is length(chains[[i]]) - min.nsamples. Not applicable to dimension reduction methods.
#' #' @export
#' #' @examples
#' #' chains <- list(rnorm(100),cumsum(rnorm(100)))
#' #' treess(chains,dist)
#' x <- function(trees,) {
#'   # Make trees coordinates
#'   coords <- matrix()
#'   if ( "multiPhylo" %in% class(trees) ) {
#'     coords <- trees2Coords(trees)
#'   } else if ( "matrix" %in% class(trees) ) {
#'     coords <- trees
#'   } else {
#'     stop("Argument \"trees\" must either be a multiPhylo or an RF coordinate matrix.")
#'   }
#'   
#'   ntrees <- dim(coords)[1]
#'   nsplits <- dim(coords)[2]
#'   n_unique_topologies <- NA
#'   
#' }
#' 
#' 
#' 
#' 
