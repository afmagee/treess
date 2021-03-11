#' #' Calculate a single effective sample size measure for a tree, or other multivariate/non-Euclidean object.
#' #'
#' #' @param trees The sample of trees from the posterior. A multiPhylo object, or a coordinate matrix as from as.RFcoords(trees)$coords.
#' #' @param summary Should intervals be computed for split probabilities ("split") tree probabilties ("topology").
#' #' @param ESS The effective sample size of the MCMC chain as a single numeric. Alternately, "per.split" to allow per-split ESS (invalid for summary="topology").
#' #' @param type Either "prediction" for prediction intervals or "confidence" for confidence intervals.
#' #' @param interval.width Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95% confidence/prediction intervals.
#' #' @param method For confidence intervals, method of interval constriction, "Jefreys"|"Wilson"|"ContinuityCorrectedWilson:.
#' #' @export
#' #' @seealso \link{binomialProportionCI}, \link{binomialProportionPI}, \link{treeStability}
#' constructTreeIntervals <- function(trees, stat, ESS, type, interval.width=0.95, method="Jeffreys") {
#' 
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
