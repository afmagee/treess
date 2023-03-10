#' Lists all available methods for computing the ESS of trees/multivariate objects.
#'
#' @param recomended If TRUE, just returns list of the best-performing methods.
#' @export
#' @examples
#' getESSMethods()
getESSMethods <- function(recommended=FALSE) {
  all_methods <- c(
    # ESS = n (not generally useful, but helpful if testing new ESS functions)
    "fixedN",
    # Dimension reduction ESS
    "minPseudoESS",
    "medianPseudoESS",
    "totalDistanceESS",
    "foldedRankMedioidESS",
    "CMDSESS",
    # ESS by Frechet generalizations of 1-D metrics
    "frechetCorrelationESS",
    "splitFrequencyESS",
    # ad-hoc ESS
    "jumpDistanceBootstrapESS",
    "jumpDistanceBootstrapUnsmoothedESS",
    "approximateESS")
  
  if ( recommended ) {
    all_methods <- c("minPseudoESS",
                     "medianPseudoESS",
                     "frechetCorrelationESS")
  }
  
  return(sort(all_methods))
}

# A function that always says the ESS is the number of samples taken
fixedN <- function(dmat,...) {
  return(dim(dmat)[1])
}

#' Calculate one or more effective sample size measures for an MCMC log file containing phylogenetic trees or other multivariate/non-Euclidean objects.
#'
#' @param x A list of chains containing tree-valued (or multivariate or non-Euclidean) variable for which to compute ESS.
#' @param dist.fn A function suitable to calculate distances between all values in x.
#' @param method The method for computing ESS (see details).
#' @param alpha Type I error rate for methods using CIs/hypothesis tests, the proportion of the asymptote used in the approximateESS.
#' @param nsim For methods using permutation/bootstrap resampling, the number of resampling iterations.
#' @param min.nsamples The minimum number of samples do be used in calculating summaries (median distance, correlation, etc.). Essentially the maximum time lag considered is length(chains[[i]]) - min.nsamples. Not applicable to dimension reduction methods.
#' @param max.approximateESS.timelag The maximum time lag considered in the approximateESS, overrides min.nsamples.
#' @export
#' @examples
#' chains <- list(rnorm(100),cumsum(rnorm(100)))
#' treess(chains,dist)
treess <- function(x,dist.fn,methods=getESSMethods(recommended=TRUE),alpha=0.05,nsim=1000,min.nsamples=5,max.approximateESS.timelag=100) {
  
  # Check for valid inputs
  if ( !(any(c("list","mcmc.list") %in% class(x))) ) {
    stop("Argument 'x' must be a list.")
  }
  
  nchains <- length(x)
  
  # ensure all chains are of same dimension(s)
  lens <- unlist(lapply(x,length))
  
  if ( length(unique(lens)) != 1 ) {
    stop("All elements of `x` must be of same length/dimension")
  }
  
  if ( "maxtrix" %in% class(x[[1]]) ) {
    ngen <- dim(x[[1]])[1]
  } else {
    ngen <- lens[1]
  }
  
  # Make sure the chosen method(s) are valid
  if ( any( !(methods %in% getESSMethods()) ) ) {
    stop("Invalid option to argument 'methods'.")
  }
  
  # Compute distance matrix, ESS
  dist_list <- lapply(x,dist.fn)
  
  if ( !all("dist" %in% unlist(lapply(dist_list,class))) ) {
    stop("dist.fn must return object of class dist")
  }
  
  dmat_list <- lapply(dist_list,as.matrix)
  
  ess <- lapply(1:length(dmat_list),function(i){
    dmat <- dmat_list[[i]]
    phy <- x[[i]]
    these_ess <- sapply(methods,function(ess_method){
      eval(call(ess_method,dmat=dmat,trees=phy,min.nsamples=min.nsamples,alpha=alpha,nsim=nsim,max.approximateESS.timelag=max.approximateESS.timelag))
    })
    names(these_ess) <- methods
    return(these_ess)
  })
  
  names(ess) <- paste0("chain_",1:nchains)
  
  return(ess)
}