#' Constructs intervals for summaries of trees.
#' 
#' For one or more MCMC chains, constructs confidence/prediction intervals for tree and split probabilities.
#'
#' @param x A list of chains containing the MCMC samples of trees as multiPhylo objects.
#' @param ESS Vector of ESS of each of the MCMC chains in x.
#' @param type Either "prediction" for prediction intervals or "confidence" for confidence intervals.
#' @param interval.width Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95\% confidence/prediction intervals.
#' @param method For confidence intervals, method of interval constriction, "Jefreys"|"Wilson"|"ContinuityCorrectedWilson".
#' @details The function returns a list of lists.
#' The top level of lists is $split and $topology, which separates split and topology probabilities.
#' Within each of these is a list, each element containing the intervals for each chain as a matrix.
#' Each matrix has splits (or trees) in rows (rows are comparable across chains), and the columns are the point estimate, lower, and upper CIs.
#' Note that prediction intervals are generated once per chain, and do not account for differences in per-chain ESS.
#' That is, they use the defaule n.new=n in \link{binomialProportionPI}.
#' @return A list of intervals for each chain, see details.
#' @export
#' @seealso \link{binomialProportionCI}, \link{binomialProportionPI}, \link{treeStability}, \link{plotTreeIntervals}
constructTreeIntervals <- function(x, ESS, type, interval.width=0.95, method="Jeffreys") {
  # recover()
  
  # Check for valid inputs
  if ( !(any(c("list","mcmc.list") %in% class(x))) ) {
    stop("Argument 'x' must be a list.")
  }
  
  if ( !("multiPhylo" %in% class(x[[1]])) ) {
    stop("Argument \"x\" must be a list of multiPhylo objects.")
  }
  
  nchains <- length(x)
  
  # ensure all chains are of same dimension(s)
  lens <- unlist(lapply(x,length))
  
  if ( length(unique(lens)) != 1 ) {
    stop("All elements of `x` must be of same length/dimension")
  }
  ngen <- lens[1]
  
  # Make trees coordinates
  all_trees <- do.call(c,x)
  coords <- trees2Coords(all_trees)
  
  nsplits <- dim(coords)[2]
  n_unique_topologies <- NA
  
  # Per-chain split probabilities
  split_probs <- lapply(1:nchains,function(i){
    chain_coords <- coords[(ngen*(i-1)+1):(i*ngen),]
    return(colMeans(chain_coords))
  })
  
  # Per-chain tree probabilities
  tree_strings <- apply(coords,1,paste0,collapse="")
  trees <- as.integer(as.factor(tree_strings))
  n_unique_topologies <- max(trees)
  tree_probs <- lapply(1:nchains,function(i){
    chain_trees <- trees[(ngen*(i-1)+1):(i*ngen)]
    chain_probs <- sapply(1:n_unique_topologies,function(k){sum(chain_trees == k)/ngen})
    return(chain_probs)
  })

  interval_fun <- NULL
  if ( grepl("con",type) ) {
    interval_fun <- function(p,n) {
      binomialProportionCI(p,n,method=method,ci.width=interval.width)
    }
  } else if ( grepl("pre",type) ) {
    interval_fun <- function(p,n) {
      binomialProportionPI(p,n,n.new=n,pi.width=interval.width)
    }
  } else {
    stop("Unrecognized \"type\" option.")
  }
  
  chains_splits <- lapply(1:nchains, function(i){
    tmp <- interval_fun(split_probs[[i]],ESS[i])
    tmp <- cbind(split_probs[[i]],tmp)
    colnames(tmp)[1] <- "point.est"
    return(tmp)
  })
  names(chains_splits) <- paste0("chain_",1:nchains)
  
  chains_trees <- lapply(1:nchains, function(i){
    tmp <- interval_fun(tree_probs[[i]],ESS[i])
    tmp <- cbind(tree_probs[[i]],tmp)
    colnames(tmp)[1] <- "point.est"
    return(tmp)
  })
  names(chains_trees) <- paste0("chain_",1:nchains)
  
  res <- list(chains_splits,chains_trees)
  names(res) <- c("split","topology")
  
  return(res)
}




