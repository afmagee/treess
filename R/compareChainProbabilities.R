#' Constructs intervals for differences (between chains) of summaries of trees.
#' 
#' For two or more MCMC chains, constructs confidence intervals for differences between tree and split probabilities.
#'
#' @param x A list of chains containing the MCMC samples of trees as multiPhylo objects.
#' @param ESS Vector of ESS of each of the MCMC chains in x.
#' @param interval.width Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95\% confidence intervals.
#' @param method Method of confidence interval constriction, "AgrestiCoffo"|"JeffreysPerks" or anything available in DescTools::BinomDiffCI.
#' @details The function returns a list of lists.
#' The top level of lists is $split and $topology, which separates split and topology probabilities.
#' Within each of these is a list, each element containing the intervals for each chain as a matrix.
#' Each matrix has splits (or trees) in rows (rows are comparable across chains), and the columns are the point estimate, lower, and upper CIs, and the effective sample size.
#' @return A list of intervals for each chain, see details.
#' @export
#' @seealso \link{binomialProportionDifferencePI}, \link{plotTreeIntervals}
compareChainProbabilities <- function(x, ESS, interval.width=0.95, method="AgrestiCoffo") {
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
  
  interval_fun <- function(p1,p2,n1,n2) {
    binomialProportionDifferenceCI(p1,p2,n1,n2,method=method,ci.width=interval.width)
  }
  
  chains_splits <- vector("list",choose(length(x),2))
  idx <- 0
  for (i in 1:(nchains-1)) {
    for ( j in (i+1):nchains) {
      idx <- idx + 1
      ess_i <- rep(ESS[i],length(split_probs[[i]]))
      ess_j <- rep(ESS[j],length(split_probs[[j]]))
      tmp <- interval_fun(p1=split_probs[[i]],p2=split_probs[[j]],n1=ess_i,n2=ess_j)
      tmp <- cbind(split_probs[[i]]-split_probs[[j]],tmp,rep(ESS[i],length(split_probs[[i]])),rep(ESS[j],length(split_probs[[j]])))
      colnames(tmp)[1] <- "point.est"
      colnames(tmp)[4] <- "ESS.1"
      colnames(tmp)[5] <- "ESS.2"
      chains_splits[[idx]] <- tmp
      names(chains_splits)[idx] <- paste0("chains_",i,"_vs_",j)
    }
  }
  

  chains_trees <- vector("list",choose(length(x),2))
  idx <- 0
  for (i in 1:(nchains-1)) {
    for ( j in (i+1):nchains) {
      idx <- idx + 1
      ess_i <- rep(ESS[i],length(tree_probs[[i]]))
      ess_j <- rep(ESS[j],length(tree_probs[[j]]))
      tmp <- interval_fun(p1=tree_probs[[i]],p2=tree_probs[[j]],n1=ess_i,n2=ess_j)
      tmp <- cbind(tree_probs[[i]]-tree_probs[[j]],tmp,rep(ESS[i],length(tree_probs[[i]])),rep(ESS[j],length(tree_probs[[j]])))
      colnames(tmp)[1] <- "point.est"
      colnames(tmp)[4] <- "ESS.1"
      colnames(tmp)[5] <- "ESS.2"
      chains_trees[[idx]] <- tmp
      names(chains_trees)[idx] <- paste0("chains_",i,"_vs_",j)
    }
  }
  
  res <- list(chains_splits,chains_trees)
  names(res) <- c("split","topology")
  
  return(res)
}
