# Calculates (\hat{p} - E[\hat{P}])^2 for split probabilities for an object of class simulatedPosterior
# Takes argument dmat only to allow use of call()
splitProbSquaredError <- function(simulated.samples) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  split_refs <- constructSplitReferences(simulated.samples$trees)
  
  # split probabilities in each chain
  per_chain <- apply(simulated.samples$indices,2,function(indices){
    tree_probs <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
    tree_probs <- tree_probs/sum(tree_probs)
    est <- colSums(split_refs * tree_probs)
    return(est)
  })
  
  # \hat{E}, aka the average average
  E_hat <- rowMeans(per_chain)
  
  squared_errors <- (per_chain - E_hat)^2

  # Splits that exist for our purposes
  # Any split that has not been seen at least once does not count
  is_in_samples <- rowSums(per_chain) > 0
  
  return(squared_errors[is_in_samples,])
}


# Calculates (\hat{p} - p)^2 for phylogeny probabilities for an object of class simulatedPosterior
treeProbSquaredError <- function(simulated.samples) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  split_refs <- constructSplitReferences(simulated.samples$trees)
  
  # split probabilities in each chain
  per_chain <- apply(simulated.samples$indices,2,function(indices){
    tree_probs <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
    tree_probs <- tree_probs/sum(tree_probs)
    return(tree_probs)
  })
  
  # \hat{E}, aka the average average
  E_hat <- rowMeans(per_chain)
  
  squared_errors <- (per_chain - E_hat)^2
  
  # Splits that exist for our purposes
  # Any split that has not been seen at least once does not count
  is_in_samples <- rowSums(per_chain) > 0
  
  return(squared_errors[is_in_samples,])
}

# Calculates d(MRC(split frequencies from this chain),MRC(pooled split frequencies from all chains)) for an object of class simulatedPosterior.
# The distance measure is RF distance
MRCSquaredError <- function(simulated.samples) {
  # recover()
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  split_refs <- constructSplitReferences(simulated.samples$trees)
  
  # per-chain MRC trees
  per_chain <- t(sapply(1:nchains,function(j){
    indices <- simulated.samples$indices[,j]
    indices <- indices[!is.na(indices)]
    probs <- sapply(1:ntrees,function(i){
      sum(indices == i)
    })/length(indices)
    MRC.from.coords(split_refs,probs)
  }))
  
  # pooled split probs
  tree_probs <- sapply(1:ntrees,function(i){sum(simulated.samples$indices == i,na.rm=TRUE)})
  tree_probs <- tree_probs/sum(tree_probs)
  pooled_split_probs <- colSums(split_refs * tree_probs)

  best_estimate <- as.numeric(pooled_split_probs > 0.5)
  
  # squared distances to best summary
  all_dists <- sapply(1:nchains,function(j){
    dist(rbind(best_estimate,per_chain[j,]),method="manhattan")^2
  })
  return(all_dists)
}



# Constructs a reference for the set of trees
# tree[[i]] contains every split where constructSplitReferences(trees)[i,] == 1
constructSplitReferences <- function(trees) {
  # recover()
  
  ntrees <- length(trees)
  taxa <- trees[[1]]$tip.label
  ntax <- length(taxa)
  
  # Get master list of all splits
  all_splits <- as.matrix(phangorn::as.splits(trees))
  
  # Order taxa alphabetically
  all_splits <- all_splits[,order(colnames(all_splits))]
  
  # Remove trivial splits
  trivial <- rowSums(all_splits) == 1 | rowSums(all_splits) == ntax
  
  all_splits <- all_splits[!trivial,]
  
  # Polarize, our rule here is that sort(taxa)[1] must be in the split
  to_polarize <- all_splits[,1] != 1
  
  all_splits[to_polarize,] <- -1 * (all_splits[to_polarize,] - 1)
  
  # Collapse to strings
  all_splits <- apply(all_splits,1,paste0,collapse="")
  
  # Find which splits are in which trees
  contains_splits <- matrix(0,nrow=length(trees),ncol=length(all_splits))
  for (i in 1:length(trees)) {
    # splits objects are annoying
    these_splits <- as.matrix(phangorn::as.splits(trees[[i]]))
    
    # alphabetize
    these_splits <- these_splits[,order(colnames(these_splits))]
    
    # remove trivial splits (only one taxon or all taxa)
    trivial <- rowSums(these_splits) == 1 | rowSums(these_splits) == ntax
    these_splits <- these_splits[!trivial,]

    # Polarize, our rule here is that splits should be <= 50% 1s, and if 50% the first element should be a 0
    to_polarize <- these_splits[,1] != 1
    these_splits[to_polarize,] <- -1 * (these_splits[to_polarize,] - 1)
    
    # Find which splits are in this tree
    these_splits <- apply(these_splits,1,paste0,collapse="")
    seen <- all_splits %in% these_splits
    for (j in which(seen)) {
      contains_splits[i,j] <- 1
    }
  }
  
  return(contains_splits)
  
}

# Calculates the MRC tree from an RF coordinate matrix and the probabilities of the trees in the matrix
# Returns a vector of the splits in the MRC tree
MRC.from.coords <- function(coords,probs=rep(1/dim(coords)[1],dim(coords)[1])) {
  if ( sum(probs) != 1) {
    probs <- probs/sum(probs)
  }
  split_probs <- colSums(coords * probs)
  return(as.numeric(split_probs > 0.5 ))
}
