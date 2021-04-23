# Calculates (\hat{p} - E[\hat{P}])^2 for split probabilities for an object of class simulatedPosterior
splitProbSquaredError <- function(simulated.samples) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  # Not strictly an element of this class of object, but will be added before we ever call splitProbSquaredError
  split_refs <- simulated.samples$coords
  
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

  return(squared_errors)
}


# Calculates (\hat{p} - p)^2 for phylogeny probabilities for an object of class simulatedPosterior
treeProbSquaredError <- function(simulated.samples) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  # split probabilities in each chain
  per_chain <- apply(simulated.samples$indices,2,function(indices){
    tree_probs <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
    tree_probs <- tree_probs/sum(tree_probs)
    return(tree_probs)
  })
  
  # \hat{E}, aka the average average
  E_hat <- rowMeans(per_chain)
  
  squared_errors <- (per_chain - E_hat)^2
  
  return(squared_errors)
}

# Calculates d(MRC(split frequencies from this chain),MRC(pooled split frequencies from all chains)) for an object of class simulatedPosterior.
# The distance measure is RF distance
MRCSquaredError <- function(simulated.samples) {
  # recover()
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  split_refs <- simulated.samples$coords
  
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

# Calculates the MRC tree from an RF coordinate matrix and the probabilities of the trees in the matrix
# Returns a vector of the splits in the MRC tree
MRC.from.coords <- function(coords,probs=rep(1/dim(coords)[1],dim(coords)[1])) {
  if ( sum(probs) != 1) {
    probs <- probs/sum(probs)
  }
  split_probs <- colSums(coords * probs)
  return(as.numeric(split_probs > 0.5 ))
}
