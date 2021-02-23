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

# Calculates (\hat{p} - p)^2 for split probabilities for an object of class simulatedPosterior
# Takes argument dmat only to allow use of call()
splitProbSquaredError <- function(simulated.samples,dmat=NULL) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  split_refs <- constructSplitReferences(simulated.samples$trees)
  
  real_probs <- colSums(split_refs*simulated.samples$probs)
  
  squared_errors <- apply(simulated.samples$indices,2,function(indices){
    tree_probs <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
    tree_probs <- tree_probs/sum(tree_probs)
    est <- colSums(split_refs * tree_probs)
    return((est - real_probs)^2)
  })
  return(squared_errors)
}


# Calculates (\hat{p} - p)^2 for phylogeny probabilities for an object of class simulatedPosterior
# Takes argument dmat only to allow use of call()
treeProbSquaredError <- function(simulated.samples,dmat=NULL) {
  # recover()
  ntrees <- length(simulated.samples$trees)
  
  squared_errors <- apply(simulated.samples$indices,2,function(indices){
    est <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
    est <- est/sum(est)
    return((est - simulated.samples$probs)^2)
  })
  return(squared_errors)
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


# Calculates d(MRC(posterior samples),MRC(true posterior)) for an object of class simulatedPosterior, where MRC(trees) is the MRC for those trees
# The distance measure is either RF distance or SPR distance
distanceToTrueMRC <- function(simulated.samples,coords,tree.dist=c("RF","SPR")) {
  # recover()
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
    all_mrc <- t(sapply(1:nchains,function(j){
      indices <- simulated.samples$indices[,j]
      indices <- indices[!is.na(indices)]
      probs <- sapply(1:ntrees,function(i){
        sum(indices == i)
      })/length(indices)
      MRC.from.coords(coords,probs)
    }))
    true_mrc <- MRC.from.coords(coords,simulated.samples$probs)
    
    if ( tree.dist == "RF" ) {
      all_dists <- dist(rbind(true_mrc,all_mrc),method="manhattan")  
      return(all_dists[1:nchains])
    } else {
      stop("Cannot compute distanceToTrueMRC for option SPR")
      true_mrc <- coords2tree(true_mrc)
      all_dists <- sapply(1:nchains,function(j){
        this_mrc <- coords2tree(all_mrc[j,])
        return(phangorn::SPR.dist(this_mrc,true_mrc))
      })
      return(all_dists)
    }
}

# Computes the expected average squared pairwise distance given a distance matrix and tree probabilities
expectedAveragePairwiseSquaredDistance <- function(dmat,probs) {
  sum(outer(probs,probs) * dmat^2)
}

treeVarianceSquaredError <- function(simulated.samples,dmat) {
  # recover()
  nchains <- dim(simulated.samples$indices)[2]
  
  all_ssd <- sapply(1:nchains,function(j){
    indices <- simulated.samples$indices[,j]
    indices <- indices[!is.na(indices)]
    if ( length(indices) >= 2 ) {
      dm <- expandDistanceMatrix(dmat,indices)
      npairs <- choose(length(indices),2)
      return(sum(dm[upper.tri(dm)]^2)/npairs/2)
    } else {
      return(NA)
    }
  })
  
  expected <- expectedAveragePairwiseSquaredDistance(dmat,simulated.samples$probs)
  
  return((all_ssd - expected)^2)
}

# Calculates the index of the medioid (summary = median) or mean-based analog (summary = mean) from a weighted distance matrix
summaryTreeIndex <- function(dmat,weights,summary=c("mean","median")) {
  if ( summary == "mean" ) {
    dmat <- dmat^2
  } else {
    if ( !(summary == "medioid") ) {
      stop("Can only calculate true summary tree for mean or median")
    }
  }
  dists <- colSums(dmat * weights)
  return(which(dists == min(dists)))
}

# # Calculated d(summary(posterior samples),summary(true posterior)) for an object of class simulatedPosterior
# # The summary is either the medioid (summary=median) tree or the mean-based analog the Karcher mean (summary == mean)
# # The distance measure is either RF distance or SPR distance
# distanceToTrueSummary <- function(simulated.samples,dmat,tree.dist=c("RF","SPR"),summary=c("mean","median")) {
#   # recover()
#   
#   ntrees <- length(simulated.samples$trees)
#   nchains <- dim(simulated.samples$indices)[2]
#   
#   true_summary_index <- summaryTreeIndex(dmat,weights=simulated.samples$probs,summary=summary)
#   
#   if ( length(true_summary_index) > 1 ) {
#     warning("There is no unique true summary tree for the true posterior, taking first.")
#     true_summary_index <- true_summary_index[1]
#   }
#   true_summary <- simulated.samples$trees[[true_summary_index]]
#   
#   estimated_summary_indices <- sapply(1:nchains,function(j){
#     est_probs <- sapply(1:ntrees,function(i){sum(simulated.samples$indices[,j] == i,na.rm=TRUE)})
#     est_probs <- est_probs/sum(est_probs)
#     idx <- summaryTreeIndex(dmat,est_probs,summary=summary)
#     if ( length(idx) > 1 ) {
#       warning(paste0("There is no unique true summary tree for chain ",j,", taking first."))
#       idx <- idx[1]
#     }
#     return(idx)
#   })
#   
#   d <- numeric(nchains)
#   if ( tolower(tree.dist) == "rf" ) {
#     d <- phangorn::RF.dist(simulated.samples$trees[estimated_summary_indices],true_summary)
#   } else if ( tolower(tree.dist) == "spr" ) {
#     d <- phangorn::SPR.dist(simulated.samples$trees[estimated_summary_indices],true_summary)
#   } else {
#     stop("Invalid tree distance metric.")
#   }
#   
#   return(d)
# }
# 
# # Wrapper for distanceToTrueSummary(tree.dist=RF,summary="mean")
# RFDistanceToTrueMean <- function(simulated.samples,dmat) {
#   distanceToTrueSummary(simulated.samples,dmat,tree.dist="RF",summary="mean")
# }


# # Computes KL between two sets of trees using probabilities of trees
# treeProbKL <- function(simulated.samples,dmat) {
#   # recover()
#   ntrees <- length(simulated.samples$trees)
#   
#   kl <- apply(simulated.samples$indices,2,function(indices){
#     est <- sapply(1:ntrees,function(i){sum(indices == i,na.rm=TRUE)})
#     return(KL(est,simulated.samples$probs))
#   })
#   
#   return(kl)
#   
# }
# 
# # Convenience function for KL for discrete distributions
# KL <- function(p,q) {
#   if ( any(p == 0) ) {
#     q <- q[p > 0]
#     p <- p[p > 0]
#   }
#   return(-sum(p*log(q/p)))
# }
