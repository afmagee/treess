#' Obtains estimates of the distribution of error and variance measures from simulated phylogenetic MCMC.
#'
#' @param simulated.samples An object of class simulatedPosterior (output of simulatePhylogeneticMCMC).
#' @param tree.dist The distance measure for trees (RF or SPR, only one).
#' @param measures The error or variance measure(s) (see details).
#' @return A named list with the distribution for each measure.
#' @export
#' @examples
#' countTrees(rmtree(100,10))
bruteForceMCMCSE <- function(simulated.samples,tree.dist=c("RF","SPR"),measures=c("treeProbSquaredError","treeVarianceSquaredError","splitProbSquaredError","distanceToTrueMRC")) {
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  # Start with the distance matrix for all unique topologies (allows us to reduce compute time)
  dmat <- matrix(nrow=ntrees,ncol=ntrees)
  if ( tolower(tree.dist) == "rf" ) {
    dmat <- as.matrix(phangorn::RF.dist(simulated.samples$trees))
  } else if ( tolower(tree.dist) == "spr" ) {
    dmat <- as.matrix(phangorn::SPR.dist(simulated.samples$trees))
  } else {
    stop("Invalid input for argument 'tree.dist'")
  }
  
  res <- lapply(measures,function(this_measure) {
    eval(call(this_measure,simulated.samples=simulated.samples,dmat=dmat))
  })
  names(res) <- measures
  return(res)
}

#' Obtains estimates of the distribution of error and variance measures from simulated phylogenetic MCMC.
#'
#' @param simulated.samples An object of class simulatedPosterior (output of simulatePhylogeneticMCMC).
#' @param tree.dist The distance measure for trees (RF or SPR, only one).
#' @param measures The error or variance measure(s) (see details).
#' @param ess.methods The ESS calculation method(s) (see details).
#' @param return.ess Should the returned lists include the calculated ESS for each chain? 
#' @param verbose Should progress be printed to screen?
#' @return The first layer are the different ess methods used, the second the performance measures.
#' So $CMDSESS$RFDistanceToTrueMean contains the distribution of RF distances to the distribution of RF distances to the true mean tree when using CMDSESS for calculating ESS.
#' @export
#' @examples
#' countTrees(rmtree(100,10))
equivalentError <- function(simulated.samples,tree.dist=c("RF","SPR"),measures=c("treeProbSquaredError","treeVarianceSquaredError","splitProbSquaredError","distanceToTrueMRC"),ess.methods=c("fixedN",getESSMethods()),return.ess=TRUE,verbose=TRUE) {
  # recover()
  
  if ( !("simulatedPosterior" %in% class(simulated.samples) )) {
    stop("'simulated.samples' must be an object of class simulatedPosterior")
  }
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  # Start with the distance matrix for all unique topologies (allows us to reduce compute time)
  dmat <- matrix(nrow=ntrees,ncol=ntrees)
  if ( tolower(tree.dist) == "rf" ) {
    dmat <- as.matrix(phangorn::RF.dist(simulated.samples$trees))
  } else if ( tolower(tree.dist) == "spr" ) {
    dmat <- as.matrix(phangorn::SPR.dist(simulated.samples$trees))
  } else {
    stop("Invalid input for argument 'tree.dist'")
  }
  
  # get ESS for every method for every chain
  if ( verbose ) {
    cat("Computing ESS for all measures and all chain\n")
    pb <- txtProgressBar(max=nchains,style=3)
  }
  all_ess <- t(sapply(1:nchains,function(j){
    dm <- expandDistanceMatrix(dmat,simulated.samples$indices[,j])
    phy <- simulated.samples$trees[simulated.samples$indices[,j]]
    these_ess <- sapply(ess.methods,function(ess_method){
      eval(call(ess_method,dmat=dm,trees=phy))
    })

    if ( verbose ) {
      setTxtProgressBar(pb,j)
    }

    return(these_ess)
  }))
  
  if ( length(ess.methods) == 1 ) {
    all_ess <- t(all_ess)
  }
  
  colnames(all_ess) <- ess.methods
  
  # draw n_eff samples for each chain and each ESS calculated
  if ( verbose ) {
    cat("\nDrawing n_eff iid samples for all measures and all chains\n")
  }
  iid_n_samples_max <- round(max(all_ess))
  iid_samples <- lapply(ess.methods,function(this_method) {
    these_ess <- round(all_ess[,colnames(all_ess) == this_method])
    iid <- simulated.samples
    drawn_indices <- lapply(these_ess,function(n_eff){
      # Some ESS may be in [0,1], but we can't compute any measures without trees
      n_eff <- max(1,n_eff)
      drawn <- sample.int(ntrees,n_eff,replace=TRUE,prob=iid$probs)
      if (n_eff < iid_n_samples_max) {
        drawn <- c(drawn,rep(NA,iid_n_samples_max - n_eff))
      }
      return(drawn)
    })
    iid$indices <- do.call(cbind,drawn_indices)
    return(iid)
  })
  
  # compute the error/variance measures for the iid samples
  if ( verbose ) {
    cat("Computing all measures for all chains\n")
  }
  res <- lapply(iid_samples,function(iid) {
    these_measures <- lapply(measures,function(this_measure) {
      eval(call(this_measure,simulated.samples=iid,dmat=dmat))
    })
    names(these_measures) <- measures
    return(these_measures)
  })
  names(res) <- ess.methods
  
  if ( return.ess ) {
    for (i in 1:length(ess.methods)) {
      res[[i]]$ESS <- all_ess[,i]
    }
  }
  
  return(res)
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

MRC <- function(trees,probs) {
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
  
}

# Calculates d(MRC(posterior samples),MRC(true posterior)) for an object of class simulatedPosterior, where MRC(trees) is the MRC for those trees
# The distance measure is either RF distance or SPR distance
distanceToTrueMRC <- function(simulated.samples,tree.dist=c("RF","SPR"),...) {
  # recover()
  
  # ntrees <- length(simulated.samples$trees)
  # nchains <- dim(simulated.samples$indices)[2]
  # 
  # if ( length(true_summary_index) > 1 ) {
  #   warning("There is no unique true summary tree for the true posterior, taking first.")
  #   true_summary_index <- true_summary_index[1]
  # }
  # true_summary <- simulated.samples$trees[[true_summary_index]]
  # 
  # estimated_summary_indices <- sapply(1:nchains,function(j){
  #   est_probs <- sapply(1:ntrees,function(i){sum(simulated.samples$indices[,j] == i,na.rm=TRUE)})
  #   est_probs <- est_probs/sum(est_probs)
  #   idx <- summaryTreeIndex(dmat,est_probs,summary=summary)
  #   if ( length(idx) > 1 ) {
  #     warning(paste0("There is no unique true summary tree for chain ",j,", taking first."))
  #     idx <- idx[1]
  #   }
  #   return(idx)
  # })
  # 
  # d <- numeric(nchains)
  # if ( tolower(tree.dist) == "rf" ) {
  #   d <- phangorn::RF.dist(simulated.samples$trees[estimated_summary_indices],true_summary)
  # } else if ( tolower(tree.dist) == "spr" ) {
  #   d <- phangorn::SPR.dist(simulated.samples$trees[estimated_summary_indices],true_summary)
  # } else {
  #   stop("Invalid tree distance metric.")
  # }
  # 
  # return(d)
  return(-Inf)
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
