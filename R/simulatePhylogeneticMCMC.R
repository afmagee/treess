#' Counts occurences of each tree in a set of trees, returns unique trees and their counts
#'
#' @param adjacency.graph An object of class treeAdjacencyGraph.
#' @param nchains Number of (independent) MCMC chains to run.
#' @param ngen Number of MCMC generations to run.
#' @param thin Samples will be thinned to one every thin generations.
#' @param move.weights Relative probability at each step of each move type (NNI,SPR).
#' @param verbose Whether to print summaries of MCMC performance.
#' @details No burnin is required, the chain is initialized from the stationary distribution.
#' @return Object of class simulatedPosterior, a list. 
#' $trees contains the unique trees in the target posterior (adjacency.graph$trees).
#' $probs contains the probabilities of the trees in the target posterior (adjacency.graph$probs).
#' $indices contains a matrix of the MCMC chains as the indices of the trees at each generation, chains in columns, generations in rows.
#' @export
#' @examples
#' data(DS3adj)
#' post <- simulatePhylogeneticMCMC(DS3adj,ngen=50,nchains=1)
#' post$trees[post$indices[,1]]
simulatePhylogeneticMCMC <- function(adjacency.graph,ngen=1000,nchains=100,thin=1,move.weights=c(2,1),verbose=TRUE) {
  # recover()
  
  if ( !("treeAdjacencyGraph" %in% class(adjacency.graph)) ) {
    stop("Argument 'adjacency.graph' must be of class treeAdjacencyGraph")
  }
  n_trees <- length(adjacency.graph$trees)
  n_tips <- length(adjacency.graph$trees[[1]]$tip.label)
  
  n_nni_neighbors <- 2*(n_tips-3)
  n_spr_neighbors <- 2*(n_tips-3)*(2*n_tips-7)
  
  use_nni <- ("nni" %in% tolower(names(adjacency.graph)))
  use_spr <- ("spr" %in% tolower(names(adjacency.graph)))
  
  n_tried_nni <- 0
  n_accepted_nni <- 0

  n_tried_spr <- 0
  n_accepted_spr <- 0
  
  move_weights <- move.weights
  moves <- c("NNI","SPR")
  if ( use_nni && !use_spr ) {
    move_weights[2] <- 0
  } else if ( use_spr && !use_nni ) {
    move_weights[1] <- 0
  } else if ( !(use_nni || use_spr) ) {
    stop("Invalid treeAdjacencyGraph: either $NNI or $SPR must be present")
  }
  move_probs <- move.weights/sum(move.weights)
  
  # We store all the trees as their indices
  mcmc <- matrix(nrow=ngen,ncol=nchains)
  
  if (verbose) {
    pb <- txtProgressBar(min=0,max=nchains,style=3)
  }
  
  for (j in 1:nchains) {
    current <- sample.int(n_trees,1,prob=adjacency.graph$probs)
    for (i in 1:ngen) {
      move_type <- sample(moves,1,prob=move_probs)
      
      proposed <- Inf
      
      if ( move_type == "NNI" ) {
        n_tried_nni <- n_tried_nni + 1
        
        # first we determine if we propose a tree in support (otherwise we would reject the tree and we're done)
        p_in_support <- length(adjacency.graph$NNI[[current]])/n_nni_neighbors
        if ( rbinom(1,1,p_in_support) ) {
          # We have proposed a tree in support, which one?
          proposed <- sample(adjacency.graph$NNI[[current]],1)
          # Do we accept the move or do we reject it?
          AR <- adjacency.graph$probs[proposed]/adjacency.graph$probs[current]
          if ( runif(1) < AR ) {
            n_accepted_nni <- n_accepted_nni + 1
          } else {
            # Reject
            proposed <- current
          }
        } else {
          proposed <- current
        }
      } else {
        n_tried_spr <- n_tried_spr + 1
        
        # first we determine if we propose a tree in support (otherwise we would reject the tree and we're done)
        p_in_support <- length(adjacency.graph$SPR[[current]])/n_spr_neighbors
        if ( rbinom(1,1,p_in_support) ) {
          # We have proposed a tree in support, which one?
          proposed <- sample(adjacency.graph$SPR[[current]],1)
          # Do we accept the move or do we reject it?
          AR <- adjacency.graph$probs[proposed]/adjacency.graph$probs[current]
          if ( runif(1) < AR ) {
            n_accepted_spr <- n_accepted_spr + 1
          } else {
            # Reject
            proposed <- current
          }
        } else {
          proposed <- current
        }
      }
      current <- proposed
      mcmc[i,j] <- current
    }
    if ( verbose ) {
      setTxtProgressBar(pb,j)
    }
  }
  if ( verbose ) {
    cat("\n")
    if ( use_nni ) {
      cat("NNI: #tried = ",n_tried_nni,"; %accepted = ",n_accepted_nni/n_tried_nni*100,"\n",sep="")
    }
    if ( use_spr ) {
      cat("SPR: #tried = ",n_tried_spr,"; %accepted = ",n_accepted_spr/n_tried_spr*100,"\n",sep="")
    }
  }
  mcmc <- mcmc[seq(1,ngen,thin),]
  
  res <- list(trees=adjacency.graph$trees,probs=adjacency.graph$probs,indices=mcmc)
  class(res) <- "simulatedPosterior"
  
  return(res)
}
