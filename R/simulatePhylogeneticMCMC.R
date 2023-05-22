#' NNI-based exploration of tree graph.
#' 
#' Takes in set of tree probabilities and an NNI adjacency graph on trees. 
#' Runs MCMC targeting the distribution on tree topologies.
#'
#' @param adjacency.graph An object of class treeAdjacencyGraph.
#' @param nchains Number of (independent) MCMC chains to run.
#' @param ngen Number of MCMC generations to run.
#' @param thin Samples will be thinned to one every thin generations.
#' @param verbose Whether to print summaries of MCMC performance.
#' @details No burnin is required, the chain is initialized from the stationary distribution.
#' @return Object of class simulatedPosterior, a list. 
#' $trees contains the unique trees in the target posterior (adjacency.graph$trees).
#' $probs contains the probabilities of the trees in the target posterior (adjacency.graph$probs).
#' $indices contains a matrix of the MCMC chains as the indices of the trees at each generation, chains in columns, generations in rows.
#' @export
#' @examples
#' post <- simulatePhylogeneticMCMC(constructAdjacencyGraph(phangorn::allTrees(5)),ngen=50,nchains=1)
#' post$trees[post$indices[,1]]
simulatePhylogeneticMCMC <- function(adjacency.graph,ngen=1000,nchains=100,thin=1,verbose=TRUE) {
  # recover()
  
  if ( !("treeAdjacencyGraph" %in% class(adjacency.graph)) ) {
    stop("Argument 'adjacency.graph' must be of class treeAdjacencyGraph")
  }
  if ( !("nni" %in% tolower(names(adjacency.graph))) ) {
    stop("Invalid treeAdjacencyGraph: $NNI must be present")
  }
  
  n_trees <- length(adjacency.graph$trees)
  n_tips <- length(adjacency.graph$trees[[1]]$tip.label)
  
  n_nni_neighbors <- 2*(n_tips-3)
  if ( attr(adjacency.graph,"rooted") == TRUE ) {
    n_nni_neighbors <- 2*(n_tips-2)
  }

  n_tried_nni <- 0
  n_accepted_nni <- 0
  
  # future-proofing for new move-types
  move_weights <- c(1)
  moves <- c("NNI")
  move_probs <- c(1)
  
  # We store all the trees as their indices
  mcmc <- matrix(nrow=ngen,ncol=nchains)
  
  if (verbose) {
    pb <- txtProgressBar(min=0,max=nchains,style=3)
  }
  
  for (j in 1:nchains) {
    current <- sample.int(n_trees,1,prob=adjacency.graph$probs)
    for (i in 1:ngen) {
      move_type <- "NNI" # future-proofing for new move-types
      
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
    cat("NNI: #tried = ",n_tried_nni,"; %accepted = ",n_accepted_nni/n_tried_nni*100,"\n",sep="")
  }
  mcmc <- mcmc[seq(1,ngen,thin),]
  
  res <- list(trees=adjacency.graph$trees,probs=adjacency.graph$probs,indices=mcmc)
  class(res) <- "simulatedPosterior"
  
  return(res)
}
