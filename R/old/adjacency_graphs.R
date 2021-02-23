#' Counts occurences of each tree in a set of trees, returns unique trees and their counts
#'
#' @param trees A multiPhylo object or list of phylo objects
#' @details Allowable methods are "pseudoPSRF", "totalDistancePSRF", "foldedRankMedioidPSRF", and "medianDistanceWilcoxTest"
#' @return A list. $trees contains all unique topologies, and $counts their counts
#' @export
#' @examples
#' countTrees(rmtree(100,10))
countTrees <- function(trees) {
  # recover()
  remaining_trees <- trees
  
  utrees <- vector("list")
  counts <- numeric()
  
  iter <- 0
  while( length(remaining_trees) > 0 ){
    iter <- iter + 1
    
    # Add this topology
    utrees[[iter]] <- remaining_trees[[1]]
    
    # count matches
    dists <- suppressWarnings(phangorn::RF.dist(remaining_trees,remaining_trees[[1]]))
    counts <- c(counts,sum(dists == 0))
    
    # Remove matches
    remaining_trees <- remaining_trees[!(dists == 0)]
  }
  
  class(utrees) <- "multiPhylo"
  
  return(list(trees=utrees,counts=counts))
}

#' Computes the NNI adjacency graph for a set of (unique) tree topologies.
#'
#' @param trees A multiPhylo object or list of phylo objects
#' @param weights The probabilities for each tree (need not sum to one)
#' @param trim Should we trim to the largest connected subset of trees in the posterior?
#' @return A list. $trees contains all topologies, $probs their probabilities, $NNI the NNI connectivity graph.
#' @export
#' @examples
#' constructAdjacencyGraph(rmtree(100,5))
constructAdjacencyGraph <- function(trees,weights=rep(1,length(trees)),trim=FALSE) {
  # recover()
  
  res <- list(trees=trees,probs=weights/sum(weights))
  
  n <- length(trees)
  
  # The RF distance is in general a lower bound on the NNI distance
  # However, for the special case of 1 NNI move, an RF distance of 2 and an NNI distance of 1 are equivalent
  dmat <- as.matrix(phangorn::RF.dist(trees))
  
  # Error checking
  if ( any(dmat %% 2 != 0) ) {
    stop("Provided trees contain polytomies. Please provide fully resolved trees.")
  }
  if ( any(dmat[upper.tri(dmat)] == 0) ) {
    stop("Trees must contain only unique trees.")
  }
  
  # Which trees are neighbors? 
  nni_adjacencies <- lapply(1:n,function(i){
    which(dmat[i,] == 2)
  })
  
  res$NNI <- nni_adjacencies
  
  # Handling of potentially nonconnected posteriors
  if (!requireNamespace("igraph",quietly=TRUE)) {
    warning("Cannot confirm that posterior distribution is connected without package igraph. If not connected problems will result in simulatePhylogeneticMCMC. Continue at own peril.")
    if ( trim ) {
      stop("Cannot trim posterior to connected subset without package igraph. Please install igraph to continue.")
    }
  } else {
    adj_mat <- dmat
    adj_mat[adj_mat > 2] <- Inf
    adj_mat[is.finite(adj_mat)] <- 1
    adj_mat[is.infinite(adj_mat)] <- 0
    graph <- igraph::graph_from_adjacency_matrix(adj_mat)
    
    if ( !igraph::is_connected(graph) ) {
      # Find largest connected subset
      comp <- igraph::components(graph)
      biggest <- which.max(comp$csize)
      
      # Prune out nonconnected trees
      keep <- which(comp$membership == biggest)
      trees <- trees[keep]
      dmat <- dmat[keep,keep]
      n <- length(trees)
      nni_adjacencies <- lapply(1:n,function(i){
        which(dmat[i,] == 2)
      })
      res$trees <- trees
      res$probs <- weights[keep]/sum(weights[keep])
      res$NNI <- nni_adjacencies
      
      # Report to user
      cat(paste0("Provided trees (",sum(comp$csize)," trees total) not connected, trimming down to largest connected subset (",n," trees total)\n"))
    }
  }
  
  class(res) <- "treeAdjacencyGraph"
  
  return(res)
}
