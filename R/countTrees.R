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