#' Counts occurrences of each tree in a set of trees, returns unique trees and their counts
#'
#' @param trees A multiPhylo object or list of phylo objects
#' @return A list. $trees contains all unique topologies, and $counts their counts
#' @export
#' @examples
#' countTrees(rmtree(100,10))
countTrees <- function(trees,rooted=FALSE) {
  # recover()
  if ( class(trees) != "multiPhylo" ) {
    if ( class(trees) == "list" && all(sapply(trees,class) == "phylo") ) {
      class(trees) <- "multiPhylo"
    } else {
      stop("Invalid input.")
    }
  }
  
  ntrees <- length(trees)
  res <- countTreesRecursive(trees = trees,
                             indices = 1:ntrees,
                             unique.indices = rep(0,ntrees),
                             counts = rep(0,ntrees),
                             idx = 0,
                             rooted = rooted)
  ordering <- order(res$counts,decreasing=TRUE)
  res$counts <- res$counts[ordering]
  res$unique.indices <- res$unique.indices[ordering]
  res$unique.indices <- res$unique.indices[res$counts > 0]
  res$counts <- res$counts[res$counts > 0]
  return(list(trees = trees[res$unique.indices], counts = res$counts))
}

countTreesRecursive <- function(trees,indices,unique.indices,counts,idx,rooted=FALSE) {
  # recover()
  ntrees <- length(trees)
  
  unique_indices <- unique.indices
  idx <- idx + 1
  
  dists <- suppressWarnings(phangorn::RF.dist(trees[indices],trees[[indices[1]]],rooted=rooted))
  
  tab <- table(dists)
  
  # Count (and pull off pile) trees matching this current tree
  matches_this_tree <- dists == 0
  
  counts[idx] <- sum(matches_this_tree)
  unique_indices[idx] <- indices[1]
  
  indices <- indices[!matches_this_tree]
  tab <- tab[names(tab) != "0"]
  dists <- dists[dists > 0]
  
  # Count (and pull off pile) any newly-discovered unique trees
  one_off_dists <- as.integer(names(which(tab == 1)))
  if ( length(one_off_dists) > 0 ) {
    # tau_i == tau_j implies d(tau_i,tau_ref) == d(tau_j,tau_ref)
    # So any unique distance to a given reference tree must be a unique tree
    is_one_off_dist <- dists %in% one_off_dists
    one_off_tree_indices <- indices[is_one_off_dist]
    
    # Remove trees matching this tree and any newly-found unique trees from consideration
    indices <- indices[!is_one_off_dist]
    tab <- tab[tab > 1]
    dists <- dists[!is_one_off_dist]
    
    # Count newly-found unique trees
    counts[idx + (1:length(one_off_tree_indices))] <- 1
    unique_indices[idx + (1:length(one_off_tree_indices))] <- one_off_tree_indices
    idx <- idx + length(one_off_tree_indices)
  }
  
  # If there are still uncounted trees, they are trees that do not match the current tree and are not unique
  # Now we handle them
  if ( length(tab) > 0 ) {
    for (i in 1:length(tab)) {
      target_distance <- as.integer(names(tab)[i])
      to_consider <- indices[dists == target_distance]
      tmp <- countTreesRecursive(trees=trees,indices=to_consider,unique.indices=unique_indices,counts=counts,idx=idx,rooted=rooted)
      unique_indices <- tmp$unique.indices
      counts <- tmp$counts
      idx <- tmp$idx
    }
  }
  
  return(list(unique.indices = unique_indices,
              counts = counts,
              idx = idx))
}
