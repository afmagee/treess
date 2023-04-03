#' Trees as collections of splits.
#'
#' Transforms sample of trees into their component splits (or clades).
#' Works like a version of ape::prop.part which returns per-tree (and which considers rooting).
#'
#' @param trees list of trees or multiPhylo object
#' @param useClades for rooted trees, if useClades=TRUE, clades are considered rather than splits. Otherwise trees are treated as unrooted.
#' @param return.names should we return the taxa names?
#' @return Output: a list, see details
#' @details If return.names=TRUE, the list has two elements: $splits is a list where each element is all the splits(/clades) in a tree, and $taxa is the taxa ordered to match the splits. 
#' If return.names=FALSE, the returned list is simply the contents of $splits.
#' In the splits, the taxa are numbered (as in ape::prop.part) 
#' @export
#' @seealso \link{splitProbs}, \link{nSplits}, \link{ape::prop.part}
perTreeSplits <- function(trees,useClades=FALSE,return.names=FALSE) {
  # recover()
  
  # trees <- ape::rmtree(10,7)
  
  ntax <- ape::Ntip(trees[[1]])
  short_taxa <- 1:ntax
  taxa <- trees[[1]]$tip.label
  
  if ( useClades && !ape::is.rooted(trees[[1]]) ) {
    stop("Trees are not rooted, cannot examine clades.")
  }
  
  # Unrooted trees are rooted to a particular (but arbitrary) taxon
  # This allows us to look at clades as if they were splits, ie it polarizes all splits
  # This is much faster than polarizing after the fact
  class(trees) <- "multiPhylo"
  if ( !useClades ) {
    trees <- ape::root(ape::unroot(trees),taxa[1],resolve.root=TRUE)
    # trees <- lapply(trees,function(phy){
    #   ape::root(ape::unroot(phy),taxa[1],resolve.root=TRUE)
    # })
  }
  
  # Trivial splits and trivial clades are slightly different
  trivial_lens <- c(1,ntax)
  if ( !useClades ) {
    trivial_lens <- c(1,ntax,ntax-1)
  }
  
  per_tree <- lapply(trees,function(phy){
    # tmpphy <- phy
    key <- match(phy$tip.label,taxa)
    is_tip <- phy$edge[,2] <= ntax
    # Taxon names are now meaningless, but all tip numbers match reference
    phy$edge[is_tip,2] <- short_taxa[key[phy$edge[is_tip,2]]]
    # phy$tip.label <- taxa # We would need these tip labels
    # if (!all.equal(phy,tmpphy)){stop("Error handling trees")}
    splits <- ape::prop.part(phy)
    lens <- lengths(splits)
    # No trivial splits
    splits <- splits[!(lens == 1 | lens == ntax | lens == ntax - 1 + useClades)]
    char_splits <- lapply(splits,function(split){
      split <- sort(split)
      return(paste0(split,collapse=";"))
      # return(paste0(sort(taxa[split]),collapse=";"))
    })
    return(unlist(char_splits))
  })
  
  # length(unique(unlist(per_tree)))
  # XX <- ape::prop.part(trees)
  # length(XX)
  # sort(unique(unlist(per_tree)))
  # sort(unlist(lapply(XX,function(xx){paste0(sort(attributes(XX)$labels[xx]),collapse=";")})))
    
  if ( return.names ) {
    return(list(splits=per_tree,taxa=taxa))
  } else {
    return(per_tree)
  }
  
  # per_tree <- lapply(trees,function(phy){
  #   key <- match(phy$tip.label,taxa)
  #   is_tip <- phy$edge[,1] <= ntax
  #   # Taxon names are now meaningless, but all tip numbers match reference
  #   phy$edge[is_tip,1] <- short_taxa[key[phy$edge[is_tip,1]]]
  #   splits <- ape::prop.part(phy)
  #   lens <- lengths(splits)
  #   splits <- splits[!(lens == 1 | lens == ntax)]
  #   char_splits <- lapply(splits,function(split){
  #     split <- sort(split)
  #     # Polarize
  #     if ( split[1] != 1) {
  #       split <- short_taxa[!(short_taxa %in% split)]
  #     }
  #     return(paste0(split,collapse=";"))
  #   })
  #   return(unlist(char_splits))
  # })
  
}

# should maybe allow conversion to prop.part class?
# should maybe make this a proper class

#' Split probabilities from collections of splits.
#'
#' Computes split probabilities from output of \link{perTreeSplits}.
#' Works like ape::prop.clade to perTreeSplit's ape::prop.part.
#'
#' @param tree.splits output of perTreeSplits.
#' @param range optional argument to define which trees to compute probabilities for, defaults to all trees.
#' @param return.names should we return the taxa names?
#' @return Mamed vector of split (or clade) probabilities.
#' @export
#' @seealso \link{perTreeSplits}
splitProbs <- function(tree.splits,range=NULL) {
  if ( !is.null(range) ) {
    tree.splits <- tree.splits[range[1]:range[2]]
  }
  splitFrechetMean(tree.splits)
}

nSplits <- function(tree.splits) {
  if( "taxa" %in% names(tree.splits) ) {
    tree.splits <- tree.splits$splits
  }
  length(unique(unlist(tree.splits)))
}

# Very fast
splitFrechetMean <- function(tree.splits) {
  all_splits <- unlist(tree.splits)
  names(all_splits) <- NULL
  tab <- table(all_splits)
  res <- as.numeric(tab/length(tree.splits))
  names(res) <- names(tab)
  return(res)
}

# Not nearly as efficient as the unweighted mean!
# Here mainly for completeness
weightedSplitFrechetMean <- function(tree.split.probs) {
  all_split_names <- unique(names(unlist(tree.split.probs)))
  all_splits <- numeric(length(all_split_names))
  for (i in 1:length(tree.split.probs)) {
    key <- match(names(tree.split.probs[[i]]),all_split_names)
    all_splits[key] <- all_splits[key] + tree.split.probs[[i]]
  }
  names(all_splits) <- all_split_names
  return(all_splits/length(tree.split.probs))
}


splitFrechetVariance <- function(tree.splits,mean=NULL,bessel=TRUE) {
  # recover()
  m <- 0
  n <- length(tree.splits)
  if ( is.numeric(mean) ) {
    m <- mean
  } else {
    m <- splitFrechetMean(tree.splits)
  }
  
  # We make use of the following identity, true for any j in 1:length(x)
  # sum((x - mean(x))^2) = sum((x - x[j])^2) - n*(x[j]^2 - mean(x)^2) - 2*n*mean(x)*(mean(x)-x[j])
  # This decomposes additively, so we can mix splitwise and treewise evalutations
  # The first component is all distances to the splits in reference tree j
  summand1 <- sapply(2:length(tree.splits),function(i){
    length(setdiff(tree.splits[[i]],tree.splits[[1]])) + length(setdiff(tree.splits[[1]],tree.splits[[i]]))
  })
  
  # The rest of the requisite sumis computed splitwise between the reference tree and the mean
  ref_tree_expanded <- numeric(length(m))
  ref_tree_expanded[match(tree.splits[[1]],names(m))] <- 1
  
  summand2 <- - n * (ref_tree_expanded^2 - m^2) - 2 * n * m * (m - ref_tree_expanded)
  
  return((sum(summand1)+sum(summand2))/(length(tree.splits)-bessel))
  
  # summand <- sapply(1:length(tree.splits),function(i){
  #   seen <- names(m) %in% tree.splits[[i]]
  #   unseen_sum_sq <- sum(m[!seen]^2)
  #   # key <- match(tree.splits[[i]],names(m))
  #   seen_sum_sq <- sum((1 - m[seen])^2)
  #   return(seen_sum_sq + unseen_sum_sq)
  # })
  # 
  # return(sum(summand)/(length(tree.splits)-bessel))
  
}

weightedSplitFrechetVariance <- function(tree.split.probs,mean=NULL,bessel=TRUE) {
  recover()
  m <- 0
  n <- length(tree.split.probs)
  if ( is.numeric(mean) ) {
    m <- mean
  } else {
    m <- weightedSplitFrechetMean(tree.split.probs)
  }
  
  # We make use of the following identity, true for any j in 1:length(x)
  # sum((x - mean(x))^2) = sum((x - x[j])^2) - n*(x[j]^2 - mean(x)^2) - 2*n*mean(x)*(mean(x)-x[j])
  # This decomposes additively, so we can mix splitwise and treewise evalutations
  # The first component is all distances to the splits in reference tree j
  ref_tree_splits <- names(tree.split.probs[[1]])
  ref_tree_split_probs <- tree.split.probs[[1]]
  one_ns <- 1:length(ref_tree_splits)
  summand1 <- sapply(2:length(tree.split.probs),function(i){
    key <- match(names(tree.split.probs[[i]]),ref_tree_splits)
    unseen <- is.na(key)
    seen <- !unseen
    key <- key[seen]
    seen_sum_sq <- sum((tree.split.probs[[i]][seen] - ref_tree_split_probs[key])^2)
    unseen_sum_sq <- sum(tree.split.probs[[i]][unseen]^2) + sum(ref_tree_split_probs[!(one_ns %in% key)]^2)
    return(seen_sum_sq+unseen_sum_sq)
  })
  
  # The rest of the requisite sumis computed splitwise between the reference tree and the mean
  ref_tree_expanded <- numeric(length(m))
  key <- match(ref_tree_splits,names(m))
  ref_tree_expanded[key] <- ref_tree_split_probs
  
  summand2 <- -n * (ref_tree_expanded^2 - m^2) - 2 * n * m * (m - ref_tree_expanded)
  
  res1 <- (sum(summand1)+sum(summand2))
  

  cat(sum(summand1),sum(summand2),"\n",sep=",")

  return((sum(summand1)+sum(summand2))/(length(tree.split.probs)-bessel))
  
  summand <- sapply(1:length(tree.split.probs),function(i){
    seen <- names(m) %in% names(tree.split.probs[[i]])
    unseen_sum_sq <- sum(m[!seen]^2)
    key <- match(names(tree.split.probs[[i]]),names(m))
    seen_sum_sq <- sum((tree.split.probs[[i]] - m[key])^2)
    return(seen_sum_sq + unseen_sum_sq)
  })
  
  res2 <- sum(summand)
  
  return(res1 - res2)
  
}
