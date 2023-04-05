# should consider making this a proper class

#' Trees as collections of splits.
#'
#' Gets list of splits/clades in trees, ignoring trivial splits/clades.
#'
#' @param trees List of trees or multiPhylo object
#' @param rooted If rooted=TRUE and trees are rooted, considers clades instead of splits. If rooted=FALSE, considers splits regardless.
#' @return A list, with attribute "labels", where element [[i]] is a vector of splits in trees[[i]]
#' @details 
#' Splits are represented as (comma separated, string-concatenated) lists of taxa.
#' Taxa are numbered according to the attribute "labels".
#' Splits are polarized such that they all exclude taxon trees[[1]]$tip.label[1].
#' @export
#' @seealso \link{splitProbs}, \link{nSplits}, \link{ape::prop.part}
treeSplits <- function(trees,rooted=FALSE) {
  ntaxa <- length(trees[[1]]$tip.label)
  
  trees <- checkRootedOption(trees,rooted)
  
  trees <- ape::.compressTipLabel(trees)
  trivial_length <- ntaxa - ifelse(rooted, 0, 1)
  if ( !rooted ) {
    trees <- ape::root.multiPhylo(trees,1,resolve.root=TRUE)
  }
  
  splits <- lapply(trees,function(phy){
    tmp <- unclass(ape::prop.part(phy,check.labels=FALSE))
    tmp <- tmp[lengths(tmp) < trivial_length]
    lens <- lengths(tmp)
    res <- sapply(tmp,paste,collapse=",")
    return(res)
  })
  
  attr(splits,"labels") <- trees[[1]]$tip.label

  return(splits)
}

#' Split probabilities from collections of splits.
#'
#' Computes split probabilities from output of \link{treeSplits}.
#'
#' @param tree.splits output of treeSplits
#' @param range optional argument to define which trees to compute probabilities for, defaults to all trees.
#' @return Named vector of split (or clade) probabilities.
#' @export
#' @seealso \link{perTreeSplits}
splitProbs <- function(tree.splits,range=NULL) {
  if ( !is.null(range) ) {
    tree.splits <- tree.splits[range[1]:range[2]]
  }
  res <- splitFrechetMean(tree.splits)
  names(res) <- attr(res,"labels")
  attr(res,"labels") <- NULL
  return(res)
}

# Computes split-based Frechet mean for a set of trees (fed in as splits)
# This is the same as the split probability
# Surprisingly fast, even with very large sets of splits
splitFrechetMean <- function(tree.splits) {
  all_splits <- unlist(tree.splits)
  tab <- table(all_splits)
  res <- as.numeric(tab/length(tree.splits))
  attr(res,"labels") <- names(tab)
  return(res)
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
  
  summand <- n * m * (1-m)
  
  sum(summand/(n - bessel))
  
}










# # Not nearly as efficient as the unweighted mean!
# # Here mainly for completeness
# weightedSplitFrechetMean <- function(tree.split.probs) {
#   all_split_names <- unique(names(unlist(tree.split.probs)))
#   all_splits <- numeric(length(all_split_names))
#   for (i in 1:length(tree.split.probs)) {
#     key <- match(names(tree.split.probs[[i]]),all_split_names)
#     all_splits[key] <- all_splits[key] + tree.split.probs[[i]]
#   }
#   names(all_splits) <- all_split_names
#   return(all_splits/length(tree.split.probs))
# }
# 
# splitFrechetVariance <- function(tree.splits,mean=NULL,bessel=TRUE) {
#   # recover()
#   m <- 0
#   n <- length(tree.splits)
#   if ( is.numeric(mean) ) {
#     m <- mean
#   } else {
#     m <- splitFrechetMean(tree.splits)
#   }
#   
#   # We make use of the following identity, true for any j in 1:length(x)
#   # sum((x - mean(x))^2) = sum((x - x[j])^2) - n*(x[j]^2 - mean(x)^2) - 2*n*mean(x)*(mean(x)-x[j])
#   # This decomposes additively, so we can mix splitwise and treewise evalutations
#   # The first component is all distances to the splits in reference tree j
#   summand1 <- sapply(2:length(tree.splits),function(i){
#     length(setdiff(tree.splits[[i]],tree.splits[[1]])) + length(setdiff(tree.splits[[1]],tree.splits[[i]]))
#   })
#   
#   # The rest of the requisite sumis computed splitwise between the reference tree and the mean
#   ref_tree_expanded <- numeric(length(m))
#   ref_tree_expanded[match(tree.splits[[1]],names(m))] <- 1
#   
#   summand2 <- - n * (ref_tree_expanded^2 - m^2) - 2 * n * m * (m - ref_tree_expanded)
#   
#   return((sum(summand1)+sum(summand2))/(length(tree.splits)-bessel))
#   
#   # summand <- sapply(1:length(tree.splits),function(i){
#   #   seen <- names(m) %in% tree.splits[[i]]
#   #   unseen_sum_sq <- sum(m[!seen]^2)
#   #   # key <- match(tree.splits[[i]],names(m))
#   #   seen_sum_sq <- sum((1 - m[seen])^2)
#   #   return(seen_sum_sq + unseen_sum_sq)
#   # })
#   # 
#   # return(sum(summand)/(length(tree.splits)-bessel))
#   
# }
# 
# weightedSplitFrechetVariance <- function(tree.split.probs,mean=NULL,bessel=TRUE) {
#   # recover()
#   m <- 0
#   n <- length(tree.split.probs)
#   if ( is.numeric(mean) ) {
#     m <- mean
#   } else {
#     m <- weightedSplitFrechetMean(tree.split.probs)
#   }
#   
#   # We make use of the following identity, true for any j in 1:length(x)
#   # sum((x - mean(x))^2) = sum((x - x[j])^2) - n*(x[j]^2 - mean(x)^2) - 2*n*mean(x)*(mean(x)-x[j])
#   # This decomposes additively, so we can mix splitwise and treewise evalutations
#   # The first component is all distances to the splits in reference tree j
#   ref_tree_splits <- names(tree.split.probs[[1]])
#   ref_tree_split_probs <- tree.split.probs[[1]]
#   one_ns <- 1:length(ref_tree_splits)
#   summand1 <- sapply(2:length(tree.split.probs),function(i){
#     key <- match(names(tree.split.probs[[i]]),ref_tree_splits)
#     unseen <- is.na(key)
#     seen <- !unseen
#     key <- key[seen]
#     seen_sum_sq <- sum((tree.split.probs[[i]][seen] - ref_tree_split_probs[key])^2)
#     unseen_sum_sq <- sum(tree.split.probs[[i]][unseen]^2) + sum(ref_tree_split_probs[!(one_ns %in% key)]^2)
#     return(seen_sum_sq+unseen_sum_sq)
#   })
#   
#   # The rest of the requisite sumis computed splitwise between the reference tree and the mean
#   ref_tree_expanded <- numeric(length(m))
#   key <- match(ref_tree_splits,names(m))
#   ref_tree_expanded[key] <- ref_tree_split_probs
#   
#   summand2 <- -n * (ref_tree_expanded^2 - m^2) - 2 * n * m * (m - ref_tree_expanded)
#   
#   res1 <- (sum(summand1)+sum(summand2))
#   
#   
#   cat(sum(summand1),sum(summand2),"\n",sep=",")
#   
#   return((sum(summand1)+sum(summand2))/(length(tree.split.probs)-bessel))
#   
#   summand <- sapply(1:length(tree.split.probs),function(i){
#     seen <- names(m) %in% names(tree.split.probs[[i]])
#     unseen_sum_sq <- sum(m[!seen]^2)
#     key <- match(names(tree.split.probs[[i]]),names(m))
#     seen_sum_sq <- sum((tree.split.probs[[i]] - m[key])^2)
#     return(seen_sum_sq + unseen_sum_sq)
#   })
#   
#   res2 <- sum(summand)
#   
#   return(res1 - res2)
#   
# }
