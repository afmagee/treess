#' Trees as collections of splits.
#'
#' Transforms sample of trees into their component splits (or clades).
#' Works like a version of ape::prop.part which returns per-tree (and which considers rooting).
#'
#' @param trees list of trees or multiPhylo object
#' @param rooted For rooted trees, if rooted=TRUE, clades are considered rather than splits. Otherwise trees are treated as unrooted.
#' @param ref Advanced option. Allows access to ref argument of \link{ape::.compressTipLabel}, so as to make splits/clades comparable across calls to different sets of trees on the same taxa.
#' @return Output an object of class perTreeSplits (a list of splits in each tree with attributes for the taxon labels and whether the tree is considered rooted or not).
#' @export
#' @seealso \link{splitProbs}, \link{nSplits}, \link{ape::prop.part}
perTreeSplits <- function(trees,rooted=FALSE,ref=NULL) {
  # recover()
  
  # We use rooted later assuming it behaves as 0 for FALSE, 1 for TRUE, so we check this now
  if ( !(rooted == TRUE || rooted == FALSE) ) {
    stop("Invalid option for rooted")
  }
  
  ntax <- ape::Ntip(trees[[1]])
  # short_taxa <- 1:ntax
  taxa <- trees[[1]]$tip.label
  if ( !is.null(ref) ) {
    taxa <- ref
  }

  trees <- checkRootedOption(trees,rooted=rooted)

  # Unrooted trees are rooted to a particular (but arbitrary) taxon
  # This allows us to look at clades as if they were splits, ie it polarizes all splits
  # This is much faster than polarizing after the fact
  # Rooted trees being analyzed for splits are re-rooted the same way
  class(trees) <- "multiPhylo"
  if ( !rooted ) {
    trees <- ape::root(ape::unroot(trees),taxa[1],resolve.root=TRUE)
    # trees <- lapply(trees,function(phy){
    #   ape::root(ape::unroot(phy),taxa[1],resolve.root=TRUE)
    # })
  }
  trees <- ape::.compressTipLabel(trees,ref=ref)
  
  per_tree <- lapply(trees,function(phy){
    # tmpphy <- phy
    # key <- match(phy$tip.label,taxa)
    # is_tip <- phy$edge[,2] <= ntax
    # # Taxon names are now meaningless, but all tip numbers match reference
    # phy$edge[is_tip,2] <- short_taxa[key[phy$edge[is_tip,2]]]
    # phy$tip.label <- taxa # We would need these tip labels
    # if (!all.equal(phy,tmpphy)){stop("Error handling trees")}
    splits <- ape::prop.part(phy)
    lens <- lengths(splits)
    # No trivial splits/clades
    splits <- splits[!(lens == 1 | lens == ntax | lens == ntax - 1 + rooted)]
    char_splits <- lapply(splits,function(split){
      split <- sort(split)
      return(paste0(split,collapse=","))
      # return(paste0(sort(taxa[split]),collapse=";"))
    })
    return(unlist(char_splits))
  })
  
  attr(per_tree,"labels") <- taxa
  attr(per_tree,"rooted") <- rooted
  class(per_tree) <- "perTreeSplits"
  
  return(per_tree)
}

#' Split probabilities from collections of splits.
#'
#' Computes split probabilities from output of \link{perTreeSplits}.
#' Works like ape::prop.clade to perTreeSplit's ape::prop.part.
#'
#' @param tree.splits output of perTreeSplits.
#' @param range optional argument to define which trees to compute probabilities for, defaults to all trees.
#' @return Named vector of split (or clade) probabilities.
#' @export
#' @seealso \link{perTreeSplits}
splitProbs <- function(tree.splits,range=NULL) {
  if ( class(tree.splits) != "perTreeSplits" ) {
    stop("Cannot compute split probabilities for object not of class \"perTreeSplits\"")
  }
  if ( !is.null(range) ) {
    tree.splits <- tree.splits[range[1]:range[2]]
  }
  splitFrechetMean(tree.splits)
}

#' Number of unique splits in a collection of \link{perTreeSplits}.
#'
#' @param tree.splits output of perTreeSplits.
#' @return Number of unique splits.
#' @export
nSplits <- function(tree.splits) {
  tree.splits <- unclass(tree.splits)
  length(unique(unlist(tree.splits)))
}

#' Split probabilities from collections of splits.
#'
#' Internal, optionless, error-checking-free version of \link{splitProbs}.
#'
#' @param tree.splits output of perTreeSplits.
#' @return Named vector of split (or clade) probabilities.
#' @keywords internal
#' @seealso \link{splitProbs}, \link{perTreeSplits}
splitFrechetMean <- function(tree.splits) {
  all_splits <- unlist(unclass(tree.splits))
  names(all_splits) <- NULL
  tab <- table(all_splits)
  res <- as.numeric(tab/length(tree.splits))
  names(res) <- names(tab)
  return(res)
}

#' Split Frechet variance from collections of splits.
#'
#' @param tree.splits output of perTreeSplits.
#' @param mean output of splitFrechetMean
#' @param bessel Should Bessel correction (divide by n-1 not n) be used? Default TRUE
#' @return Named vector of split (or clade) probabilities.
#' @details This function interchanges the order of summation.
#' The Frechet variance is the (average) squared deviation from each split vector to the mean.
#' Assuming all splits are collected into an nTrees x nSplits matrix X, this is sum_i sum_j (x_ij - mu_j)^2.
#' By interchanging the order of summation, the problem becomes a sum of variances of vectors of Bernoulli trials.
#' Thus, the Frechet split variance can be obtained directly from the mean and the number of trees.
#' @keywords internal
#' @seealso \link{splitProbs}, \link{perTreeSplits}, \link{splitFrechetMean}.
splitFrechetVariance <- function(tree.splits,mean=NULL,bessel=TRUE) {
  # recover()
  class(tree.splits) <- NULL
  p <- 0
  n <- length(tree.splits)
  if ( is.numeric(mean) ) {
    p <- mean
  } else {
    p <- splitFrechetMean(tree.splits)
  }
  
  summand <- length(tree.splits) * p * (1 - p)
  
  return((sum(summand))/(length(tree.splits)-bessel))
}

#' Useful quantities from collections of split Frechet means.
#'
#' Computes the grand mean from a collection of split Frechet means.
#' Can also return the sum of squared deviations to that grand mean and standard deviations of split frequencies.
#' Assumes equal weighting!
#'
#' @param mean.list list of output of \link{splitFrechetMean} or \link{splitProbs}.
#' @return List with grand mean ($mean), and optionally sum of squared deviations to the grand mean ($ssq) 
#' and optionally the standard deviations of all split frequencies.
#' @keywords internal
#' @seealso \link{splitProbs}, \link{perTreeSplits}
perTreeSplitsConvergenceHelper <- function(mean.list,weights=1,return.ssq=FALSE,return.sdsf=FALSE) {
  # recover()
  
  all_split_names <- unique(unlist(lapply(mean.list,names)))
  mat <- matrix(0,nrow=length(mean.list),ncol=length(all_split_names))
  for (i in 1:length(mean.list)) {
    key <- fastmatch::fmatch(names(mean.list[[i]]),all_split_names)
    mat[i,key] <- mean.list[[i]]
  }
  m <- colMeans(mat)
  names(m) <- all_split_names
  
  res <- list(mean=m)
  
  if (return.ssq) {
    ssq <- apply(mat,1,function(x){
      sum((x - m)^2)
    })
    res$ssq <- sum(ssq)
  }
  
  if (return.sdsf) {
    sdsf <- apply(mat,2,sd)
    res$sdsf <- sdsf
  }
  
  return(res)
  
}

#' Split Frechet variance from collections of splits which came from tree samples with uneven lengths.
#' 
#' For tree samples of even length, the unweighted mean is faster and should be used.
#'
#' @param split.prob.list List of outputs of perTreeSplits.
#' @param w Weights (e.g., lengths of each set of trees)
#' @return Named vector of split (or clade) probabilities.
#' @keywords internal
weightedSplitFrechetMean <- function(split.prob.list,w) {
  all_split_names <- unique(names(unlist(split.prob.list)))
  all_splits <- numeric(length(all_split_names))
  for (i in 1:length(split.prob.list)) {
    key <- fastmatch::fmatch(names(split.prob.list[[i]]),all_split_names)
    all_splits[key] <- all_splits[key] + split.prob.list[[i]] * w[i]
  }
  names(all_splits) <- all_split_names
  return(all_splits/sum(w))
}

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

# weightedSplitFrechetVariance <- function(tree.split.probs,mean=NULL,bessel=TRUE) {
#   recover()
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

###############
# Functions to make perTreeSplits a simple class
###############
# should maybe allow conversion to prop.part class?
#' @export
print.perTreeSplits <- function(x,...) {
  cat(length(x), "phylogenies processed into",ifelse(attr(x,"rooted"),"clades","splits"))
}

#' @export
print.treeSplits <- function(x,...) {
  cat("A phylogeny processed into",ifelse(attr(x,"rooted"),"clades","splits"))
}

#' @export
"[[.perTreeSplits" <- function(x, i) {
  rooted <- attr(x,"rooted")
  labels <- attr(x,"labels")
  class(x) <- NULL
  s <- x[[i]]
  attr(s,"rooted") <- rooted
  attr(s,"labels") <- labels
  class(s) <- "treeSplits"
  return(s)
}

#' @export
`$.perTreeSplits` <- function(x, name) {
  x[[name]]
}

#' @export
"[.perTreeSplits" <- function(x, i) {
  rooted <- attr(x,"rooted")
  labels <- attr(x,"labels")
  class(x) <- NULL
  s <- x[i]
  attr(s,"rooted") <- rooted
  attr(s,"labels") <- labels
  class(s) <- "perTreeSplits"
  return(s)
}

#' @export
str.perTreeSplits <- function(object,...)
{
  cat("Object of class \"perTreeSplits\"\n")
  str(unclass(object),...)
}

#' @export
c.perTreeSplits <- function(...,recursive=FALSE) {
  res <- list(...)
  classes <- lapply(res,class)
  if ( !all(classes == "perTreeSplits") ) {
    stop("Cannot concatenate provided objects.")
  }
  
  all_compressed_labels <- lapply(res,function(x){
    paste0(attr(x,"labels"),collapse=",")
  })
  
  if ( length(unique(all_compressed_labels)) > 1 ) {
    stop("Cannot combine \"perTreeSplits\" objects based on different reference labels.")
  }
  
  all_rooted_options <- lapply(res,function(x){
    attr(x,"rooted")
  })
  if ( length(unique(all_rooted_options)) > 1 ) {
    stop("Cannot combine \"perTreeSplits\" objects when some are rooted and others are unrooted.")
  }
  
  rooted <- attr(res[[1]],"rooted")
  labels <- attr(res[[2]],"labels")
  
  res <- lapply(res,function(x){
    x <- unclass(x)
    return(x[1:length(x)])
  })
  res <- do.call(c,res)
  
  attr(res,"rooted") <- rooted
  attr(res,"labels") <- labels
  class(res) <- "perTreeSplits"
  
  return(res)
}
