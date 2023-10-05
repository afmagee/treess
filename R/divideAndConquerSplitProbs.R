#' Approximate split probabilities from collections of splits.
#'
#' When there are many trees in the sample, the trees are large, and/or the posterior distribution is diffuse, the number of splits can be very large.
#' There may be a more tractable number of splits above some probability cutoff.
#' This function approximates that set of splits, and their probabilities.
#' The accuracy depends on several settings.
#'
#' @param tree.splits Output of \link{perTreeSplits}.
#' @param range Optional argument to define which trees to compute probabilities for, defaults to all trees.
#' @param cutoff Minimum split probability target. See details.
#' @param alpha Modifies target threshold from cutoff using a one-sided confidence interval with width alpha. Set to NULL if undesired. See details.
#' @param second.pass Gets correct probabilities for all splits identified at cost of extra computation. See details.
#' @return Named vector of split (or clade) probabilities.
#' @details 
#' The function operates on many batches of size sqrt(length(tree.splits)) to make the problem of identifying splits with probability >= cutoff easier.
#' 
#' The function first finds the set of target splits, which is all splits which are above a threshold probability on average across all batches.
#' If second.pass == FALSE, the function then returns these average probabilities.
#' These will not be exactly correct, and error is likely larger closer to the cutoff.
#' If second.pass == TRUE, the function then loops back over all batches and compute the actual probabilities of the splits.
#' This is slower, and may find more splits than second.pass == FALSE.
#' However, the algorithm is still not guaranteed to find all splits above the cutoff
#'
#' The thoroughness of the search for the set of splits depends on the actual threshold used, which is a function of alpha and the cutoff.
#' Specifically, during the first pass the threshold is either cutoff (if alpha == NULL) or the alpha x 100\% confidence interval for proportion = cutoff with sample size sqrt(n).
#' The smaller alpha is, the more splits which will be retained.
#' @export
#' @seealso \link{splitProbs}, \link{binomialProportionCI}
divideAndConquerSplitProbs <- function(tree.splits,range=NULL,cutoff=0.1,alpha=0.05,second.pass=TRUE) {
  if ( class(tree.splits) != "perTreeSplits" ) {
    stop("Cannot compute split probabilities for object not of class \"perTreeSplits\"")
  }
  if ( !is.null(range) ) {
    tree.splits <- tree.splits[range[1]:range[2]]
  }
  
  n <- length(tree.splits)
  b <- floor(sqrt(n))
  nbatches <- floor(n/b)
  batches <- lapply(1:nbatches,function(i){
    seq(i,n,nbatches)
  })
  
  threshold <- cutoff
  if ( is.numeric(alpha) ) {
    threshold <- binomialProportionCI(cutoff,b,"Jeffreys",ci.width=(1-2*alpha))[1,1]
  }
  
  # recover()
  
  probs <- NULL
  if ( second.pass ) {
    per_batch <- lapply(batches,function(idx){
      splitFrechetMean(tree.splits[idx])
    })
    
    per_batch_thinned <- lapply(per_batch,function(batch){
      return(batch[batch > threshold])
    })
    
    probs <- weightedSplitFrechetMean(per_batch_thinned,lengths(per_batch_thinned))
    
    split_names <- names(probs)
    probs <- numeric(length(probs))
    for (i in 1:length(per_batch)) {
      where <- fastmatch::fmatch(split_names,names(per_batch[[i]]))
      here <- !is.na(where)
      where <- where[here]
      probs[here] <- probs[here] + per_batch[[i]][where] * lengths(batches)[i]
    }
    probs <- probs/n
    names(probs) <- split_names
    probs <- probs[probs > cutoff]
  } else {
    per_batch <- lapply(batches,function(idx){
      probs <- treess:::splitFrechetMean(tree.splits[idx])
      return(probs[probs > threshold])
    })
    # cat("Had to check",mean(lengths(per_batch))," splits on average per batch.\n")
    
    probs <- weightedSplitFrechetMean(per_batch,lengths(per_batch))
    
    probs <- probs[probs > cutoff]
    # # We can compute 
    # split_names <- names(probs)
    # std_errors <- numeric(length(probs))
    # for (i in 1:length(per_batch)) {
    #   where <- match(split_names,names(per_batch[[i]]))
    #   here <- !is.na(where)
    #   where <- where[here]
    #   std_errors[here] <- std_errors[here] + (probs[here]-per_batch[[i]][where])^2
    #   not_here <- is.na(where)
    #   std_errors[not_here] <- std_errors[not_here] + probs[not_here]^2
    # }
    # std_errors <- sqrt(std_errors/nbatches)
    # return(cbind(probs,std_errors))
  }
  return(probs)
}

divideAndConquerSplitProbsFromFile <- function(file,burnin.frac,rooted,cutoff=0.1,alpha=0.05,second.pass=TRUE,format="beast",ref=NULL) {
  if ( !file.exists(file) ) {
    stop("Cannot find file.")
  }
  
  
  if ( !is.null(ref) ) {
    stop("Cannot currently handle labeling/relabeling of trees")
  }
  
  recover()
  text <- scan(file,what=character(),sep="\n")
  
  # TODO: if !is.null(ref), parse label block in NEXUS tree files, label trees, and relabel as needed
  if ( tolower(format) == "beast" ) {
    tree_lines <- grepl("tree STATE",text)
    if ( sum(tree_lines) == 0 ) {
      stop("Cannot find any trees.")
    }
    text <- text[tree_lines]
  } else {
    stop("Unrecognized tree file format.")
  }
  burnin <- round(burnin.frac * length(text))
  text <- text[-c(1:burnin)]
  
  n <- length(text)
  b <- floor(sqrt(n))
  nbatches <- floor(n/b)
  batches <- lapply(1:nbatches,function(i){
    seq(i,n,nbatches)
  })
  
  threshold <- cutoff
  if ( is.numeric(alpha) ) {
    threshold <- binomialProportionCI(cutoff,b,"Jeffreys",ci.width=(1-2*alpha))[1,1]
  }
  
  # recover()
  per_batch <- lapply(batches,function(idx){
    trees <- ape::read.tree(text=text[idx])
    probs <- treess:::splitFrechetMean(perTreeSplits(trees,rooted=rooted,ref=ref))
    return(probs[probs > threshold])
  })
  
  probs <- treess:::weightedSplitFrechetMean(per_batch,lengths(per_batch))
  first_pass_probs <- probs[probs > cutoff]
  if ( second.pass ) {
    split_names <- names(probs)
    probs <- numeric(length(probs))
    for (i in 1:length(per_batch)) {
      trees <- ape::read.tree(text=text[batches[[i]]])
      this_batch <- treess:::splitFrechetMean(perTreeSplits(trees,rooted=rooted,ref=ref))
      where <- fastmatch::fmatch(split_names,names(this_batch))
      here <- !is.na(where)
      where <- where[here]
      probs[here] <- probs[here] + this_batch[where] * lengths(batches)[i]
    }
    probs <- probs/n
    names(probs) <- split_names
    probs <- probs[probs > cutoff]
  } else {
    probs <- probs[probs > cutoff]
  }
  return(probs)
}