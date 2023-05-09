#' Fit Conditional Clade Distribution to a collection of trees.
#'
#' A Conditional Clade Distribution is a distribution on phylogenies.
#' It can produce better estimates of tree probabilities than simply counting trees in the posterior.
#'
#' @param trees List of trees or multiPhylo object
#' @param weights Optional vector of weights for the trees.
#' @return An SBN object.
#' @export
fitCCD <- function(trees,weights=1) {
  ntrees <- length(trees)
  ntaxa <- ape::Ntip(trees[[1]])
  
  unweighted <- FALSE
  if ( length(weights) < ntrees || !is.numeric(weights) ) {
    weights <- rep(1,ntrees)
    unweighted <- TRUE
  }
  
  subsplit_list <- treeSubsplits(trees,TRUE)
  
  # recover()
  # TODO this is not being done efficiently and will likely not scale well
  subsplits <- NULL
  parent_clades <- NULL
  parents_and_children <- NULL
  if (unweighted) {
    subsplits <- table(unlist(subsplit_list))
    weights <- subsplits
    subsplits <- names(subsplits)
  } else {
    nsubsplits_per <- length(subsplit_list[[1]])
    weights <- unlist(lapply(weights,function(w){rep(w,nsubsplits_per)}))
    subsplits <- unlist(subsplit_list)
    tmp_df <- data.frame(subsplits=subsplits,weights=weights)
    tmp <- aggregate(weights ~ subsplits,data=tmp_df,FUN=sum)
    subsplits <- tmp$subsplits
    weights <- tmp$weights
  }
  parents_and_children <- strsplit(subsplits,split="|",fixed=TRUE)
  parent_clades <- c(as.character(1:ntaxa),
                     unique(sapply(parents_and_children,function(x){x[[1]]})))
  
  conditional_probabilities <- vector("list",length=length(parent_clades))
  children <- vector("list",length=length(parent_clades))
  for (i in 1:length(parents_and_children)) {
    parent_index <- which(parent_clades == parents_and_children[[i]][1])
    child_indices <- sort(c(which(parent_clades == parents_and_children[[i]][2]),
                            which(parent_clades == parents_and_children[[i]][3])))
    child <- complex(real=child_indices[1],imaginary=child_indices[2])
    
    if ( (child %in% children[[parent_index]]) ) {
      where <- which(children[[parent_index]] == child)
      conditional_probabilities[[where]] <- conditional_probabilities[[where]] + weights[i]
    } else {
      children[[parent_index]] <- c(children[[parent_index]],child)
      conditional_probabilities[[parent_index]] <- c(conditional_probabilities[[parent_index]],weights[i])
    }
  }
  
  for (i in 1:length(conditional_probabilities)) {
    conditional_probabilities[[i]] <- log(conditional_probabilities[[i]]/sum(conditional_probabilities[[i]]))
  }
  
  ccd <- list(clades=parent_clades,
              child.splits=children,
              cpd=conditional_probabilities)
  attr(ccd,"labels") <- attr(subsplit_list,"labels")
  attr(ccd,"rooted") <- TRUE
  
  return(ccd)
}

#' Computes the density of a phylogeny from an SBN.
#'
#' @param phy An ape-style phylogeny
#' @param CCD The Conditional Clade Distribution.
#' @param log Should the log-probability be returned? (FALSE is not recommended.)
#' @return the (log-)probability mass
#' @export
dCCD <- function(phy,CCD,log=TRUE) {
  if ( attr(CCD,"rooted") == TRUE && !ape::is.rooted.phylo(phy) ) {
    stop("CCD is rooted but provided phylogeny is not.")
  } else if ( attr(CCD,"rooted") == FALSE && ape::is.rooted.phylo(phy) ) {
    stop("CCD is unrooted but provided phylogeny is rooted.")
  }
  
  subsplits <- NULL
  if ( attr(CCD,"rooted") ) {
    subsplits <- getRootedSubsplits(phy,check.tree=TRUE,precompressed=FALSE,ref=attr(CCD,"labels"))
  } else {
    stop("Unrooted CCDs not supported.")
  }
  
  # recover()
  
  logd <- 0.0
  parents_and_children <- strsplit(subsplits,split="|",fixed=TRUE)
  for (i in 1:length(subsplits)) {
    parent_child <- strsplit(subsplits[[i]],split="|",fixed=TRUE)[[1]]
    
    parent_index <- NULL
    if (parents_and_children[[i]][1] %in% CCD$clades) {
      parent_index <- which(CCD$clades == parents_and_children[[i]][1])
    } else {
      logd <- -Inf
      break
    }
    
    child_indices <- sort(c(which(CCD$clades == parents_and_children[[i]][2]),
                            which(CCD$clades == parents_and_children[[i]][3])))
    child <- complex(real=child_indices[1],imaginary=child_indices[2])
    
    if (child %in% CCD$child.splits[[parent_index]]) {
      logd <- logd + CCD$cpd[[parent_index]][which(CCD$child.splits[[parent_index]] == child)]
    } else {
      logd <- -Inf
      break
    }
  }
  names(logd) <- NULL
  
  if (!log) {
    logd <- exp(logd)
  }
  
  return(logd)
}
