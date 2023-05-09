# should consider making this a proper class

#' Trees as collections of subsplits.
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
getCladePairs <- function(trees,rooted=TRUE) {
  if ( !rooted ) {
    stop("Rooted trees only.")
  }
  
  if ( !all(ape::is.rooted.multiPhylo(trees)) ) {
    stop("Trees must be rooted")
  }
  if ( !all(ape::is.binary.multiPhylo(trees)) ) {
    stop("Trees must be bifurcating.")
  }
  
  # recover()
  ntaxa <- length(trees[[1]]$tip.label)
  
  trees <- checkRootedOption(trees,rooted)
  
  trees <- ape::.compressTipLabel(trees)
  
  clade_pairs <- lapply(trees,getRootedSubsplits,check.tree=FALSE)
  
  attr(clade_pairs,"labels") <- trees[[1]]$tip.label
  attr(clade_pairs,"rooted") <- TRUE
  
  return(clade_pairs)
}

#' Get subsplits from a tree, assuming it is rooted and bifurcating
#'
#' @param phy An ape-style tree
#' @param precompressed Whether or not the tree has already had its tip labels compressed by ape::.compressTipLabel
#' @param ref If the tip labels need compressing, this reference will be used, must be provided if precompressed == FALSE
#' @param check.tree Should rooting and bifurcation be checked? This should not be set to false unless trees have been pre-checked or bad things may happen!
#' @return A character vector of subsplits.
#' @details 
#' Subsplits are represented as |-separated lists of taxa, parent|child1|child2.
#' Each group is a comma separated, string-concatenated list of taxa.
#' @internal
getRootedCladePairs <- function(phy,check.tree=TRUE,precompressed=TRUE,ref=NA) {
  # recover()
  
  if ( !precompressed ) {
    if ( !is.character(ref) || length(ref) != ape::Ntip(phy) ) {
      stop("Must provide valid reference labels for compressing tip labels.")
    }
    phy <- ape::.compressTipLabel(list(phy),ref=ref)[[1]]
  }
  
  if ( check.tree ) {
    if ( !ape::is.rooted.phylo(phy) ) {
      stop("Tree must be rooted")
    }
    if ( !ape::is.binary(phy) ) {
      stop("Tree must be bifurcating.")
    }
  }

  phy <- ape::reorder.phylo(phy,order="postorder")
  
  descendants <- phangorn::Descendants(phy,type="tips")
  descendants <- lapply(descendants,sort)
  
  clade_pairs <- character(ape::Nnode(phy))
  for (i in 1:(dim(phy$edge)[1]/2)) {
    parent <- phy$edge[2*i,1]
    children <- phy$edge[c(2*(i-1)+1,2*i),2]
    clade_pairs[i] <- paste0(c(paste0(descendants[[parent]],collapse=","),
                             paste0(descendants[[children[1]]],collapse=","),
                             paste0(descendants[[children[2]]],collapse=",")),
                           collapse="|")
  }
  
  return(clade_pairs)
}
