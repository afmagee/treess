#' Trees as vectors of splits.
#'
#' Transforms sample of trees into their coordinates in (reduced) RF space. Reduction is to the set of coordinates (splits/bipartitions) seen in at least one tree.
#'
#' @param trees list of trees or multiPhylo object
#' @return Output: a list. $coords is a matrix of trees as coordinates (see details), and $taxa is a list of taxa in each split in the matrix (ordered 1:n).
#' @details In the output$coords, every row is a tree, every column the presence/absence of a split.
#' @export
as.RFcoords <- function(trees) {
  splits <- trees2Coords(trees,namesplits=TRUE)
  split_taxa <- colnames(splits)
  split_taxa <- lapply(1:dim(splits)[2],function(i){
    strsplit(split_taxa[i],";")[[1]]
  })
  colnames(splits) <- NULL
  return(list(coords=splits,taxa=split_taxa))
}

#' Tree distance
#'
#' Calculates unrooted RF distance between a set of trees.
#'
#' @param coords matrix of trees as coordinates, as.RFcoords(trees)$coords.
#' @return Distances as an object of class dist.
#' @export
RF.from.coords <- function(coords) {
  if ( (length(dim(coords)) !=2) || any(!(as.numeric(coords) %in% c(0,1))) ) {
    stop("Input must be coordinate matrix")
  }
  return(dist(coords,method="manhattan"))
}

# Transforms sample of trees into their coordinates in (reduced) KF space.
# Reduction is to the set of coordinates/bipartitions seen in at least one tree
# Arguments: 
#   trees: list of trees or multiPhylo object
#   namesplits: should we name the splits in the matrix?
# Output: trees as coordinates
trees2Coords <- function(trees,namesplits=FALSE) {
  # recover()
  
  taxa <- trees[[1]]$tip.label
  
  ntax <- length(taxa)
  
  # Get master list of all splits
  all_splits <- as.matrix(phangorn::as.splits(trees))
  
  # Order alphabetically
  all_splits <- all_splits[,order(colnames(all_splits))]
  
  # Remove trivial splits
  trivial <- rowSums(all_splits) == 1 | rowSums(all_splits) == ntax  | rowSums(all_splits) == (ntax - 1)
  
  all_splits <- all_splits[!trivial,]
  
  # Polarize, our rule here is that all splits should include the first taxon
  to_polarize <- all_splits[,1] == 0
  
  all_splits[to_polarize,] <- -1 * (all_splits[to_polarize,] - 1)
  
  # Collapse to strings
  split_names <- c()
  if ( namesplits ) {
    taxa <- colnames(all_splits)
    split_names <- apply(all_splits,1,function(x){
      paste0(taxa[as.logical(x)],collapse=";")
    })
  }
  # all_splits <- apply(all_splits,1,paste0,collapse="")
  # This is faster. I don't know why, but it is
  all_splits <- apply(all_splits,1,function(x){paste0(as.raw(x),collapse="")})
  
  coords <- matrix(0,nrow=length(trees),ncol=length(all_splits))
  
  for (i in 1:length(trees)) {
    these_splits <- as.matrix(phangorn::as.splits(trees[[i]]))
    
    # alphabetize
    these_splits <- these_splits[,order(colnames(these_splits))]
    
    # remove trivial splits (only one taxon or all taxa)
    trivial <- rowSums(these_splits) == 1 | rowSums(these_splits) == ntax  | rowSums(these_splits) == (ntax - 1)
    
    these_splits <- these_splits[!trivial,]
    
    # Polarize, our rule here is that splits should be <= 50% 1s, and if 50% the first element should be a 0
    to_polarize <- these_splits[,1] == 0
    
    these_splits[to_polarize,] <- -1 * (these_splits[to_polarize,] - 1)
    
    # these_splits <- apply(these_splits,1,paste0,collapse="")
    # This is faster. I don't know why, but it is
    these_splits <- apply(these_splits,1,function(x){paste0(as.raw(x),collapse="")})
    
    # seen <- all_splits %in% these_splits
    
    coords[i,all_splits %in% these_splits] <- 1
    
  }
  
  if ( namesplits ) {
    colnames(coords) <- split_names
  }
  
  return(coords)
  
}