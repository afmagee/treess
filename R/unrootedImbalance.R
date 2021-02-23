#' Computes the Colless imbalance measure marginalized over all rootings of the tree.
#'
#' @param phy The tree as a class phylo object.
#' @return The imbalance measure.
#' @export
unrootedImbalance <- function(phy) {
  # recover()
  if ( ape::is.rooted(phy) ) {
    phy <- ape::unroot(phy)
  }
  
  ui <- NaN
  if ( phy$Nnode == (length(phy$tip.label) - 2) ) {
    ui <- unrootedImbalanceResolved(phy)
  } else {
    stop("Tree must be fully resolved.")
    # ui <- unrootedImbalanceUnesolved(phy)
  }
  
  return(ui)
}

unrootedImbalanceResolved <- function(phy) {
  # recover()
  
  ntax <- length(phy$tip.label)
  
  tip_descendants <- phangorn::Descendants(phy,(ntax+1):(ntax+phy$Nnode),type="tips")
  
  abs_diffs <- sapply(1:phy$Nnode,function(i){
    node <- i + ntax
    children <- phy$edge[phy$edge[,1] == node,2]
    s <- sapply(children,function(j){
      if ( j <= ntax ) {
        return(1)
      } else {
        return(length(tip_descendants[[j - ntax]]))
      }
    })
    if ( sum(s) < ntax ) { #unless we're at the root, we're missing all the taxa in the complement of the clade
      s <- c(s,ntax-sum(s))
    }
    orientations <- sapply(1:length(s),function(j){
      (2*s[j]-1) * (max(s[-j]) - min(s[-j]))
    })
    return(sum(orientations))
  })
  
  return(1/(2*ntax-3)*sum(abs_diffs))
}

# unrootedImbalanceUnresolved <- function(phy) {
#   ntax <- length(phy$tip.label)
#   
#   tip_descendants <- phangorn::Descendants(phy,(ntax+1):(ntax+phy$Nnode),type="tips")
#   
#   abs_diffs <- sapply(1:phy$Nnode,function(i){
#     node <- i + ntax
#     children <- phy$edge[phy$edge[,1] == node,2]
#     n <- sapply(children,function(j){
#       if ( j <= ntax ) {
#         return(1)
#       } else {
#         return(length(tip_descendants[[j - ntax]]))
#       }
#     })
#     if ( sum(n) < ntax ) { #unless we're at the root, we're missing all the taxa in the complement of the clade
#       n <- c(n,ntax-sum(n))
#     }
#     orientations <- sapply(1:length(n),function(i){
#       (2*n[i]-1) * (max(n[-i]) - min(n[-i]))
#     })
#     return(sum(orientations))
#   })
#   
#   return(1/(2*ntax-3)*sum(abs_diffs))
# }