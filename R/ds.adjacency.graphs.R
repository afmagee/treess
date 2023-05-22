#' Compact summaries of posterior distributions.
#'
#' Specifically, these are sets of topologies (there are no branch lengths) for standard test datasets collectively referred to as DS1, DS3, and DS5.
#' These three datasets capture strongly peaked (DS3), flat (DS5), and multimodal (DS1) posterior distributions.
#' The sets of trees represent the set of NNI-connected subset of the "best trees" from Whidden et al. (2020).
#' The best trees are either all trees in the 95\% HPD set or the 4096 with the highest posterior support (whichever is smaller).
#' Here, only the largest subset of the "best trees" which form a connected graph (with respect to NNI tree moves) are kept.
#' 
#' @docType data
#' 
#' @usage data(ds.adjacency.graphs)
#' 
#' @format A list of objects of class treeAdjacencyGraph. 
#' This class is a list containing the phylogenies, their (renormalized) posterior probabilities, and the NNI connectivity graph linking the trees.
#' 
#' @keywords datasets
#' 
#' @references
#' Magee, Karcher, Matsen, and Minin (2023). "How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error." Bayesian Analysis, 1(1), 1-29.
#' 
#' Whidden et al. (2020) "Systematic exploration of the high likelihood set of phylogenetic tree topologies." Systematic Biology 69.2 (2020): 280-293.
#' 
#' @source Original DS analyses come from Whidden et al. (2020). NNI-connected subsets were determined by Magee et al. (2023).
#' 
#' @examples 
#' data(ds.adjacency.graphs)
#' \dontrun{
#' simulatePhylogeneticMCMC(ds.adjacency.graphs$ds1.adjacency.graph)
#' }
"ds.adjacency.graphs"