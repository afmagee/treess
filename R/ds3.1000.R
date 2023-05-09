#' An MCMC sample of 1000 trees.
#'
#' Specifically, this is a set of 1000 topologies (there are no branch lengths) produced by using \link{simulatePhylogeneticMCMC} on the standard test dataset known as DS3.
#' 
#' @docType data
#' 
#' @usage data(ds3.1000)
#' 
#' @format An ape-style multiPhylo object (without branch lengths).
#' 
#' @keywords datasets
#' 
#' @references
#' Magee, Karcher, Matsen, and Minin (2023). "How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error." Bayesian Analysis, 1(1), 1-29.
#' Whidden et al. (2020) "Systematic exploration of the high likelihood set of phylogenetic tree topologies." Systematic Biology 69.2 (2020): 280-293.
#' 
#' @source Whidden et al. (2020)
#' 
#' @examples 
#' \dontrun{
#' data(ds3.1000)
#' }