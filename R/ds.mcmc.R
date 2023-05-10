#' MCMC samples of 1000 trees each.
#'
#' Specifically, these are sets of 1000 topologies (there are no branch lengths) produced by using \link{simulatePhylogeneticMCMC} on several standard test datasets collectively referred to as DS1, DS3, and DS5.
#' These three datasets capture strongly peaked (DS3), flat (DS5), and multimodal (DS1) posterior distributions.
#' For each, three chains were run resulting in low, medium, and high ESS values
#' 
#' @docType data
#' 
#' @usage data(ds.mcmc)
#' 
#' @format A list of ape-style multiPhylo objects (without branch lengths).
#' 
#' @keywords datasets
#' 
#' @references
#' Magee, Karcher, Matsen, and Minin (2023). "How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error." Bayesian Analysis, 1(1), 1-29.
#' 
#' Whidden et al. (2020) "Systematic exploration of the high likelihood set of phylogenetic tree topologies." Systematic Biology 69.2 (2020): 280-293.
#' 
#' @source Original DS analyses come from Whidden et al. (2020).
#' 
#' @examples 
#' data(ds.mcmc)
#' \dontrun{
#' treess(ds.mcmc,phangorn::RF.dist)
#' }
"ds.mcmc"