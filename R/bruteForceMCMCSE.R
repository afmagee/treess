#' Obtains estimates of the distribution of error and variance measures from simulated phylogenetic MCMC.
#'
#' @param simulated.samples An object of class simulatedPosterior (output of simulatePhylogeneticMCMC).
#' @param tree.dist The distance measure for trees (RF or SPR, only one).
#' @param measures The error or variance measure(s) (see details).
#' @return A named list with the distribution for each measure.
#' @export
#' @examples
#' countTrees(rmtree(100,10))
bruteForceMCMCSE <- function(simulated.samples,tree.dist=c("RF","SPR"),measures=c("treeProbSquaredError","splitProbSquaredError","MRCSquaredError")) {
  # recover()
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  # Start with the distance matrix for all unique topologies (allows us to reduce compute time)
  dmat <- matrix(nrow=ntrees,ncol=ntrees)
  if ( tolower(tree.dist) == "rf" ) {
    dmat <- as.matrix(phangorn::RF.dist(simulated.samples$trees))
  } else if ( tolower(tree.dist) == "spr" ) {
    dmat <- as.matrix(phangorn::SPR.dist(simulated.samples$trees))
  } else {
    stop("Invalid input for argument 'tree.dist'")
  }
  
  res <- lapply(measures,function(this_measure) {
    eval(call(this_measure,simulated.samples=simulated.samples))
  })
  names(res) <- measures
  return(res)
}
