#' Computes quantities needed for calculating Monte Carlo error for trees. 
#'
#' @param simulated.samples An object of class simulatedPosterior (output of simulatePhylogeneticMCMC).
#' @param measures The error or variance measure(s) (see details).
#' @return A named list with the distribution for each measure. Same as a single element in \link{effectiveSizeEquivalentError}.
#' @export
#' @seealso \link{effectiveSizeEquivalentError}
bruteForceMCMCSE <- function(simulated.samples,measures=c("treeProbSquaredError","splitProbSquaredError","MRCSquaredError")) {
  # recover()
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  simulated.samples$coords <- trees2Coords(simulated.samples$trees)
  
  res <- lapply(measures,function(this_measure) {
    eval(call(this_measure,simulated.samples=simulated.samples))
  })
  names(res) <- measures
  return(res)
}
