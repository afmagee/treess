#' Estimates distribution of error if we drew ESS trees IID from the (known) posterior distribution.
#' 
#' For one or more ESS methods, this function loops over chains for a pseudo-MCMC run and for each chain draws round(ESS) samples iid from the posterior distribution.
#' Then it computes a variety of summary measures useful for determining whether or not the ESS method works.
#' 
#' @param simulated.samples An object of class simulatedPosterior (output of \link{simulatePhylogeneticMCMC}).
#' @param tree.dist The distance measure for trees, only used for ESS computations (RF or SPR, use only one, default RF).
#' @param measures The error or variance measure(s) (see details).
#' @param ess.methods The ESS calculation method(s) (see details).
#' @param observed.trees.only If TRUE, restricts the drawing of iid trees to only draw from the portion of the posterior seen in the MCMC chains. Otherwise can draw any tree in the posterior.
#' @param return.ess Should the returned lists include the calculated ESS for each chain? 
#' @param verbose Should progress be printed to screen?
#' @return The first layer are the different ESS methods used, the second the performance measures.
#' So $CMDSESS$MRCSquaredError contains the distribution of RF distances to the expected MRC tree when drawing ESS samples IID from the true posterior (using CMDSESS for calculating ESS).
#' @details There are three options for measuring Monte Carlo error.
#' 1) treeProbSquaredError, yields a vector of squared differences between the per-chain estimate of each topology probability and the average topology probability (averaged over all chains).
#' 2) splitProbSquaredError yields a vector of squared differences between the per-chain estimate of each split probability and the average split probability (averaged over all chains).
#' 3) MRCSquaredError yields a vector of squared RF distances from the MRC of each MCMC chain to the MRC obtained by pooling all chains to compute split frequencies.
#' Both (1) and (2) return matrices, with chains in columns and trees/splits in rows. Summarizing these is up to the user.
#' The trees in treeProbSquaredError are ordered the same as they are in simulated.samples.
#' The splits in splitProbSquaredError are ordered the same as calling as.RFcoords(simulated.samples$trees).
#' The option observed.trees.only is useful for comparing summaries of the Monte Carlo error between ESS iid samples and the MCMC run by guaranteeing no trees are drawn here that are not present in the MCMC chains.
#' @export
effectiveSizeEquivalentError <- function(simulated.samples,tree.dist="RF",measures=c("treeProbSquaredError","splitProbSquaredError","MRCSquaredError"),ess.methods=getESSMethods(),observed.trees.only=TRUE,return.ess=TRUE,verbose=TRUE) {
  # recover()
  
  if ( !("simulatedPosterior" %in% class(simulated.samples) )) {
    stop("'simulated.samples' must be an object of class simulatedPosterior")
  }
  
  ntrees <- length(simulated.samples$trees)
  nchains <- dim(simulated.samples$indices)[2]
  
  # Restrict to sample only observed trees, if desired
  topo_probs <- simulated.samples$probs
  if ( observed.trees.only ) {
    all_visited_trees <- unique(as.integer(simulated.samples$indices))
    seen <- (1:ntrees) %in% all_visited_trees
    unseen <- !seen
    topo_probs[unseen] <- 0
    topo_probs <- topo_probs/sum(topo_probs)
  }
  
  # Start with the distance matrix for all unique topologies (allows us to reduce compute time)
  dmat <- matrix(nrow=ntrees,ncol=ntrees)
  if ( tolower(tree.dist) == "rf" ) {
    dmat <- as.matrix(phangorn::RF.dist(simulated.samples$trees))
  } else if ( tolower(tree.dist) == "spr" ) {
    dmat <- as.matrix(phangorn::SPR.dist(simulated.samples$trees))
  } else {
    stop("Invalid input for argument 'tree.dist'")
  }
  
  # Add the topologies as RF coordinates (allows us to reduce compute time)
  simulated.samples$coords <- trees2Coords(simulated.samples$trees)
  
  # get ESS for every method for every chain
  if ( verbose ) {
    cat("Computing ESS for all measures and all chain\n")
    pb <- txtProgressBar(max=nchains,style=3)
  }
  all_ess <- t(sapply(1:nchains,function(j){
    dm <- expandDistanceMatrix(dmat,simulated.samples$indices[,j])
    phy <- simulated.samples$trees[simulated.samples$indices[,j]]
    these_ess <- sapply(ess.methods,function(ess_method){
      eval(call(ess_method,dmat=dm,trees=phy))
    })
    
    if ( verbose ) {
      setTxtProgressBar(pb,j)
    }
    
    return(these_ess)
  }))
  
  if ( length(ess.methods) == 1 ) {
    all_ess <- t(all_ess)
  }
  
  colnames(all_ess) <- ess.methods
  
  # draw n_eff samples for each chain and each ESS calculated
  if ( verbose ) {
    cat("\nDrawing n_eff iid samples for all measures and all chains\n")
  }
  iid_n_samples_max <- round(max(all_ess))
  iid_samples <- lapply(ess.methods,function(this_method) {
    these_ess <- round(all_ess[,colnames(all_ess) == this_method])
    iid <- simulated.samples
    drawn_indices <- lapply(these_ess,function(n_eff){
      # Some ESS may be in [0,1], but we can't compute any measures without trees
      n_eff <- max(1,n_eff)
      drawn <- sample.int(ntrees,n_eff,replace=TRUE,prob=topo_probs)
      if (n_eff < iid_n_samples_max) {
        drawn <- c(drawn,rep(NA,iid_n_samples_max - n_eff))
      }
      return(drawn)
    })
    iid$indices <- do.call(cbind,drawn_indices)
    return(iid)
  })
  
  # compute the error/variance measures for the iid samples
  if ( verbose ) {
    cat("Computing all measures for all chains\n")
  }
  res <- lapply(iid_samples,function(iid) {
    these_measures <- lapply(measures,function(this_measure) {
      eval(call(this_measure,simulated.samples=iid))
    })
    names(these_measures) <- measures
    return(these_measures)
  })
  names(res) <- ess.methods
  
  if ( return.ess ) {
    for (i in 1:length(ess.methods)) {
      res[[i]]$ESS <- all_ess[,i]
    }
  }
  
  return(res)
}