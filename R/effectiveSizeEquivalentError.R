#' Estimates distribution of error if we drew ESS trees IID from the (known) posterior distribution.
#'
#' @param simulated.samples An object of class simulatedPosterior (output of simulatePhylogeneticMCMC).
#' @param tree.dist The distance measure for trees (RF or SPR, only one).
#' @param measures The error or variance measure(s) (see details).
#' @param ess.methods The ESS calculation method(s) (see details).
#' @param return.ess Should the returned lists include the calculated ESS for each chain? 
#' @param verbose Should progress be printed to screen?
#' @return The first layer are the different ess methods used, the second the performance measures.
#' So $CMDSESS$distanceToTrueMRC contains the distribution of RF distances to the true MRC tree when drawing ESS samples IID from the true posterior (using CMDSESS for calculating ESS).
#' @export
effectiveSizeEquivalentError <- function(simulated.samples,tree.dist=c("RF","SPR"),measures=c("treeProbSquaredError","treeVarianceSquaredError","splitProbSquaredError","distanceToTrueMRC"),ess.methods=getESSMethods(),return.ess=TRUE,verbose=TRUE) {
  # recover()
  
  if ( !("simulatedPosterior" %in% class(simulated.samples) )) {
    stop("'simulated.samples' must be an object of class simulatedPosterior")
  }
  
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
      drawn <- sample.int(ntrees,n_eff,replace=TRUE,prob=iid$probs)
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
      eval(call(this_measure,simulated.samples=iid,dmat=dmat))
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