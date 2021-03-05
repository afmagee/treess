#' Examines stability of tree-based summaries from a sample of trees using a block-bootstrap to resample the MCMC run.
#'
#' This is a convenience wrapper for treeStabilityConvergence.
#' It allows you to compute uncertainty about tree summaries for the full chain and only the full chain.
#' For more, see \link{treeStabilityConvergence}.
#'
#' @param trees The sample of trees from the posterior. A multiPhylo object, or a coordinate matrix as from as.RFcoords(trees)$coords.
#' @param nrep The number of bootstrap replicates used to assess the tree stability.
#' @param probs The quantiles in the results summary.
#' @param consensus.threshold The consensus tree will contain only splits above this probability (p > 0.5). Can be a vector of probabilities.
#' @param ACT Autocorrelation time for computing the block size for bootstrapping the chain. Default ("sqrt") conservatively scales ACT with sqrt of subsample length. Can also be provided as numeric.
#' @param min.split.freq For ASDSF only. Minimum split frequency to be included in ASDSF.
#' @return One element of the list returned by \link{treeStabilityConvergence}.
#' @seealso \link{treeStabilityConvergence},\link{plotTreeStabilityConvergence}
#' @export
treeStability <- function(trees,consensus.threshold=c(0.5,0.75,0.95),nrep=100,min.split.freq=0.01,probs=c(0.05,0.5,0.95),ACT="sqrt") {
  # recover()
  
  cts <- treeStabilityConvergence(trees,consensus.threshold=consensus.threshold,nrep=nrep,probs=probs,ACT=ACT,full.chain.only=TRUE)[[1]]
  return(cts)
}

#' Examines stability of tree-based summaries from a sample of trees using a block-bootstrap to resample the MCMC run.
#'
#' By performing this resampling at a variety of chain lengths, provides information about how estimates have converged over the run.
#' Summaries include ASDSF, tree topology probabilities, and MRC-type trees.
#' This function can be slow with large samples of trees or when computing values at many subsample sizes.
#'
#' @param trees The sample of trees from the posterior. A multiPhylo object, or a coordinate matrix as from as.RFcoords(trees)$coords.
#' @param nsizes The number of points to compute the stability, the number of chain subsample lengths.
#' @param nrep The number of bootstrap replicates used to assess the tree stability at each chain subsample length.
#' @param probs The quantiles in the results summary.
#' @param consensus.threshold The consensus tree will contain only splits above this probability (p > 0.5). Can be a vector of probabilities.
#' @param spacing Options "linear|logarithmic". When determining the chain subsample sizes, should they be spaced linearly or logarithmically?
#' @param ACT Autocorrelation time for computing the block size for bootstrapping the chain. Default ("sqrt") conservatively scales ACT with sqrt of subsample length. Can also be provided as numeric.
#' @param min.split.freq For ASDSF only. Minimum split frequency to be included in ASDSF.
#' @param treess.method If treess.method is an option in getESSMethods(), for each block size we compute its treess. Default ("none") is to not compute treess to save time. Only one measure may be supplied.
#' @return An object of class conTreeStability (a nested list). See details.
#' @details 
#' Stability is assessed by comparing summaries computed from bootstrapped samples of the Markov chain to the summary from the full sample.
#' For each subsample length n, the summary is computed for trees[1:n].
#' Then, nrep new chains of length n are generated via a block bootstrap, and the distance between each bootstrapped chain and the original subsampled chain are computed.
#' 
#' The summaries are the (1) vector of split frequencies, (2) the vector of tree topology probabilities, and (3) the consensus trees at given thresholds.
#' Split frequencies are compared using the ASDSF (equivalent to 1/nsplits times the Manhattan distance between the two vectors).
#' Topology probabilities are compared using the Euclidean distance between the vectors.
#' Consensus trees are compared using the RF distance.
#' 
#' The real sample sizes taken will not be exactly those provided because of the block sizes for bootstrapping.
#' 
#' When computing treess for the block sizes, the RF distance is used, rather than a user-specified option, to save on computation time.
#' 
#' The output is a nested list.
#' $using_ess records whether the ESS was calculated.
#' Each other element in the list contains summaries for a given subsample length (named by the length or ESS if treess.method != "none").
#' Each of these lists is a list of two matrices.
#' $quantiles summarizes the chosen quantiles of the bootstrapped distances.
#' $boot contains the nrep bootstrapped distance measures (replicates in rows, statistics in columns).
#' $ess contains the computed treess value, if one is in fact computed.
#' $chain.length records the length of the subsampled chain (which contains trees[1:chain.length]), allowing other ess measures to be computed.
#' Both matrices have labels for the statistics.
#' @seealso \link{treeStability},\link{plotTreeStabilityConvergence}
#' @export
treeStabilityConvergence <- function(trees,
                                     nsizes=10,
                                     nrep=100,
                                     probs=c(0.05,0.5,0.95),
                                     consensus.threshold=c(0.5,0.75,0.95),
                                     ACT="sqrt",
                                     spacing="logarithmic",
                                     min.split.freq=0.01,
                                     treess.method="none",
                                     ...) {

  # recover()

  if ( any(consensus.threshold < 0.5) ) {
    stop("Argument \"consensus.threshold\" must be >= 0.5.")
  }
  
  # Make trees coordinates
  coords <- matrix()
  if ( "multiPhylo" %in% class(trees) ) {
    coords <- trees2Coords(trees)
  } else if ( "matrix" %in% class(trees) ) {
    coords <- trees
  } else {
    stop("Argument \"trees\" must either be a multiPhylo or an RF coordinate matrix.")
  }
  
  ntrees <- dim(coords)[1]
  nsplits <- dim(coords)[2]
  n_unique_topologies <- NA
  
  tree_indices <- apply(coords,1,paste0,collapse="")
  tree_indices <- as.integer(as.factor(tree_indices))
  n_unique_topologies <- max(tree_indices)
  
  # We need to know about the used-defined ACT to generate valid sample sizes, and we need integer-valued sample sizes
  if ( is.numeric(ACT) ) {
    ACT <- ceiling(ACT)
  }
  
  # Figure out target sample sizes for calculating
  sample_sizes <- c()
  if ( hasArg("full.chain.only") && list(...)$full.chain.only == TRUE ) {
    # This allows us to call this function once for sample_size == length(trees) in treeStability
    sample_sizes <- ntrees
  } else {
    # We need chains long enough to
    #   1) capture at least some of the splits observed
    #   2) compute ESS measures
    #   3) get at least a few batches out of given the ACT
    # The choice of the lower bound that meets these criteria may change in future versions
    min_chain_length <- max(nsplits/2,50)
    if ( is.numeric(ACT) ) {
      min_chain_length <- max(2*ACT,min_chain_length)
    }
    if ( spacing == "linear" ) {
      # Uniformly cover chain lengths linearly between min,max
      sample_sizes <- round(seq(min_chain_length,ntrees,length.out=nsizes))
    } else {
      # Log-uniform sampling puts more mass at smaller chain lengths where more variability is expected
      sample_sizes <- round(exp(seq(log(min_chain_length),log(ntrees),length.out=nsizes)))
    }
  }
  
  # block-bootstrap batch sizes and information
  # actual sample sizes will be b*a <= sample_sizes
  b <- floor(sqrt(sample_sizes))
  a <- floor(sample_sizes/b)
  
  if ( is.numeric(ACT) ) {
    b <- rep(ACT,length(sample_sizes))
    a <- floor(sample_sizes/b)
    if ( min(a) < 2 ) {
      stop("Value of \"ACT\" too large for specified chain subsample sizes.")
    }
  } else if ( ACT != "sqrt" ) {
    stop("Argument \"ACT\" must either be an integer or \"sqrt\".")
  }
  
  # Iterate over (sub)sample sizes
  per_size <- lapply(1:length(sample_sizes),function(k){
    nsamps <- b[k]*a[k]
    
    # Get consensus tree, make trees coordinates for easy bootstrapping
    best_tree_probs <- c()
    best_split_probs <- colMeans(coords[1:nsamps,])
    best_con_splits <- lapply(consensus.threshold,function(thresh) {
      best_split_probs > thresh
    })
    best_topo_probs <- sapply(1:n_unique_topologies,function(i){
      sum(tree_indices[1:nsamps] == i)/nsamps
    })

    # Run bootstrapping
    boot <- lapply(1:nrep,function(n){
      # Get chain indices of bootstrapped replicates
      starts <- sample.int(nsamps - b[k] + 1,a[k],replace=TRUE)
      idx <- unlist(lapply(starts,function(s){s:(s+b[k]-1)}))
      
      # Tree probs (need for all)
      tree_probs <- sapply(1:ntrees,function(i){sum(idx == i)/length(idx)})
      
      # For returning
      dists <- c()
      
      split_probs <- colSums(tree_probs * coords)
      
      asdsf <- ASDSF(list(split_probs,best_split_probs),min.freq=min.split.freq)
      mrc <- sapply(1:length(consensus.threshold),function(j) {
        thresh <- consensus.threshold[j]
        con_splits <- split_probs > thresh
        sum(con_splits) + sum(best_con_splits[[j]]) - 2*sum(con_splits & best_con_splits[[j]])
      })
      tp <- sapply(1:n_unique_topologies,function(i){
        sum(tree_indices[idx] == i)/nsamps
      })
      topo_probs <- as.numeric(dist(rbind(tp,best_topo_probs)))
      dists <- c(asdsf,topo_probs,mrc)
      return(dists)
    })
    boot <- do.call(rbind,boot)
    colnames(boot) <- c("ASDSF","TOPOPROBS",paste0("MRC_",consensus.threshold*100,"%")  )
    quants <- t(apply(boot,2,quantile,probs=probs,type=1))
    return(list(chain.length=nsamps,quantiles=quants,boot=boot))
  })
  if ( treess.method %in% getESSMethods() ) {
    # Check for user input of these parameters
    alpha <- 0.05
    if ( hasArg(alpha) ) {
      alpha <- list(...)$alpha
    }
    min.nsamples <- 5
    if ( hasArg(alpha) ) {
      min.nsamples <- list(...)$min.nsamples
    }
    nsim <- 1000
    if ( hasArg(nsim) ) {
      nsim <- list(...)$nsim
    }
    
    # Avoid computing distance matrix if possible
    dmat <- matrix(NA,ntrees,ntrees)
    if ( treess.method != "splitFrequencyESS" ) {
      # RF distance is manhattan distance between trees as points in RF coordinate space
      dmat <- as.matrix(dist(coords),method="manhattan")
    }
    for (i in 1:length(sample_sizes)) {
      ns <- (a*b)[i]
      ess <- eval(call(treess.method,dmat=dmat[1:ns,1:ns],trees=coords[1:ns,],min.nsamples=min.nsamples,alpha=alpha,nsim=nsim))
      per_size[[i]]$ess <- ess
    }
  }
  class(per_size) <- "treeStability"
  return(per_size)
}

#' Visualizes stability of the tree estimates from a sample of trees over the length of the run.
#'
#' @param tree.stability.convergence Output of treeStabilityConvergence.
#' @param stat Which summary of tree stability to ploe? ASDSF|MRC|TOPOPROBS
#' @param colors Colors for each of the consensus thresholds for plotting. Any provided transparencies will be removed. If NA defaults are used (defaults not guaranteed to be pretty).
#' @param CI.color.alpha Transparency value for plotting CI thresholds
#' @param plot.median Is there a median line to be plotted? Only works if 0.5 is included in treeStabilityConvergence(probs=c(...,0.5,...)).
#' @param use.ess If TRUE, and there are computed ESS values in tree.stability.convergence, plots ESS for x-axis. Otherwise x-axis is raw number of samples.
#' @param xlab The x-axis label, NA for defaults.
#' @param ylab The y-axis label, NA for defaults.
#' @param main The plot title.
#' @param ... Further arguments to be passed to plotting functions.
#' @return Nothing, plots the diagnostic.
#' @details 
#' Note that the option use.ess is purely cosmetic, as ESS is only computed for the unbootstrapped subchain.
#' If choosing to plot the x-axis based on ESS, and the ESS is not estimated to be monotonically increasing in subchain length, sensible plotting results cannot be guaranteed.
#' @seealso \link{treeStability},\link{treeStabilityConvergence}
#' @export
plotTreeStabilityConvergence <- function(tree.stability.convergence,stat,colors=NA,CI.color.alpha=0.5,plot.median=TRUE,use.ess=TRUE,xlab=NA,ylab=NA,main=NA,...) {
  # recover()
  stat <- toupper(stat)
  if ( !(stat %in% c("ASDSF","MRC","TOPOPROBS")) ) {
    stop("Unrecognized option to argument \"stat\".")
  }
  
  if ( !("treeStability" %in% class(tree.stability.convergence)) ) {
    stop("Invalid input to argument tree.stability.convergence")
  }
  using_ess <- use.ess && !is.null(tree.stability.convergence[[1]]$ess)
  
  if ( is.na(xlab) ) {
    if ( using_ess) {
      xlab <- "ESS"
    } else {
      xlab <- "# samples"
    }
  }
  
  if ( is.na(ylab) ) {
    if ( stat == "MRC" ) {
      ylab <- "RF distance"
    } else if ( stat == "TOPOPROBS" ) {
      ylab <- "Euclidean distance"
    } else if ( stat == "ASDSF" ) {
      ylab <- "ASDSF"
    }
  }
  
  summary_is <- ""
  if ( stat == "MRC" ) {
    summary_is <- "consensus tree"
  } else if ( stat == "TOPOPROBS" ) {
    summary_is <- "tree probabilities"
  } else if ( stat == "ASDSF" ) {
    summary_is <- "split probabilities"
  }
  
  
  if ( is.na(main) ) {
    main <- paste0("distance between ",summary_is,",\nbootstrapped MCMC chain to real MCMC chain")
  }
  
  # Strip out informational elements
  x <- tree.stability.convergence
  
  # Simplify output, keep only quantiles for the statistic of interest
  # Transposing bridges slightly more human readable output with format assumed below
  x <- lapply(x,function(x_){
    tmp <- x_$quantiles
    tmp <- tmp[grepl(stat,rownames(tmp)),]
    if ( stat == "MRC" ) {
      rownames(tmp) <- gsub("MRC_","",rownames(tmp))
      return(t(tmp))
    } else {
      tmp <- cbind(tmp)
      colnames(tmp) <- NULL
      return(tmp)
    }
  })
  
  # information about the quantiles used in CI calculations
  consensus_threshold <- colnames(x[[1]])
  n_thresh <- dim(x[[1]])[2]
  
  CI_threshold <- rownames(x[[1]])
  
  if ( plot.median ) {
    CI_threshold <- CI_threshold[CI_threshold != "50%"]
  }
  
  if ( (length(CI_threshold) %% 2) != 0 ) {
    stop("There are an odd number of CIs to plot.")
  }
  n_ci <- length(CI_threshold)/2
  
  sample_sizes <- numeric(length(x))
  if ( using_ess ) {
    sample_sizes <- unlist(lapply(tree.stability.convergence,function(tsc){tsc$ess}))
  } else {
    sample_sizes <- unlist(lapply(tree.stability.convergence,function(tsc){tsc$chain.length}))
  }
  
  if ( any(is.na(colors)) ) {
    colors <-heat.colors(n_thresh)
  } else {
    if ( !(all(lengths(strsplit(colors,"")) %in% c(7,9)) && all(grepl("#",colors,fixed=TRUE))) ) {
      stop("Argument \"colors\" must be entirely hexadecimal color values")
    }
    if ( any(lengths(strsplit(colors,"")) == 9) ) {
      colors <- lapply(colors,function(co){
        str <- strsplit(co,"")[[1]]
        if ( length(str) == 9) {
          str <- str[1:7]
        }
        return(paste0(str,collapse=""))
      })
      colors <- unlist(colors)
    }
  }
  solid.colors <- colors
  ci.colors <- adjustcolor(colors,alpha.f=CI.color.alpha)
  
  # plot range
  xl <- range(sample_sizes)
  yl <- range(unlist(x))
  if ( hasArg("ylim") ) {
    yl <- ylim
  }
  
  plot(NULL,NULL,xlim=xl,ylim=yl,xlab=xlab,ylab=ylab,main=main,...)
  xc <- c(sample_sizes,rev(sample_sizes))
  for (i in 1:n_thresh) {
    for (j in 1:n_ci) {
      ci <- sapply(1:length(sample_sizes),function(n){
        all <- x[[n]][,i]
        if ( plot.median ) {
          all <- all[names(all) != "50%"]
        }
        all[c(j,length(all)-j+1)]
      })
      yc <- c(ci[1,],rev(ci[2,]))
      polygon(xc,yc,col=ci.colors[i],border=NA)
    }
    if ( plot.median ) {
      y <- sapply(1:length(sample_sizes),function(n){
        all <- x[[n]][,i]
          all[names(all) == "50%"]
      })
      lines(sample_sizes,y,col=solid.colors[i],...)
    }
  }
  
  if ( !is.null(consensus_threshold) ) {
    legend("topright",fill=solid.colors,legend=consensus_threshold,border=NA,bty="n")
  }

}
