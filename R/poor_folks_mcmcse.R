#' Stability of MCMC estimates from trees.
#'
#' For some split probability threshold(s), examines stability of the consensus tree from a sample of trees.
#'
#' @param trees The sample of trees from the posterior. A multiPhylo object, or a coordinate matrix as from as.RFcoords(trees)$coords.
#' @param consensus.threshold The consensus tree will contain only splits above this probability (p >= 0.5). Can be a vector of probabilities.
#' @param nrep The number of (bootstrap) replicates used to assess the consensus tree stability.
#' @param tree.probs Trees are taken to be equal in probability (e.g. posterior samples) unless this argument specifies non-uniform probabilities.
#' @param block.size Determines the block size for block bootstrapping. Option "sqrt" scales the block size as the square root of the chain length. Otherwise an integer must be specified. This could be N/ESS.
#' @return A list. $quantiles contains the quantiles specified by probs for all thresholds. $boot contains a matrix of all the bootstrapped distances.
#' @details Stability is assessed as the RF distance from (bootstrap) replicate consensus trees to the consensus tree from all samples. The real sample size taken will not usually be exactly the number of trees because of the block sizes for bootstrapping.
#' @export
treeStability <- function(trees,stat,consensus.threshold=c(0.5,0.75,0.95),nrep=100,probs=c(0.05,0.5,0.95),block.size="sqrt") {
  # recover()
  
  cts <- treeStabilityConvergence(trees,stat=stat,sizes=0,consensus.threshold=consensus.threshold,nrep=nrep,probs=probs,block.size=block.size)[[1]]
  class(cts) <- "list"
  return(cts)
}

#' For some split probability threshold(s), examines stability of the consensus tree from a sample of trees over the length of the run.
#'
#' @param trees The sample of trees from the posterior. A multiPhylo object, or a coordinate matrix as from as.RFcoords(trees)$coords (only usable with stat="ASDSF|MRC")
#' @param stat The summary to be plotted. ASDSF vizualizes convergence of split frequencies, MRC of the summary tree(s), topoProbs of the tree probabilities, imbalance of the tree (im)balance.
#' @param sizes Either an integer giving the number of points to compute the stability (at least 3 and no more than the number of trees), or the sample sizes at which the stability should be computed.
#' @param consensus.threshold The consensus tree will contain only splits above this probability (p >= 0.5). Can be a vector of probabilities.
#' @param nrep The number of (bootstrap) replicates used to assess the consensus tree stability.
#' @param probs The quantiles in the results summary.
#' @param spacing Options "linear|logarithmic". When determining the sequence of sizes (sizes given as an integer), should they be linearly or logarithmically spaced?
#' @param min.split.freq For ASDSF only. Minimum split frequency to be included in ASDSF.
#' @param block.size Determines the block size for block bootstrapping. Option "sqrt" scales the block size as the square root of the chain length. Otherwise an integer must be specified, this is not advised unless only evaluating at the full chain.
#' @return An object of class conTreeStability (a list). For each sample size, $quantiles contains the quantiles specified by probs for all thresholds. $boot contains a matrix of all the bootstrapped distances.
#' @details Stability is assessed as the RF distance from (bootstrap) replicate consensus trees to the consensus tree from all samples. The real sample sizes taken will not be exactly those provided because of the block sizes for bootstrapping.
#' @export
treeStabilityConvergence <- function(trees,stat="MRC",sizes=10,consensus.threshold=c(0.5,0.75,0.95),nrep=100,probs=c(0.05,0.5,0.95),spacing="logarithmic",min.split.freq=0.01,block.size="sqrt",...) {
  #TODO: 
  #      (1) finish implenting all diagnostics
  #      (2) overhaul return to allow running all at once
  
  # recover()
  stat <- toupper(stat)
  
  if ( !(stat %in% c("ASDSF","IMBALANCE","MRC","TOPOPROBS","RFFRECHETVAR")) ) {
    stop("Unrecognized option to argument \"stat\".")
  }
  
  if ( any(consensus.threshold < 0.5) ) {
    stop("Argument \"consensus.threshold\" must be >= 0.5.")
  }
  
  if ( "matrix" %in% class(trees) && stat == "IMBALANCE" ) {
    stop("Argument \"trees\" must be a multiPhylo for stat=\"imbalance\".")
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
  
  all_imbalance <- c()
  if ( stat == "IMBALANCE" ) {
    all_imbalance <- unlist(lapply(trees,unrootedImbalance))
  }
  
  use.dmat <- stat == "TOPOPROBS" || stat == "RFFRECHETVAR"
  tree_indices <- c()
  rf_dmat <- matrix()
  if ( use.dmat ) {
    rf_dmat <- as.matrix(dist(coords,method="manhattan"))
    tree_indices <- rep(NA,ntrees)
    idx <- 0
    while( any(is.na(tree_indices)) ) {
      idx <- idx + 1
      first_unidentified <- min(which(is.na(tree_indices)))
      matches <- which(rf_dmat[first_unidentified,] == 0)
      tree_indices[matches] <- idx
    }
    n_unique_topologies <- max(tree_indices)
  }

  # Figure out target sample sizes for calculating
  sample_sizes <- c()
  if ( length(sizes) == 1 ) {
    # This allows us to call this function once for sample_size == length(trees) in treeStability
    if ( sizes < 2 ) {
      sample_sizes <- ntrees
    } else {
      if ( spacing == "linear" ) {
        sample_sizes <- round(seq(ntrees/sizes,ntrees,length.out=sizes))
      } else {
        sample_sizes <- round(exp(seq(log(nsplits/2),log(ntrees),length.out=sizes)))
      }
    }
  } else {
    if ( any(sizes < 0 ) || any(sizes > length(trees)) ) {
      stop("Cannot evaluate stability at specified sizes.")
    }
    sample_sizes <- sizes
  }
  
  # block-bootstrap batch sizes and information
  # actual sample sizes will be b*a <= sample_sizes
  b <- floor(sqrt(sample_sizes))
  a <- floor(sample_sizes/b)
  
  if ( block.size %% 1 == 0 ) {
    b <- rep(block.size,length(sample_sizes))
    a <- floor(sample_sizes/b)
    if ( min(a) < 2 ) {
      stop("Value of \"block.size\" too large for specified chain subsample sizes.")
    }
  } else if ( block.size != "sqrt" ) {
    stop("Argument \"block.size\" must either be an integer or \"sqrt\".")
  }
  
  # Iterate over (sub)sample sizes
  per_size <- lapply(1:length(sample_sizes),function(k){
    nsamps <- b[k]*a[k]
    
    # Get consensus tree, make trees coordinates for easy bootstrapping
    best_tree_probs <- c()
    best_var <- NA
    best_imbalance <- NA
    best_split_probs <- colMeans(coords[1:nsamps,])
    best_con_splits <- lapply(consensus.threshold,function(thresh) {
      best_split_probs >= thresh
    })
    if ( stat == "TOPOPROBS" ) {
      best_topo_probs <- sapply(1:n_unique_topologies,function(i){
        sum(tree_indices[1:nsamps] == i)/nsamps
      })
    }
    if ( stat == "RFFRECHETVAR") {
      stop("not implemented")
    }
    if ( stat == "IMBALANCE" ) {
      best_imbalance <- mean(all_imbalance[1:nsamps])
    }

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
      
      if ( stat == "ASDSF") {
        dists <- ASDSF(list(split_probs,best_split_probs),min.freq=min.split.freq)
      } else if ( stat == "MRC" ) {
        dists <- sapply(1:length(consensus.threshold),function(j) {
          thresh <- consensus.threshold[j]
          con_splits <- split_probs >= thresh
          sum(con_splits) + sum(best_con_splits[[j]]) - 2*sum(con_splits & best_con_splits[[j]])
        })
      } else if (stat == "TOPOPROBS" ) {
        topo_probs <- sapply(1:n_unique_topologies,function(i){
          sum(tree_indices[idx] == i)/nsamps
        })
        dists <- KL(topo_probs,best_topo_probs)
      } else if ( stat == "IMBALANCE" ) {
        dists <- (mean(all_imbalance[idx]) - best_imbalance)^2
      } else if ( stat == "RFFRECHETVAR" ) {
        
      }
      return(dists)
    })
    boot <- do.call(rbind,boot)
    if ( stat == "MRC" ) {
      colnames(boot) <- paste0(consensus.threshold*100,"%")  
    }
    quants <- apply(boot,2,quantile,probs=probs,type=1)
    return(list(quantiles=quants,boot=boot))
  })
  names(per_size) <- b*a
  class(per_size) <- "treeStability"
  return(per_size)
}

#' Visualizes stability of the consensus tree from a sample of trees over the length of the run.
#'
#' @param x Either a multiPhylo object or the output of treeStabilityConvergence.
#' @param colors Colors for each of the consensus thresholds for plotting. Any provided transparencies will be removed. If NA defaults are used (defaults not guaranteed to be pretty).
#' @param CI.color.alpha Transparency value for plotting CI thresholds
#' @param plot.median Is there a median line to be plotted? Only works if 0.5 is included in treeStabilityConvergence(probs=c(...,0.5,...)).
#' @return Nothing, plots the diagnostic.
#' @details If providing trees, use ... to pass arguments to consensusTreeStabilityConvergence().
#' @export
plotConsensusTreeStabilityConvergence <- function(x,colors=NA,CI.color.alpha=0.5,plot.median=TRUE,...) {
  # recover()
  
  if ( !("treeStability" %in% class(x)) ) {
    x <- consensusTreeStabilityConvergence(x,...)
    
    # Now we strip out the arguments that aren't for plotting
    optional_args <- names(list(...))
    non_plotting <- names(formals("consensusTreeStabilityConvergence"))
    non_plotting <- non_plotting[non_plotting != "..."]
    optional_args <- list(...)[!(optional_args %in% non_plotting)]
  }
  
  # We don't need the raw values, just the quantiles
  x <- lapply(x,function(x_){x_$quantiles})
  
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
  
  sample_sizes <- as.numeric(names(x))
  
  if ( any(is.na(colors)) ) {
    colors <-rainbow(n_thresh)
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
  
  plot(NULL,NULL,xlim=xl,ylim=yl,xlab="# iterations",ylab="RF(boot,real)",...)
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
  
  legend("topright",fill=solid.colors,legend=consensus_threshold,border=NA,bty="n")
  
}
