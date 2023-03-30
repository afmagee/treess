#' Constructs intervals for summaries of trees.
#' 
#' For one or more MCMC chains, constructs confidence/prediction intervals for tree and split probabilities.
#'
#' @param x A list of chains containing the MCMC samples of trees as multiPhylo objects.
#' @param ESS Vector of ESS of each of the MCMC chains in x.
#' @param type Either "confidence" for confidence intervals or "prediction" for prediction intervals.
#' @param interval.width Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95\% confidence/prediction intervals.
#' @param method For confidence intervals, method of interval constriction, "Jefreys"|"Wilson"|"ContinuityCorrectedWilson".
#' @details The function returns a list of lists.
#' The top level of lists is $split and $topology, which separates split and topology probabilities.
#' Within each of these is a list, each element containing the intervals for each chain as a matrix.
#' Each matrix has splits (or trees) in rows (rows are comparable across chains), and the columns are the point estimate, lower, and upper CIs, and the effective sample size.
#' Note that prediction intervals are generated once per chain, and do not account for differences in per-chain ESS.
#' That is, they use the default n.new=n in \link{binomialProportionPI}.
#' @return A list of intervals for each chain, see details.
#' @export
#' @seealso \link{binomialProportionCI}, \link{binomialProportionPI}, \link{treeStability}, \link{plotTreeIntervals}
constructTreeIntervals <- function(x, ESS, type="confidence", interval.width=0.95, method="Jeffreys") {
  # recover()
  
  # Check for valid inputs
  if ( !(any(c("list","mcmc.list") %in% class(x))) ) {
    stop("Argument 'x' must be a list.")
  }
  
  if ( !("multiPhylo" %in% class(x[[1]])) ) {
    stop("Argument \"x\" must be a list of multiPhylo objects.")
  }
  
  nchains <- length(x)
  
  # ensure all chains are of same dimension(s)
  lens <- unlist(lapply(x,length))
  
  if ( length(unique(lens)) != 1 ) {
    stop("All elements of `x` must be of same length/dimension")
  }
  ngen <- lens[1]
  
  # Make trees coordinates
  all_trees <- do.call(c,x)
  coords <- trees2Coords(all_trees)
  
  nsplits <- dim(coords)[2]
  n_unique_topologies <- NA
  
  # Per-chain split probabilities
  split_probs <- lapply(1:nchains,function(i){
    chain_coords <- coords[(ngen*(i-1)+1):(i*ngen),]
    return(colMeans(chain_coords))
  })
  
  # Per-chain tree probabilities
  tree_strings <- apply(coords,1,paste0,collapse="")
  trees <- as.integer(as.factor(tree_strings))
  n_unique_topologies <- max(trees)
  tree_probs <- lapply(1:nchains,function(i){
    chain_trees <- trees[(ngen*(i-1)+1):(i*ngen)]
    chain_probs <- sapply(1:n_unique_topologies,function(k){sum(chain_trees == k)/ngen})
    return(chain_probs)
  })

  interval_fun <- NULL
  if ( grepl("con",type) ) {
    interval_fun <- function(p,n) {
      binomialProportionCI(p,n,method=method,ci.width=interval.width)
    }
  } else if ( grepl("pre",type) ) {
    interval_fun <- function(p,n) {
      binomialProportionPI(p,n,n.new=n,pi.width=interval.width)
    }
  } else {
    stop("Unrecognized \"type\" option.")
  }
  
  chains_splits <- lapply(1:nchains, function(i){
    tmp <- interval_fun(split_probs[[i]],ESS[i])
    tmp <- cbind(split_probs[[i]],tmp,rep(ESS[i],length(split_probs[[i]])))
    colnames(tmp)[1] <- "point.est"
    colnames(tmp)[4] <- "ESS"
    return(tmp)
  })
  names(chains_splits) <- paste0("chain_",1:nchains)
  
  chains_trees <- lapply(1:nchains, function(i){
    tmp <- interval_fun(tree_probs[[i]],ESS[i])
    tmp <- cbind(tree_probs[[i]],tmp,rep(ESS[i],length(tree_probs[[i]])))
    colnames(tmp)[1] <- "point.est"
    colnames(tmp)[4] <- "ESS"
    return(tmp)
  })
  names(chains_trees) <- paste0("chain_",1:nchains)
  
  res <- list(chains_splits,chains_trees)
  names(res) <- c("split","topology")
  
  return(res)
}


#' Plots intervals for summaries of trees.
#' 
#' Plots probabilities of splits or trees estimated by two independent MCMC chains, with confidence/prediction intervals.
#' Differences in probabilities between chains are highlighted.
#' Differences may be obtained by comparing CIs/PIs between runs (default) or by using CIs on the difference in probabilities between runs.
#'
#' @param x Output of \link{constructTreeIntervals}.
#' @param summary "split" to plot split probabilities, "topology" to plot topology probabilities.
#' @param differences Either NA or output of \link{compareChainProbabilities}. See details.
#' @param threshold If given as a numeric vector, dashed lines are drawn at these values as a reference.
#' @param bar.col CIs/PIs are plotted using bar.col[1] if the probabilities are the same between chains and bar.col[2] otherwise. Default is green and red.
#' @param point.col Color of points, default is same as bar color. It is possible to specify a single color.
#' @param bad.only If TRUE, CIs/PIs are only plotted for non-overlapping intervals (bar.col[2]) to highlight chain-chain mismatches.
#' @param log.axes If TRUE, x- and y-axes are logged. Default is FALSE for splits, TRUE for trees.
#' @param xlab The x-axis label, NA for defaults.
#' @param ylab The y-axis label, NA for defaults.
#' @param main The plot title.
#' @details Colors of cross bars show whether or not the probabilities are significantly different, highlighting (dis)agreement between chains.
#' This disagreement can be assessed two ways.
#' The default (differences = NA) is to compare CIs for probabilities, and if they do not overlap declare them significantly different.
#' Alternately, if differences is given the output of \link{compareChainProbabilities}, comparison is based on whether the CI for the difference in probabilities between runs includes 0 or not.
#' Bars that cross thresholds indicate uncertainty with respect to which side of the threshold a split is on.
#' By default, thresholds are plotted at p=0.5 (the boundary for being in an MRC tree), 0.75 (moderate support for a split), and 0.95 (strong support for a split).
#' @return Nothing, generates plot.
#' @export
#' @seealso \link{binomialProportionCI}, \link{binomialProportionPI}, \link{treeStability}, \link{constructTreeIntervals}
plotTreeIntervals <- function(x,summary="split",differences=NA,chains=c(1,2),threshold=c(0.5,0.75,0.95),bar.col=c("springgreen4","firebrick1"),point.col=bar.col,bad.only=FALSE,log.axes=NA,xlab=NA,ylab=NA,main=NA) {
  # recover()
  
  if ( length(chains) !=2 ) {
    stop("Must specify exactly 2 chains.")
  }
  
  to.plot <- grep(summary,names(x))
  if ( length(to.plot) != 1 ) {
    stop("Invalid \"summary\" option.")
  }
  nchains <- length(x[[to.plot]])
  
  do_log <- ifelse(summary == "split","","xy")
  if ( !is.na(log.axes) ) {
    do_log <- ifelse(log.axes,"xy","")
  }
  
  if ( do_log == "xy" ) {
    # Avoid -Inf
    all_p <- unlist(x[[to.plot]])
    min_val <- min(all_p[all_p > 0])
    for (i in 1:length(x[[to.plot]])) {
      x[[to.plot]][[i]][x[[to.plot]][[i]] == 0] <- min_val
    }
  }
  
  if ( is.na(xlab) ) {
    xlab <- paste0("chain ",chains[1])
  }
  
  if ( is.na(ylab) ) {
    ylab <- paste0("chain ",chains[2])
  }
  
  if ( is.na(main) ) {
    main <- ifelse(grepl("split",summary),"split probabilities","topology probabilities")
  }
  
  xlim <- c(0,1)
  ylim <- c(0,1)
  if ( grepl("topo",summary) ) {
    xlim <- range(c(as.numeric(x$topology[[chains[1]]]),as.numeric(x$topology[[chains[2]]])))
    ylim <- xlim
  }
  
  plot(NULL,NULL,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,log=do_log)
  abline(a=0,b=1,col="grey")
  
  # track whether probabilities are significantly different (case 2) or not (case 1)
  cases <- rep(NA,dim(x[[to.plot]][[chains[1]]])[1])
  if ( !is.list(differences) && is.na(differences) ) {
    # track whether intervals overlap (case 1) or not (case 2)
    l1 <- x[[to.plot]][[chains[1]]][,2]
    h1 <- x[[to.plot]][[chains[1]]][,3]
    
    l2 <- x[[to.plot]][[chains[2]]][,2]
    h2 <- x[[to.plot]][[chains[2]]][,3]
    
    cases <- 1 + as.numeric((l1 < l2 & h1 < l2) | (l1 > h2 & h1 > h2))
  } else {
    if ( !("list" %in% class(differences)) || 
         !grepl(summary,names(differences)) || 
         length(differences[[to.plot]]) != choose(nchains,2) ) {
      stop("Unrecognized input to argument \"differences\", must either be NA or output of compareChainProbabilities.")
    }
    # Find comparison of chains[1] against chains[2]
    ij <- grep(paste0("chains_",min(chains),"_vs_",max(chains)),names(differences[[to.plot]]))
    l <- differences[[to.plot]][[ij]][,2]
    h <- differences[[to.plot]][[ij]][,3]
    cases <- 1 + as.numeric((l < 0 & h < 0) | (l > 0 & h > 0))
  }
  
  point_col <- NULL
  if ( length(point.col) == 1 ) {
    point_col <- rep(point.col,length(cases))
  } else {
    point_col <- point.col[cases]
  }
  
  bar_col <- NULL
  if (length(bar.col) == 1) {
    bar_col <- rep(bar.col, length(cases))
  } else {
    bar_col <- bar.col[cases]
  }

  # green for intervals that overlap
  # red otherwise
  for (i in 1:dim(x[[to.plot]][[chains[1]]])[1]) {
    if ( cases[i] == 2 || !bad.only ) {
      xc <- x[[to.plot]][[chains[1]]][i,1]
      yc <- x[[to.plot]][[chains[2]]][i,1]
      
      xr <- x[[to.plot]][[chains[1]]][i,2:3]
      yr <- x[[to.plot]][[chains[2]]][i,2:3]
      
      lines(xr,c(yc,yc),col=bar_col[i])
      lines(c(xc,xc),yr,col=bar_col[i])
    }
  }
  
  points(x[[to.plot]][[chains[1]]][,"point.est"],x[[to.plot]][[chains[2]]][,"point.est"],col=point_col,pch=16)
  
  if ( is.numeric(threshold) ) {
    abline(h=threshold,col="grey",lty=2)
    abline(v=threshold,col="grey",lty=2)
  }
  
  
  
}
