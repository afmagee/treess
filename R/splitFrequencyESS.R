#' Calculates ESS using a Frechet-like generalization of the univariate ESS of Vats and Knudson (2018).
#'
#' Internal matrix-based representation of trees not recommended for large trees.
#'
#' @param trees All trees in the MCMC chain, either as trees or as RF coordinates.
#' @references Vats and Knudson (2021). "Revisiting the Gelman-Rubin diagnostic." Statistical Science, 36(4), 518-529.
#' @keywords internal
splitFrequencyESS <- function(trees,...) {
  # recover()
  
  # As this function is purely internal, we do only minimal checking that this is a coordinate matrix
  coords <- NULL
  n <- 0
  if ( "matrix" %in% class(trees) && all(rowSums(trees) == sum(trees[1,])) && all(colMeans(trees) <= 1) && all(colMeans(trees) >= 0) ) {
    coords <- trees
    n <- dim(trees)[1]
  } else if ( "multiPhylo" %in% class(trees) ) {
    coords <- trees2Coords(trees)
    n <- length(trees)
  } else {
    stop("Cannot compute splitFrequencyESS for input other than multiPhylo or RF coordinate matrix.")
  }
  
  # Global mean
  mu <- colMeans(coords)
  
  # Global Frechet variance
  sig_sq <- sum(sapply(1:n,function(i){
    dist(rbind(mu,coords[i,]))^2
  }))/(n-1)
  
  # b = batch size, a = # batches
  b <- floor(n^(1/2))
  a <- floor(n/b)
  
  ybar <- t(sapply(1:a,function(k){
    first <- (k-1)*b+1
    last <- k*b
    colMeans(coords[first:last,])
  }))
  
  tau_sq_b <- b/(a-1) * sum(sapply(1:a,function(k){dist(rbind(ybar[k,],mu))^2}))
  
  # b = batch size, a = # batches
  b <- floor((n^(1/2))/3)
  a <- floor(n/b)
  
  ybar <- t(sapply(1:a,function(k){
    first <- (k-1)*b+1
    last <- k*b
    colMeans(coords[first:last,])
  }))
  
  tau_sq_b_3 <- b/(a-1) * sum(sapply(1:a,function(k){dist(rbind(ybar[k,],mu))^2}))
  
  tau_sq <- 2*tau_sq_b - tau_sq_b_3
  
  cat(tau_sq,tau_sq_b,tau_sq_b_3,sig_sq,"\n",sep=";")
  
  return(n*sig_sq/tau_sq)
  
}

#' Calculates ESS using a Frechet-like generalization of the univariate ESS of Vats and Knudson (2018).
#'
#' Slightly slower than matrix-coordinates version for small trees, much faster for large trees.
#'
#' @param trees All trees in the MCMC chain as a multiPhylo object.
#' @param rooted Should we be using clades (for rooted trees only) or splits
#' @references Vats, Dootika, and Christina Knudson. "Revisiting the Gelman-Rubin diagnostic." arXiv preprint arXiv:1812.09384 (2018).
#' @keywords internal
splitFrequencyESS2 <- function(trees,rooted,max.b.3=NA,...) {
  # recover()
  
  # As this function is purely internal, we do only minimal checking that this is a coordinate matrix
  if ( !("multiPhylo" %in% class(trees)) ) {
    stop("Cannot compute splitFrequencyESS for input other than multiPhylo.")
  }
  
  n <- length(trees)
  
  # b = batch size, a = # batches
  b <- floor(n^(1/2))
  a <- floor(n/b)
  
  # use <- 1:b
  # if ( !is.na(max.b.3) && max.b.3 < floor(b/3)) {
  #   max_b_3 <- round(max.b.3)
  #   thin_factor <- floor(floor(b/3)/max_b_3)
  #   b <- 3 * max_b_3
  #   trees <- trees[seq(1,n,thin_factor)]
  #   a <- floor(length(trees)/b)
  # }
  
  splits <- treeSplits(trees,rooted=rooted)

  # Global mean
  mu <- splitFrechetMean(splits)

  # Global Frechet variance
  sig_sq <- splitFrechetVariance(splits,mean=mu,bessel=TRUE)

  ybar <- lapply(1:a,function(k){
    first <- (k-1)*b+1
    last <- k*b
    splitFrechetMean(splits[first:last])
  })

  summand <- sapply(1:length(ybar),function(i){
    key <- fastmatch::fmatch(attr(ybar[[i]],"labels"),attr(mu,"labels"))
    unseen_sum_sq <- sum(mu[-key]^2)
    seen_sum_sq <- sum((ybar[[i]] - mu[key])^2)
    return(seen_sum_sq + unseen_sum_sq)
  })

  tau_sq_b <- b/(a-1) * sum(summand)
  
  # b = batch size, a = # batches
  b <- floor(b/3)
  a <- floor(n/b)
  
  ybar <- lapply(1:a,function(k){
    first <- (k-1)*b+1
    last <- k*b
    splitFrechetMean(splits[first:last])
  })
  
  summand <- sapply(1:length(ybar),function(i){
    key <- fastmatch::fmatch(attr(ybar[[i]],"labels"),attr(mu,"labels"))
    unseen_sum_sq <- sum(mu[-key]^2)
    seen_sum_sq <- sum((ybar[[i]] - mu[key])^2)
    return(seen_sum_sq + unseen_sum_sq)
  })
  
  tau_sq_b_3 <- b/(a-1) * sum(summand)

  tau_sq <- 2*tau_sq_b - tau_sq_b_3
  
  cat(tau_sq,tau_sq_b,tau_sq_b_3,sig_sq,"\n",sep=";")
  
  return(n*sig_sq/tau_sq)
  
}


#' #' Calculates ESS using a Frechet-like generalization of the univariate ESS of Vats and Knudson (2018).
#' #'
#' #' Slightly slower than matrix-coordinates version for small trees, much faster for large trees.
#' #'
#' #' @param trees All trees in the MCMC chain as a multiPhylo object.
#' #' @param useClades Should we be using clades (for rooted trees only) or splits
#' #' @references Vats, Dootika, and Christina Knudson. "Revisiting the Gelman-Rubin diagnostic." arXiv preprint arXiv:1812.09384 (2018).
#' #' @keywords internal
#' splitFrequencyESS <- function(trees,rooted,max.b.3=NA,...) {
#'   # recover()
#'   
#'   # As this function is purely internal, we do only minimal checking that this is a coordinate matrix
#'   if ( !("multiPhylo" %in% class(trees)) ) {
#'     stop("Cannot compute splitFrequencyESS for input other than multiPhylo.")
#'   }
#'   
#'   n <- length(trees)
#'   
#'   # b = batch size, a = # batches
#'   b <- floor(n^(1/2))
#'   a <- floor(n/b)
#'   
#'   use <- 1:b
#'   if ( !is.na(max.b.3) && max.b.3 < floor(b/3)) {
#'     max_b_3 <- round(max.b.3)
#'     thin_factor <- floor(floor(b/3)/max_b_3)
#'     b <- 3 * max_b_3
#'     trees <- trees[seq(1,n,thin_factor)]
#'     a <- floor(length(trees)/b)
#'   }
#'   
#'   splits <- perTreeSplits(trees,useClades=rooted,return.names=FALSE)
#'   
#'   # Global mean
#'   mu <- splitFrechetMean(splits)
#'   
#'   # Global Frechet variance
#'   sig_sq <- splitFrechetVariance(splits,mean=mu,bessel=TRUE)
#'   
#'   ybar <- lapply(1:a,function(k){
#'     first <- (k-1)*b+1
#'     last <- k*b
#'     splitFrechetMean(splits[first:last])
#'   })
#'   
#'   summand <- sapply(1:length(ybar),function(i){
#'     seen <- names(mu) %in% names(ybar[[i]])
#'     unseen_sum_sq <- sum(mu[!seen]^2)
#'     key <- match(names(ybar[[i]]),names(mu))
#'     seen_sum_sq <- sum((ybar[[i]] - mu[key])^2)
#'     return(seen_sum_sq + unseen_sum_sq)
#'   })
#'   
#'   tau_sq_b <- b/(a-1) * sum(summand)
#'   
#'   # b = batch size, a = # batches
#'   b <- floor(b/3)
#'   a <- floor(length(trees)/b)
#'   
#'   ybar <- lapply(1:a,function(k){
#'     first <- (k-1)*b+1
#'     last <- k*b
#'     splitFrechetMean(splits[first:last])
#'   })
#'   
#'   summand <- sapply(1:length(ybar),function(i){
#'     seen <- names(mu) %in% names(ybar[[i]])
#'     unseen_sum_sq <- sum(mu[!seen]^2)
#'     key <- match(names(ybar[[i]]),names(mu))
#'     seen_sum_sq <- sum((ybar[[i]] - mu[key])^2)
#'     return(seen_sum_sq + unseen_sum_sq)
#'   })
#'   
#'   tau_sq_b_3 <- b/(a-1) * sum(summand)
#'   
#'   tau_sq <- 2*tau_sq_b - tau_sq_b_3
#'   
#'   return(n*sig_sq/tau_sq)
#'   
#' }
