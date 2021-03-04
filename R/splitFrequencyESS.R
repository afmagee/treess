#' Calculates ESS using a Frechet-like generalization of the univariate ESS of Vats and Knudson (2018).
#'
#' @param dmat For compatibility with eval and call, not used here.
#' @param trees All trees in the MCMC chain, either as trees or as RF coordinates.
#' @param nsim For compatibility with eval and call, not used here.
#' @param alpha For compatibility with eval and call, not used here.
#' @param min.nsamples For compatibility with eval and call, not used here.
#' @references Vats, Dootika, and Christina Knudson. "Revisiting the Gelman-Rubin diagnostic." arXiv preprint arXiv:1812.09384 (2018).
#' @keywords internal
splitFrequencyESS <- function(dmat=NA,trees,min.nsamples=NA,alpha=NA,nsim=NA) {
  # recover()
  
  # As this function is putely internal, we do only minimal checking that this is a coordinate matrix
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
  
  tau_sq_b <- b/(a-1) * sum(sapply(1:a,function(k){dist(rbind(ybar[k,],mu))}))
  
  # b = batch size, a = # batches
  b <- floor((n^(1/2))/3)
  a <- floor(n/b)
  
  ybar <- t(sapply(1:a,function(k){
    first <- (k-1)*b+1
    last <- k*b
    colMeans(coords[first:last,])
  }))
  
  tau_sq_b_3 <- b/(a-1) * sum(sapply(1:a,function(k){dist(rbind(ybar[k,],mu))}))
  
  tau_sq <- 2*tau_sq_b - tau_sq_b_3
  
  return(n*sig_sq/tau_sq)
  
}

