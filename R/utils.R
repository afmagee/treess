#' Prediction intervals for proportions.
#' 
#' Given a sample size and an estimated proportion, returns interval in which proportions computed from new samples may lie.
#' Useful for turning estimates of the treESS into ranges of plasible split and topology probabilities from independent MCMC runs.
#'
#' @param p The (estimated) probability.
#' @param n The sample size (or effective sample size) from which the probabilities were computed.
#' @param n.new The sample size (or effective sample size) of the new dataset for which prediction intervals are to be computed.
#' @param pi.width Width of the desired prediction interval (value in [0,1]). Equivalent to 1 - 2*alpha.
#' @return A length(p)x2 matrix containing the lower and upper ends of the prediction interval.
#' @export
#' @details Method uses a Bayesian approach and is thus most similar to method="Jeffreys" in \link{binomialProportionCI}.
#' The method returns probabilities from a new sample, counts divided by the sample size.
#' @seealso \link{binomialProportionCI}
#' @examples
#' 
#' p <- rbeta(10,1,1) # Generate probabilities
#' n <- 100 # Sample size
#' x <- rbinom(length(p),n,p) # Random variables
#' 
#' # compute CIs
#' binomialProportionPI(x/n,n)

binomialProportionPI <- function(p,n,n.new=n,pi.width=0.95) {
  if (!requireNamespace("rmutil",quietly=TRUE)) {
    stop("Cannot compute binomial prediction intervals without package rmutil Please install rmutil to continue.")
  }
  
  # recover()
  
  # PI bounds
  lb <- ub <- numeric(length(p))
  
  # Ends of PIs in probability
  a_2 <- (1-pi.width)/2
  one_minus_a_2 <- 1 - a_2
  
  # With Beta(0.5,0.5) prior on p and binomial likelihood, posterior is Beta with alpha=a and beta=b
  a <- p*n + 0.5
  b <- (1-p)*n + 0.5
  
  # rmutil parameterizes m = a/(a+b) and s = a+b
  m <- a/(a+b)
  s <- a+b
  
  lb <- rmutil::qbetabinom(a_2,round(n.new),m,s)
  ub <- rmutil::qbetabinom(one_minus_a_2,round(n.new),m,s)
    
  pi <- cbind(lb,ub)
  colnames(pi) <- paste0(round(100*c(a_2,one_minus_a_2),1),"%")
  return(pi/n.new)
}

#' Confidence intervals for proportions.
#' 
#' Provides three methods for computing confidence intervals for binomial proportions.
#' Useful for turning estimates of the treESS into CIs for split and topology probabilities.
#'
#' @param p The (estimated) probability.
#' @param n The sample size (or effective sample size).
#' @param method Jefreys|Wilson|ContinuityCorrectedWilson.
#' @param ci.width Width of the desired confidence interval (value in [0,1]). Equivalent to 1 - 2*alpha.
#' @return A length(p)x2 matrix containing the lower and upper ends of the CI.
#' @export
#' @details In practice, the Jeffreys and Wilson intervals often seem to give similar results, while the Corrected Wilson interval often gives much wider intervals.
#' @references Newcombe (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods." Statistics in Medicine.
#' @seealso \link{binomialProportionPI}
#' @examples
#' 
#' p <- rbeta(100,1,1) # Generate probabilities
#' n <- 100 # Sample size
#' x <- rbinom(length(p),n,p) # Random variables
#' 
#' # compute CIs
#' jci <- binomialProportionCI(x/n,n,"jeffreys")
#' wci <- binomialProportionCI(x/n,n,"wilson")
#' cci <- binomialProportionCI(x/n,n,"wilson.corrected")
#' 
#' # compare performance
#' sum(x/n >= jci[,1] & x/n <= jci[,2])/length(p)
#' sum(x/n >= wci[,1] & x/n <= wci[,2])/length(p)
#' sum(x/n >= cci[,1] & x/n <= cci[,2])/length(p)

binomialProportionCI <- function(p,n,method,ci.width=0.95) {
  # recover()
  # CI bounds
  lb <- ub <- numeric(length(p))
  # Ends of CIs in probability
  a_2 <- (1-ci.width)/2
  one_minus_a_2 <- 1 - a_2
  
  # Bounds are slightly different at p=0,p=1 for Jeffreys and Corrected Wilson intervals
  is_0 <- p < .Machine$double.eps
  is_1 <- abs(1 - p) < .Machine$double.eps
  
  if ( (grepl("wil",tolower(method))) ) {
    z <- qnorm(one_minus_a_2,0,1)
    two_n_z_sq <- (2*(n+z^2))
    if ( grepl("cor",tolower(method)) ) {
      # Newcombe (1998) Eqn 4
      # Some negative square roots may ensue, but these will be corrected in the next step, no need to alarm the user
      suppressWarnings({
      lb <- (2*n*p + z^2 - 1 - z*sqrt(z^2 - 2 - 1/n + 4*p*(n*(1-p) + 1))) / two_n_z_sq
      ub <- (2*n*p + z^2 + 1 + z*sqrt(z^2 + 2 - 1/n + 4*p*(n*(1-p) + 1))) / two_n_z_sq
      })
      
      # boundaries, fixes invalid answers
      lb[is_0] <- 0
      ub[is_1] <- 1
    } else {
      # Newcombe (1998) Eqn 3
      lb <- (2*n*p + z^2 - z*sqrt(z^2 + 4*n*p*(1-p))) / two_n_z_sq
      ub <- (2*n*p + z^2 + z*sqrt(z^2 + 4*n*p*(1-p))) / two_n_z_sq
    }
  } else if ( (grepl("jef",tolower(method))) ) {
    # Jeffreys confidence interval, Bayesian credible interval with Beta(0.5,0.5) prior on p
    a <- p*n + 0.5
    b <- (1-p)*n + 0.5
    
    lb <- qbeta(a_2,a,b)
    ub <- qbeta(one_minus_a_2,a,b)
    
    # boundaries
    lb[is_0] <- 0
    ub[is_0] <- (1 - (a/2)^(1/n))[is_0]
    
    lb[is_1] <- ((a/2)^(1/n))[is_1]
    ub[is_1] <- 1
  } else {
    stop("Invalid \"method\"")
  }
  ci <- cbind(lb,ub)
  colnames(ci) <- paste0(round(100*c(a_2,one_minus_a_2),1),"%")
  return(ci)
}

#' Finds s0 for a piecewise constant curve by linear interpolation.
#' 
#' Curve is assumed to be nondecreasing.
#'
#' @param x The time lags at which measurements are taken (MUST INCLUDE 0!)
#' @param y The measurement associate with each x that we are using for to compute the ESS.
#' @param y.crit The critical value above which independence is achieved.
#' @keywords internal
finds0Smoothed <- function(x,y,y.crit) {
  n <- length(y)
  
  if ( (!x[1] == 0 && length(x) == n) ) {
    stop("findt0Smoothed requires x[1] = 0 and length(x) = length(y)")
  }
  if ( any(y[-1] < y[-n]) ) {
    stop("findt0Smoothed requires nondecreasing y")
  }
  
  # Find change points and remove replicate entries, track x coordinates
  # x=(1,2,3,4,5,6),y=(1,1,1,2,2,3) -> x=(1,4,6),y=(1,2,3)
  if ( (any(y[-1] == y[-n])) ) {
    first <- match(unique(y),y)
    x <- x[first]
    y <- y[first]
  }
  
  first_larger <- min(which(y > y.crit))
  
  if (first_larger == 0) {
    stop("t0 cannot be negative")
  }
  
  x2 <- x[first_larger]
  x1 <- x[first_larger-1]
  y2 <- y[first_larger]
  y1 <- y[first_larger-1]
  slope <- (y2 - y1)/(x2 - x1)
  thin <- x1 + (y.crit - y1)/slope
  
  return(thin)
}

#' KL divergence
#' 
#' Calculates Kullbeck-Leibler divergence KL(p||q) for discrete p and q.
#'
#' @param p probability mass function
#' @param q probability mass function
#' @return KL divergence.
#' @export
KL <- function(p,q) {
  summand <- p * log(q/p)
  # 0 * log(0) is 0 for KL purposes
  summand[p < .Machine$double.eps] <- 0
  return(-sum(summand))
}

expandDistanceMatrix <- function(dmat,indices) {
  # recover()
  
  # When computing the variance of the posterior for different ESS, NAs get introduced into the indices
  indices <- indices[!is.na(indices)]
  
  ntrees <- dim(dmat)[1]
  ngen <- length(indices)
  
  # If all indices are the same index, we know the new distance matrix immediately
  if ( length(unique(indices)) == 1 ) {
    return(matrix(0,ngen,ngen))
  }
  
  # Label the matrix so we know which tree is represented in each row/column
  row.names(dmat) <- 1:ntrees
  colnames(dmat) <- 1:ntrees
  
  # Delete all unneeded rows/columns
  to_remove <- sort(which(sapply(1:ntrees,function(i){sum(indices == i) == 0})),decreasing=TRUE)
  for (i in to_remove) {
    dmat <- dmat[,-i]
    dmat <- dmat[-i,]
  }
  
  n_used_trees <- dim(dmat)[1]
  unique_indices <- as.numeric(colnames(dmat))
  
  # Duplicate all the rows/columns in order
  for (i in 1:n_used_trees) {
    count <- sum(indices == unique_indices[i])
    # If count == 1, the matrix is unchanged for now. Otherwise we must expand the entries here
    if (count > 1) {
      current_dim <- dim(dmat)[1]
      new_dim <- current_dim + count - 1
      
      tmp <- matrix(nrow=new_dim,ncol=new_dim)
      
      # Expanding the matrix breaks it up into three chunks (some may have size 0), left, center, and right (or upper, middle, and lower)
      current_indices <- as.numeric(colnames(dmat))
      
      idx <- unique_indices[i]
      
      left_current <- which(current_indices < idx)
      center_current <- which(current_indices == idx)
      right_current <- which(current_indices > idx)
      
      new_names <- c(current_indices[left_current],rep(current_indices[center_current],count),current_indices[right_current])
      
      left_new <- left_current
      center_new <- center_current:(center_current+count-1)
      right_new <- c()
      if ( max(center_new) < new_dim ) {
        right_new <- (max(center_new)+1):(max(center_new)+length(right_current))
      }
      
      # Fill upper-left portion (already expanded) and get first chunk of central portion (the portion to expand)
      first <- c()
      if ( length(left_current) > 0 ) {
        tmp[left_new,left_new] <- dmat[left_current,left_current]
        first <- dmat[left_new,center_current]
      }
      
      # Fill bottom-right portion (not yet expanded) and get last chunk of central portion (the portion to expand)
      last <- c()
      if ( length(right_current) > 0 ) {
        tmp[right_new,right_new] <- dmat[right_current,right_current]
        last <- dmat[right_current,center_current]
      }
      
      # Fill upper-right and bottom-left portions
      if ( (length(left_current) > 0) && (length(right_current) > 0) ) {
        # upper-right
        tmp[left_new,right_new] <- dmat[left_current,right_current]
        # bottom-left
        tmp[right_new,left_new] <- dmat[right_current,left_current]
      }
      
      # Fill newly-expanded portion
      new_row <- c(first,rep(0,count),last)
      tmp[,center_new] <- new_row
      tmp <- t(tmp)      
      tmp[,center_new] <- new_row
      
      # Done
      dmat <- tmp
      colnames(dmat) <- new_names
      row.names(dmat) <- new_names
    }
  }
  # Re-order the matrix to match the indices
  blocked_indices <- as.numeric(colnames(dmat))
  blocked_order <- rank(blocked_indices,ties.method="first")
  target_order <- rank(indices,ties.method="first")
  key <- match(target_order,blocked_order)
  dmat <- dmat[,key]
  dmat <- dmat[key,]
  
  return(dmat)
}

# Finds all unique trees (ignoring branch lengths) in the input set
uniqueTopologies <- function(trees) {
  
  remaining_trees <- trees
  
  topos <- vector("list")
  counts <- numeric()
  
  iter <- 0
  while( length(remaining_trees) > 0 ){
    iter <- iter + 1
    # Strip branch lengths
    clado <- remaining_trees[[1]]
    clado$edge.length <- NULL
    
    # Add this topology
    topos[[iter]] <- clado
    
    # count matches
    dists <- suppressWarnings(RF.dist(remaining_trees,clado))
    counts <- c(counts,sum(dists == 0))
    
    # Remove matches
    remaining_trees <- remaining_trees[dists > 0]
  }
  
  class(topos) <- "multiPhylo"
  
  return(list(topologies=topos,counts=counts))
}

# Finds all unique trees (with branch lengths) in the input set
uniqueTrees <- function(trees) {
  
  remaining_trees <- trees
  
  utrees <- vector("list")
  counts <- numeric()
  
  iter <- 0
  while( length(remaining_trees) > 0 ){
    iter <- iter + 1
    
    # Add this topology
    utrees[[iter]] <- remaining_trees[[1]]
    
    # count matches
    dists <- suppressWarnings(KF.dist(remaining_trees,remaining_trees[[1]]))
    counts <- c(counts,sum(dists == 0))
    
    # Remove matches
    remaining_trees <- remaining_trees[!(dists == 0)]
  }
  
  class(utrees) <- "multiPhylo"
  
  return(list(trees=utrees,counts=counts))
}
