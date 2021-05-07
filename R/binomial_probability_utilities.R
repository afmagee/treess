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
#' @references 
#' Newcombe (1998). "Two-sided confidence intervals for the single proportion: comparison of seven methods." Statistics in Medicine.
#' Brown, Cai, and DasGupta (2001). "Interval estimation for a binomial proportion." Statistical science.
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
    ub[is_1] <- 1
  } else {
    stop("Invalid \"method\"")
  }
  ci <- cbind(lb,ub)
  colnames(ci) <- paste0(round(100*c(a_2,one_minus_a_2),1),"%")
  return(ci)
}

#' Confidence intervals for differences in proportions.
#' 
#' Provides two methods for computing confidence intervals for binomial proportions.
#' Additionally allows for any method available in DescTools::BinomDiffCI.
#' Useful for turning estimates of the treESS into CIs for differences between split and topology probabilities between chains.
#'
#' @param p1 The (estimated) probability in sample 1.
#' @param p2 The (estimated) probability in sample 2.
#' @param n1 The sample size (or effective sample size) in sample 1.
#' @param n2 The sample size (or effective sample size) in sample 2.
#' @param method Either AgrestiCoffo|JeffreysPerks or an option for package DescTool's function BinomDiffCI (if installed).
#' @param ci.width Width of the desired confidence interval (value in [0,1]). Equivalent to 1 - 2*alpha.
#' @return A length(p)x2 matrix containing the lower and upper ends of the CI.
#' @export
#' @details There are a variety of options for differences in binomial proportions with a variety of performance options. DescTools recommends mn, though Agresti-Coffo and Jeffreys-Perks nominally perform well (see references).
#' @references 
#' Agresti, A., & Caffo, B. (2000). Simple and effective confidence intervals for proportions and differences of proportions result from adding two successes and two failures. The American Statistician, 54(4), 280-288.
#' Newcombe, R. G. (1998). Interval estimation for the difference between independent proportions: comparison of eleven methods. Statistics in medicine, 17(8), 873-890.
#' @seealso \link{binomialProportionPI}, DescTools::BinomDiffCI
#' @examples
#' 
#' p1 <- rbeta(100,1,1) # Generate probabilities
#' p2 <- p1 * rbeta(100,5,1) # Generate probabilities
#' n1 <- n2 <- 100 # Sample size
#' x1 <- rbinom(length(p1),n1,p1) # Random variables
#' x2 <- rbinom(length(p2),n2,p2) # Random variables
#' 
#' aci <- binomialProportionDifferenceCI(x1/n1,x2/n2,rep(n1,100),rep(n2,100),"AgrestiCoffo")
#' jci <- binomialProportionDifferenceCI(x1/n1,x2/n2,rep(n1,100),rep(n2,100),"JeffreysPerks")
#' 
#' sum(aci[,1] < 0 & aci[,2] > 0)
#' sum(jci[,1] < 0 & jci[,2] > 0)


binomialProportionDifferenceCI <- function(p1,p2,n1,n2,method,ci.width=0.95) {
  # recover()
  if ( !(length(p1) == length(p2) && length(p1) == length(n1) && length(p1) == length(n2))) {
    stop("Inputs p1,p2,n1,n2 must all be same length.")
  }
  
  # CI bounds
  lb <- ub <- numeric(length(p1))
  # Ends of CIs in probability
  a_2 <- (1-ci.width)/2
  one_minus_a_2 <- 1 - a_2
  
  # estimate the difference
  delta <- p1 - p2
  
  # Normal quantile at p=alpha/2
  z <- qnorm(one_minus_a_2,0,1)
  
  # Get ESS-adjusted counts
  x1 <- p1*n1
  x2 <- p2*n2
  
  if ( (grepl("jef",tolower(method))) ) {
    psi <- 0.5 * ((x1+0.5)/(n1+1) + (x2+0.5)/(n2+1))
    u <- (1/n1 + 1/n2)/4
    v <- (1/n1 - 1/n2)/4
    
    theta_star <- (delta + (z^2)*v*(1-2*psi))/(1 + z^2 * u)
    
    w <- z/(1 + (z^2)*u) * sqrt( (u*(4*psi*(1-psi)-delta^2)) + (2*v*(1-2*psi)*delta) + (4*(z^2)*(u^2)*(1-psi)*psi) + ((z^2)*(v^2)*((1-2*psi)^2)))
    
    lb <- theta_star - w
    ub <- theta_star + w
  } else if ( (grepl("agr",tolower(method))) ) {
    # Agresti's and Coffo's confidence interval
    p1_adj <- (x1+1)/(n1+2)
    p2_adj <- (x2+1)/(n2+2)

    lb <- (p1_adj - p2_adj) - z * sqrt( (p1_adj*(1-p1_adj))/(n1+2) + (p2_adj*(1-p2_adj))/(n2+2) )
    ub <- (p1_adj - p2_adj) + z * sqrt( (p1_adj*(1-p1_adj))/(n1+2) + (p2_adj*(1-p2_adj))/(n2+2) )
  } else if ( tolower(method) %in% c("ac", "wald", "waldcc", "score", "scorecc", "mn", "mee", "blj", "ha", "hal", "jp") ) {
    if (!requireNamespace("DescTools",quietly=TRUE)) {
      stop(paste0("Cannot use metho ",method," without DescTools installed."))
    }
    tmp <- DescTools::BinomDiffCI(x1,n1,x2,n2,conf.level=ci.width,method=method,sides="two.sided")
    lb <- tmp[,2]
    ub <- tmp[,3]
  } else {
    stop("Invalid \"method\"")
  }
  ci <- cbind(lb,ub)
  colnames(ci) <- paste0(round(100*c(a_2,one_minus_a_2),1),"%")
  return(ci)
}
