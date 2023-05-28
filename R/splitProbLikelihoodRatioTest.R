#' Examine differences in split probabilities among chains
#' 
#' Takes split (or clade) probabilities from multiple MCMC runs and compares them using a likelihood ratio test, possibly accounting for differing effective sample sizes among chains.
#'
#' @param x A list of: trees, perTreeSplits, or split probability vectors.
#' @param ESS The effective sample sizes for the MCMC chains in x, or a target ESS cutoff as in Fabreti and Höhna (2022) (see details).
#' @param rooted If trees are being fed into x, should splits (FALSE) or clades (TRUE) be retrieved.
#' @param log.p If TRUE and return.p.value == TRUE, log(p) is returned instead of p.
#' @return The p-values from a LRT for differences in proportion for each observed split.
#' @references 
#' Fabreti and Höhna (2022). "Convergence assessment for Bayesian phylogenetic analysis using MCMC simulation." Methods in Ecology and Evolution, 13(1), 77-90.
#' 
#' Magee, Karcher, Matsen, and Minin (2023). "How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error." Bayesian Analysis, 1(1), 1-29.
#' @details 
#' For a given split (or clade), a LRT is conducted testing the null hypothesis that the probability of the split is the same in all chains.
#' (Splits are tested independently.)
#' Instead of using the raw number of times a split is seen in a particular chain, an adjusted count is used to accommodate Monte Carlo error.
#' Specifically, for a split with observed frequency p in a chain with ESS effective samples of the phylogeny, the adjusted count used is p * ESS.
#' 
#' If the actual effective sample size is used, this test generalizes the confidence-interval based approach used in \link{binomialProportionDifferenceCI} to >2 chains.
#' 
#' Alternately, a fixed sample size can be fed in, such as the 625 suggested by Fabreti and Höhna (2022).
#' In this case, the acceptable difference in split frequencies is fixed, and as the chains are run longer (targeting the same mode) eventually all p-values will exceed a given cutoff.
#' This is similar to the approach used in convenience::plotDiffSplits.
#' 
#' @export
#' @seealso \link{compareChainProbabilities}, \link{binomialProportionDifferenceCI}, \link{ASDSF}
#' @examples
#' \dontrun{
#' # Load data, compute ESS
#' data(ds.mcmc)
#' chains <- list(ds.mcmc$ds3.low,ds.mcmc$ds3.med,ds.mcmc$ds3.high)
#' ess <- treess(chains,phangorn::RF.dist,methods="frechetCorrelationESS")
#' 
#' # Examine differences in split probabilities accounting for Monte Carlo error
#' pvals.ess <- splitProbLikelihoodRatioTest(chains,as.numeric(ess))
#' hist(pvals.ess,breaks=seq(0,1,0.05),xlab="p-value",main="Comparison of split probabilities in DS3")
#' 
#' # Examine differences in split probabilities more like convenience::plotDiffSplits
#' alpha <- 0.05 # a p-value threshold
#' pvals.fixed <- splitProbLikelihoodRatioTest(chains,625,log.p=TRUE)
#' split.probs <- splitProbs(perTreeSplits(do.call(c,chains)))
#' key <- match(names(split.probs),names(pvals.fixed))
#' plot(split.probs,-pvals.fixed[key],xlab="split frequency",ylab="-log(LRT p-value)",main="incongruent splits are those above the red line")
#' abline(h=-log(alpha),col="red")
#' }

splitProbLikelihoodRatioTest <- function(x,ESS,rooted=FALSE,log.p=FALSE) {
  # recover()
  is_same_tolerance <- .Machine$double.xmin
  nchains <- length(x)
  xclasses <- sapply(x,class)
  split_probs <- NULL
  # If numerics fed in, make sure we can work with them
  if ( all(xclasses == "numeric") ) {
    split_probs <- x
    lens <- lengths(split_probs)
    is_named <- sapply(split_probs,function(x){!is.null(names(x))})
    if ( !all(is_named) ) {
      if ( length(unique(lens)) > 1 ) {
        stop("Cannot match uneven length split probability vectors.")
      } else {
        for (i in 1:length(split_probs)) {
          names(split_probs[[i]]) <- 1:lens[1]
        }
      }
    }
  } else { # If trees or split probabilities fed in, process them
    if ( all(xclasses == "multiPhylo") || all(xclasses == "list") && all(sapply(x[[1]],class) == "phylo") ) {
      x <- lapply(x,perTreeSplits)
    } else if ( !all(xclasses == "perTreeSplits") ) {
      stop("Invalid input to argument \"x\".")
    }
    split_probs <- lapply(x,splitFrechetMean)
  }
  
  if ( length(ESS) != length(x) ) {
    if ( length(ESS) == 1 ) {
      warning("Single value provided for ESS taken to represent all chains.")
      ESS <- rep(ESS,length(x))
    } else {
      stop("Number of ESS values provided does not match number of chains.")
    }
  }
  
  all_splits <- unique(unlist(lapply(split_probs,names)))
  stats <- sapply(all_splits,function(split){
    probs <- sapply(split_probs,function(chain){
      where <- which(names(chain) == split)
      if ( length(where) == 1 ) {
        return(chain[where])
      } else if ( length(where) == 0 ) {
        return(0)
      } else {
        stop("Severe error, same split found multiple times in chain.")
      }
    })
    adjusted_counts <- probs * ESS
    prob <- sum(adjusted_counts)/sum(ESS)
    
    # log likelihoods up to the combinatorial factor, which cancels
    # splits which are absent or fixed in one or more chains contribute 0 * log(0) terms which become NaN
    # we take 0^0 to be 1, so we force 0 * log(0) to be 0 here
    lnL0 <- sum(adjusted_counts) * log(prob) + (sum(ESS) - sum(adjusted_counts)) * log(1.0 - prob)
    if ( abs(sum(adjusted_counts) - sum(ESS)) < is_same_tolerance || sum(adjusted_counts) < is_same_tolerance ) {
      lnL0 <- 0.0
    }
    summand <- adjusted_counts * log(probs) + (ESS - adjusted_counts) * log(1.0 - probs)
    summand[abs(adjusted_counts - ESS) < is_same_tolerance | adjusted_counts < is_same_tolerance] <- 0.0
    lnL1 <- sum(summand)
    return(-2 * (lnL0 - lnL1))
  })
  
  results <- pchisq(stats,df=nchains-1,lower.tail=FALSE,log.p=log.p)

  return(results)
}