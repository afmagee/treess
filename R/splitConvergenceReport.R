#' Convergence diagnostic via differences in split probabilities
#'
#' The
#'
#' @param x One or more samples of trees (presumably MCMC chains), possibly processed as splits. See details.
#' @param frac For single-chain diagnostic, controls size of subsamples
#' @details
#' A variety of inputs are accepted.
#' A single sample of trees (presumably an MCMC output) can be fed in as a multiPhylo object, multiple sets of samples may be fed in as a list of multiPhylo objects.
#' Pre-processed perTreeSplits (or lists thereof) are acceptable.
#' When providing multiple chains, all must have the same length.
#' The single-chain diagnostic is performed by dividing the chain into nsubchains sub chains of size floor(length(x) / n).
#' (This will result in the last length(x) %% nsubchains trees going unused by the diagnostic.)
#'
#' @export
MNSS <- 

reportMNSS <- function() {
  
}


.mnss <- function(x, nsubchains = 4, alpha = 0.05, rooted = FALSE) {
  # splits <- .sanitizeInputsMNSS(x, frac, rooted)
  recover()

  if (length(nsubchains) > 1 || !is.numeric(nsubchains)) {
    stop("Invalid supplied \"nsubchains\".")
  }

  n <- 0
  nchains <- 0
  single_chain <- NA
  if (length(x) == 1) {
    n <- floor(length(x[[1]]) / nsubchains)
    nchains <- nsubchains
    x <- x[[1]]
    single_chain <- TRUE
  } else {
    if (length(unique(lengths(x))) != 1) {
      stop("Provided chains must all be same length.")
    }
    nchains <- length(x)
    n <- length(x[[1]])
    x <- do.call(c, x)
    single_chain <- FALSE
  }

  probs <- lapply(1:nchains, function(k){
    # TODO: replace with splitFrechetMean?
    treess::splitProbs(x[(k - 1) * n + (1:n)])
  })

  helper <- treess:::perTreeSplitsConvergenceHelper(probs, return.gcomp = TRUE, return.mat = TRUE)

  gtilde <- helper$gcomp

  gcrit <- pchisq(nchains - 1, 1 - alpha)
  
  nss <- gcrit / gtilde
  
  spans <- apply(helper$mat, 2, function(f) {
    max(f) - min(f)
  })
  
  trends <- rep(NA, length(nss))
  if (single_chain) {
    trends <- helper$mat[nsubchains,] - helper$mat[1,]
  }
  
  sdsf <- helper$sdsf
  
}

.sanitizeInputsMNSS <- function(x, frac = 0.25, rooted = FALSE) {
  if (class(x) == "multiPhylo") {
    return(list(perTreeSplits(x, rooted)))
  } else if (class(x) == "perTreeSplits") {
    return(list(x))
  } else if (class(x) == "list") {
    stop("Not yet implemented")
  } else {
    stop("Unrecognized input to \"x\".")
  }
}