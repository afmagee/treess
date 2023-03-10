% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareChainProbabilities.R
\name{compareChainProbabilities}
\alias{compareChainProbabilities}
\title{Constructs intervals for differences (between chains) of summaries of trees.}
\usage{
compareChainProbabilities(
  x,
  ESS,
  interval.width = 0.95,
  method = "AgrestiCoffo"
)
}
\arguments{
\item{x}{A list of chains containing the MCMC samples of trees as multiPhylo objects.}

\item{ESS}{Vector of ESS of each of the MCMC chains in x.}

\item{interval.width}{Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95\% confidence intervals.}

\item{method}{Method of confidence interval constriction, "AgrestiCoffo"|"JeffreysPerks" or anything available in DescTools::BinomDiffCI.}
}
\value{
A list of intervals for each chain, see details.
}
\description{
For two or more MCMC chains, constructs confidence intervals for differences between tree and split probabilities.
}
\details{
The function returns a list of lists.
The top level of lists is $split and $topology, which separates split and topology probabilities.
Within each of these is a list, each element containing the intervals for each chain as a matrix.
Each matrix has splits (or trees) in rows (rows are comparable across chains), and the columns are the point estimate, lower, and upper CIs, and the effective sample size.
}
\seealso{
\link{binomialProportionDifferencePI}, \link{plotTreeIntervals}
}
