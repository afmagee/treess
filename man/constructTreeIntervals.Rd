% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructTreeIntervals.R
\name{constructTreeIntervals}
\alias{constructTreeIntervals}
\title{Constructs intervals for summaries of trees.}
\usage{
constructTreeIntervals(
  x,
  ESS,
  type = "confidence",
  interval.width = 0.95,
  method = "Jeffreys"
)
}
\arguments{
\item{x}{A list of chains containing the MCMC samples of trees as multiPhylo objects.}

\item{ESS}{Vector of ESS of each of the MCMC chains in x.}

\item{type}{Either "confidence" for confidence intervals or "prediction" for prediction intervals.}

\item{interval.width}{Width of the interval, value in (0,1), 1 - 2*alpha. 0.95 for 95\% confidence/prediction intervals.}

\item{method}{For confidence intervals, method of interval constriction, "Jefreys"|"Wilson"|"ContinuityCorrectedWilson".}
}
\value{
A list of intervals for each chain, see details.
}
\description{
For one or more MCMC chains, constructs confidence/prediction intervals for tree and split probabilities.
}
\details{
The function returns a list of lists.
The top level of lists is $split and $topology, which separates split and topology probabilities.
Within each of these is a list, each element containing the intervals for each chain as a matrix.
Each matrix has splits (or trees) in rows (rows are comparable across chains), and the columns are the point estimate, lower, and upper CIs, and the effective sample size.
Note that prediction intervals are generated once per chain, and do not account for differences in per-chain ESS.
That is, they use the default n.new=n in \link{binomialProportionPI}.
}
\seealso{
\link{binomialProportionCI}, \link{binomialProportionPI}, \link{treeStability}, \link{plotTreeIntervals}
}
