% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treeStability.R
\name{plotTreeStabilityConvergence}
\alias{plotTreeStabilityConvergence}
\title{Visualizes stability of the tree estimates from a sample of trees over the length of the run.}
\usage{
plotTreeStabilityConvergence(
  tree.stability.convergence,
  stat,
  colors = NA,
  CI.color.alpha = 0.5,
  plot.median = TRUE,
  use.ess = TRUE,
  xlab = NA,
  ylab = NA,
  main = NA,
  ...
)
}
\arguments{
\item{tree.stability.convergence}{Output of treeStabilityConvergence.}

\item{stat}{Which summary of tree stability to ploe? ASDSF|MRC|TOPOPROBS}

\item{colors}{Colors for each of the consensus thresholds for plotting. Any provided transparencies will be removed. If NA defaults are used (defaults not guaranteed to be pretty).}

\item{CI.color.alpha}{Transparency value for plotting CI thresholds}

\item{plot.median}{Is there a median line to be plotted? Only works if 0.5 is included in treeStabilityConvergence(probs=c(...,0.5,...)).}

\item{use.ess}{If TRUE, and there are computed ESS values in tree.stability.convergence, plots ESS for x-axis. Otherwise x-axis is raw number of samples.}

\item{xlab}{The x-axis label, NA for defaults.}

\item{ylab}{The y-axis label, NA for defaults.}

\item{main}{The plot title.}

\item{...}{Further arguments to be passed to plotting functions.}
}
\value{
Nothing, plots the diagnostic.
}
\description{
Visualizes stability of the tree estimates from a sample of trees over the length of the run.
}
\details{
Note that the option use.ess is purely cosmetic, as ESS is only computed for the unbootstrapped subchain.
If choosing to plot the x-axis based on ESS, and the ESS is not estimated to be monotonically increasing in subchain length, sensible plotting results cannot be guaranteed.
}
\seealso{
\link{treeStability},\link{treeStabilityConvergence}
}
