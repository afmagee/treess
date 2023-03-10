% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/countTrees.R
\name{countTrees}
\alias{countTrees}
\title{Counts occurences of each tree in a set of trees, returns unique trees and their counts}
\usage{
countTrees(trees)
}
\arguments{
\item{trees}{A multiPhylo object or list of phylo objects}
}
\value{
A list. $trees contains all unique topologies, and $counts their counts
}
\description{
Counts occurences of each tree in a set of trees, returns unique trees and their counts
}
\details{
Allowable methods are "pseudoPSRF", "totalDistancePSRF", "foldedRankMedioidPSRF", and "medianDistanceWilcoxTest"
}
\examples{
countTrees(rmtree(100,10))
}
