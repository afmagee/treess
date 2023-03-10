% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructAdjacencyGraph.R
\name{constructAdjacencyGraph}
\alias{constructAdjacencyGraph}
\title{Computes the NNI adjacency graph for a set of (unique) tree topologies.}
\usage{
constructAdjacencyGraph(trees, weights = rep(1, length(trees)), trim = FALSE)
}
\arguments{
\item{trees}{A multiPhylo object or list of phylo objects}

\item{weights}{The probabilities for each tree (need not sum to one)}

\item{trim}{Should we trim to the largest connected subset of trees in the posterior?}
}
\value{
A list. $trees contains all topologies, $probs their probabilities, $NNI the NNI connectivity graph.
}
\description{
Computes the NNI adjacency graph for a set of (unique) tree topologies.
}
\examples{
constructAdjacencyGraph(rmtree(100,5))
}
