# treess

A package for assessing the MCMC effective sample size of phylogenetic trees.

Installation requires the package `devtools`:

    devtools::install_github("afmagee/treess")

Note that while treess depends on `ape`, `coda`, and `phangorn`, it does not load any of them.

## Citation
If using `treess`, please cite: Magee, A F, Karcher, M D, Matsen, E M, and Minin, V M. How trustworthy is your tree? Bayesian phylogenetic effective sample size through the lens of Monte Carlo error. Bayesian Analysis, 1(1), 1-29.

If using the approximateESS, medianPseudoESS, or minPseudoESS, please also cite: Lanfear, R, Hua, X, & Warren, D L (2016). Estimating the effective sample size of tree topologies from Bayesian phylogenetic analyses. Genome Biology and Evolution, 8(8), 2319-2332.

## Basic use
The basic purpose of the `treess` package is to compute the effective sample size of a set of MCMC samples of a phylogenetic tree (or other non-Euclidean object).
This is accomplished using the `treess()` function, which is called on a list of MCMC chains and will compute the ESS for each chain separately.
The following pseudo-code may be useful.

    library(treess)
    library(ape)
    library(phangorn)

    # Read in samples
    tree_samples_chain_1 <- read.tree("/path/to/chain1.trees")
    tree_samples_chain_2 <- read.tree("/path/to/chain2.trees")

    # List all chains
    all_chains <- list(tree_samples_chain_1,tree_samples_chain_2)

    # Compute all ESS measures for all chains
    all_ess_measures <- treess(all_chains,dist.fn=RF.dist)

To get the ESS of the concatenated MCMC chain, the following code may be helpful.

    # Concatenate chains
    all_chains_concatenated <- unlist(all_chains,recursive=FALSE)
    class(all_chains_concatenated) <- "multiPhylo"

    # Compute all ESS measures for concatenated chain
    all_ess_measures_concatenated <- treess(list(all_chains_concatenated),dist.fn=RF.dist)

## Comparing multiple chains with confidence

Tree ESS measures can be used to construct confidence intervals on split probabilities (and differences between probabilities).
This can be useful for plotting split probabilities to compare the estimates between two chains.

    # Extract a specific ESS measure
    per_chain_frechet_ess <- unlist(lapply(all_ess_measures,function(x){x$frechetCorrelationESS}))

    # Compute CIs for each run
    tree_intervals <- constructTreeIntervals(all_chains,per_chain_frechet_ess,type="confidence")

    # Examine differences in split probabilities between chains (used to color CIs in plot, not strictly needed)
    chain_comparisons <- compareChainProbabilities(all_chains,per_chain_frechet_ess)

    # Plot!
    plotTreeIntervals(tree_intervals,chain_comparisons)

## Advanced use
The `treess` package also allows for testing ESS measures and generating MCMC samples of tree topologies on known posterior distributions.
For examples of how to do this, please see [this repository](https://bitbucket.org/afmagee/tree_convergence_code/).
