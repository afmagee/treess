test_that("Compare simulation study output values to cached values.", {
  
  # Get stored values
  file <- system.file("testdata","precomputed_errors.rda",package="treess")
  load(file)

  # Get adjacency graphs
  data("ds.adjacency.graphs")
  
  # run simulation pipeline
  simulated <- withr::with_seed(42,simulatePhylogeneticMCMC(ds.adjacency.graphs$ds3.adjacency.graph,ngen=500,nchains=10,verbose=FALSE))
  error <- withr::with_seed(42,effectiveSizeEquivalentError(simulated,"rf",ess.methods="medianPseudoESS",verbose=FALSE))
  
  testthat::expect_equal(error$medianPseudoESS$ESS,precomputed.error$medianPseudoESS$ESS)
  testthat::expect_equal(error$medianPseudoESS$MRCSquaredError,precomputed.error$medianPseudoESS$MRCSquaredError)
  testthat::expect_equal(error$medianPseudoESS$treeProbSquaredError,precomputed.error$medianPseudoESS$treeProbSquaredError)
  # testthat::expect_equal(error$medianPseudoESS$splitProbSquaredError,precomputed.error$medianPseudoESS$splitProbSquaredError)
  
  testthat::expect_equal(error$logPosteriorESS$ESS,precomputed.error$logPosteriorESS$ESS)
  testthat::expect_equal(error$logPosteriorESS$MRCSquaredError,precomputed.error$logPosteriorESS$MRCSquaredError)
  testthat::expect_equal(error$logPosteriorESS$treeProbSquaredError,precomputed.error$logPosteriorESS$treeProbSquaredError)
  # testthat::expect_equal(error$logPosteriorESS$splitProbSquaredError,precomputed.error$logPosteriorESS$splitProbSquaredError)
  
  # Ordering of splits is essentially arbitrary, and distinct between trees2coords and perTreeSplits
  # Here we check for equivalence ignoring order
  testthat::expect_true(areMatricesEquivalentWithReordering(error$medianPseudoESS$splitProbSquaredError,precomputed.error$medianPseudoESS$splitProbSquaredError,tolerance=1e-10))
  testthat::expect_true(areMatricesEquivalentWithReordering(error$logPosteriorESS$splitProbSquaredError,precomputed.error$logPosteriorESS$splitProbSquaredError,tolerance=1e-10))
  
})