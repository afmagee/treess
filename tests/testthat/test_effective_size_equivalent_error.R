test_that("Compare simulation study output values to cached values.", {
  
  # Get stored values
  file <- system.file("testdata","precomputed_errors.rda",package="treess")
  load(file)

  # Get adjacency graphs
  data("ds.adjacency.graphs")
  
  # run simulation pipeline
  simulated <- withr::with_seed(42,simulatePhylogeneticMCMC(ds.adjacency.graphs$ds3.adjacency.graph,ngen=500,nchains=10,verbose=FALSE))
  error <- withr::with_seed(42,effectiveSizeEquivalentError(simulated,"rf",ess.methods="medianPseudoESS",verbose=FALSE))

  testthat::expect_equivalent(error,precomputed.error)
})