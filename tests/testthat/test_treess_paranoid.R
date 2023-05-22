test_that("Compare computed ESS values to cached values.", {

  # Get stored ESS values (all approaches, )
  file <- system.file("testdata","precomputed_ess.csv",package="treess")
  precomputed <- as.matrix(read.csv(file,row.names=1))

  # Get MCMC chains, will load to be an object called ds.mcmc
  data("ds.mcmc")

  # Compute ESS
  ess <- treess(ds.mcmc,phangorn::RF.dist,methods=getESSMethods(FALSE))
  ess <- do.call(cbind,ess)

  # Check all ESS methods
  for (i in 1:dim(ess)[1]) {
    ess_method <- rownames(precomputed)[i]
    expect_equal(ess[i,],precomputed[i,],tolerance=1e-10,label=paste0("Newly computed ESS values computed using ",ess_method),expected.label=paste0("stored ESS values using ",ess_method))
  }
})