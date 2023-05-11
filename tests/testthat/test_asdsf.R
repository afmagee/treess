test_that("Test ASDSF from tree inputs.", {
  data("ds.mcmc")
  
  asdsf <- ASDSF(list(ds.mcmc$ds3.high[1:500],ds.mcmc$ds3.high[501:1000]))
  expect_equal(asdsf,0.0046044162495868,tolerance=1e-10)

  msdsf <- ASDSF(list(ds.mcmc$ds3.high[1:500],ds.mcmc$ds3.high[501:1000]),summary.fn=max)
  expect_equal(msdsf,0.0268700576850888,tolerance=1e-10)
  
  weighted_asdsf <- ASDSF(list(ds.mcmc$ds3.high[1:500],ds.mcmc$ds3.high[501:1000]),weights="entropy")
  expect_equal(weighted_asdsf,0.0111919464205538,tolerance=1e-10)
  
  higher_cutoff_asdsf <- ASDSF(list(ds.mcmc$ds3.high[1:500],ds.mcmc$ds3.high[501:1000]),min.freq=0.1)
  expect_equal(higher_cutoff_asdsf,0.0038221988172246,tolerance=1e-10)
})