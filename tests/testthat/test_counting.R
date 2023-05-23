test_that("Test tree counting.", {
  # Count the trees
  data("ds.mcmc")
  counts <- lapply(ds.mcmc,countTrees)
  
  # Make sure there are as many trees as expected
  ntrees <- as.integer(sapply(counts,function(x){length(x$trees)}))
  expect_identical(ntrees,c(28L,245L,680L,
                             10L,55L,55L,
                             66L,823L,835L))
  
  # Make sure there are as many counts as trees
  ncounts <- as.integer(sapply(counts,function(x){length(x$counts)}))
  expect_identical(ncounts,ntrees)
  
  # Make sure we counted all the trees
  ntot <- as.integer(sapply(counts,function(x){sum(x$counts)}))
  expect_identical(ntot,rep(1000L,9))
  
  # Make sure we counted the trees _right_
  expect_identical(as.integer(counts[[1]]$counts),c(218L, 202L, 75L, 65L, 48L, 43L, 41L, 
                                                    37L, 30L, 30L, 25L, 25L, 21L, 20L, 
                                                    18L, 17L, 12L, 12L, 11L, 11L, 7L, 
                                                    6L, 6L, 6L, 5L, 4L, 3L, 2L))
})