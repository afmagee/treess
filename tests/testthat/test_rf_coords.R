test_that("Test getting splits from DS3 best trees.", {
  
  # Get stored values
  file <- system.file("testdata","ds3_as_split_matrix.csv",package="treess")
  precomputed <- read.csv(file,check.names=FALSE)
  
  # Compute splits and reorder/reindex to match
  data("ds.adjacency.graphs")
  phy <- ds.adjacency.graphs$ds3.adjacency.graph$trees
  rc <- as.RFcoords(phy)
  
  ntaxa <- ape::Ntip(phy[[1]])
  taxa <- 1:ntaxa
  
  splits <- rc$taxa
  splits <- lapply(splits,as.integer)
  splits <- lapply(splits,function(split){
    bitsplit <- rep(0,ntaxa)
    bitsplit[split] <- 1
    if ( bitsplit[1] == 1 ) {
      bitsplit <- 1 - bitsplit
    }
    bitsplit <- as.logical(bitsplit)
    return(taxa[bitsplit])
  })
  splits <- lapply(splits,paste0,collapse=",")
  splits <- unlist(splits)
  key <- order(splits)
  
  x <- rc$coords[,key]
  
  colnames(x) <- gsub(",",";",splits[order(splits)])
  
  testthat::expect_equivalent(colnames(x),colnames(precomputed))
  testthat::expect_equivalent(x,as.matrix(precomputed))
})