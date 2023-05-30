areMatricesEquivalentWithReordering <- function(x,y,tolerance=testthat::testthat_tolerance()) {
  # recover()
  d <- dim(x)[2]
  if ( dim(y)[2] != d ) {
    return(FALSE)
  }
  
  n <- dim(x)[1]
  foundAll <- TRUE
  for (i in 1:n) {
    where <- 0
    foundThis <- FALSE
    m <- dim(y)[1]
    for (j in 1:m) {
      if ( all.equal(x[1,],y[j,],tolerance=tolerance) == TRUE ) {
        foundThis <- TRUE
        where <- j
        break
      }
    }
    if ( !foundThis ) {
      foundAll <- FALSE
      break
    }
    if ( m == 2 ) {
      x <- matrix(x[-1,],ncol=d,nrow=1,byrow=TRUE)
      y <- matrix(y[-where,],ncol=d,nrow=1,byrow=TRUE)
    } else {
      x <- x[-1,]
      y <- y[-where,]
    }
  }
  
  if ( dim(x)[1] != 0 || dim(y)[1] != 0 ) {
    foundAll <- FALSE
  }
  
  return(foundAll)
}
