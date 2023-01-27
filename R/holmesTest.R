#' Holmes' 2-sample tree test
#' 
#' Holmes' tree test is a permutation test of equivalence of two (or more) tree distributions.
#' The input are sets of samples of trees, from which a minimum spanning tree is constructed to link all trees.
#' The test statistic is the number of edges which connect trees from the same set of samples.
#' The null distribution is obtained by permuting the labels (the assignment of trees to sets of samples).
#' 
#' @param x Either a list of multiPhylo objects or a distance matrix for distances between trees, if the latter argument "labels" is needed.
#' @param dist.fn Function for computing tree to tree distances, must be capable of being called on a multiPhylo. Not needed if x is a distance matrix.
#' @param labels If x is a distance matrix, this argument assigns each column and row to a bootstrap replicate.
#' @param B number of bootstrap replicates.
#' @param setColor colors for treesets for plotting.
#' @param betweenColor color to distinguish edges between sets for plotting.
#' @return See details.
#' @details 
#' holmesTest returns a list with $p.value the (estimated) p-value from the permutation test, $S_0 the observed test statistic, and optionally $nullDistribution the null distribution of the test statistic obtained by permutation.
#' 
#' holmesPlot returns NULL and plots the visualization.
#' @export
#' @references 
#' Holmes (2005). "Statistical approach to tests involving phylogenies." In Mathematics of Evolution and Phylogeny, ed. Gascuel.
#' @examples
#' \dontrun{
#' # Get a bunch of trees from two distinct distributions
#' true1 <- ape::rtree(10,FALSE,br=rexp,rate=170)
#' true2 <- ape::rtree(10,FALSE,br=rexp,rate=170)
#' 
#' set1 <- lapply(1:100,function(i){
#'   aln <- phangorn::simSeq(true1,l=200)
#'   dm <- phangorn::dist.ml(aln)
#'   ape::bionj(dm)
#' })
#' class(set1) <- "multiPhylo"
#' 
#' set2 <- lapply(1:100,function(i){
#'   aln <- phangorn::simSeq(true2,l=200)
#'   dm <- phangorn::dist.ml(aln)
#'   ape::bionj(dm)
#' })
#' 
#' holmesTest(list(set1,set2),phangorn::RF.dist)
#' holmesPlot(list(set1,set2),phangorn::RF.dist)
#' }
holmesTest <- function(x,
                       dist.fn=NULL,
                       labels=NULL,
                       B="automatic",
                       R="automatic",
                       returnNullDistribution=FALSE,
                       bootstrap.se.tolerance=0.005,
                       ...) {
  
  # Check inputs before doing costly things
  if ( (!is.numeric(B)) && (!grepl("auto",B)) ) {
    stop("Invalid input for number of bootstrap replicates B.")
  }

  if ( (!is.numeric(R)) && (!grepl("auto",R)) ) {
    stop("Invalid input for number of replicate non-unique MST trees R.")
  }
  
  tmp <- prepForHolmes(x=x,dist.fn=dist.fn,labels=labels)
  x <- tmp$x
  labels <- tmp$labels
  
  # recover()
  
  all_dists <- x[upper.tri(x)]
  ties_present <- length(unique(all_dists)) < length(all_dists)
  
  res <- NULL
  sparse_spanning_trees <- list()
  if (ties_present) {
    warning("There are non-unique distances between trees. MST may not be unique, averaging over R = ",R," MSTs. Note that p-value may change on re-run.")
    sparse_spanning_trees <- getSparseSpanningTreeListHolmes(x,R,labels)
    R <- length(sparse_spanning_trees)
  } else {
    sparse_spanning_trees[[1]] <- getSparseSpanningTree(x,shuffle=FALSE)
  }
  n_spanning_trees <- length(sparse_spanning_trees)
  
  S_0 <- holmesTestStat(sparse_spanning_trees,labels)
  
  p_se <- NA
  S_star <- numeric()
  if ( is.numeric(B) ) {
    S_star <- sapply(1:B,function(b){
      holmesTestStat(sparse_spanning_trees,sample(labels))
    })
  } else {
    
    S_star <- numeric(1000)
    S_star[1:100] <- sapply(1:100,function(b){
      holmesTestStat(sparse_spanning_trees,sample(labels))
    })
    
    idx <- 100
    se <- sd(S_star[1:100] > S_0)/sqrt(100)
    
    while (se > bootstrap.se.tolerance) {
      idx <- idx + 1
      if ( idx > length(S_star) ) {
        S_star <- c(S_star,vector("numeric",1000))
      }
      S_star[idx] <- holmesTestStat(sparse_spanning_trees,sample(labels))
      se <- sd(S_star[1:idx] > S_0)/sqrt(idx)
    }
    B <- idx
    S_star <- S_star[1:idx]
    p_se <- sd(S_star > S_0)/sqrt(B)
  }
  
  p <- sum(S_star > S_0)/B
  res <- list(
    p.value = p,
    S_0 = S_0,
    R = R,
    B = B,
    p.value.standard.error = p_se
  )
  if (returnNullDistribution) {
    res$nullDistribution = S_star
  }
  
  return(res)
}

#' @describeIn holmesTest
#' Uses classical multidimensional scaling to plot trees in a 2D plane and shows edges connecting the minimum spanning tree.
#' Optionally reports p-value.
#' @export
holmesPlot <- function(x,dist.fn=NULL,
                       labels=NULL,
                       B=NULL,
                       setColors=1+(1:min(c(length(x),length(unique(labels))))),
                       betweenColor=1) {
  tmp <- prepForHolmes(x=x,dist.fn=dist.fn,labels=labels)
  x <- tmp$x
  labels <- as.integer(as.factor(tmp$labels))
  ntrees <- tmp$ntrees
  
  p_value <- NULL
  if ( is.numeric(B) && B > 0 ) {
    p_value <- holmesTest(x=x, dist.fn=dist.fn, labels=tmp$labels, B=B)$p.value
  }
  
  sparse_spanning_tree <- getSparseSpanningTree(x,shuffle=TRUE)
  
  # recover()
  
  y <- stats::cmdscale(x,k=2)
  
  cols <- setColors[labels]
  
  plot(NULL,NULL,xlim=range(y[,1]),ylim=range(y[,2]),pch=16,col=cols,xlab="x1",ylab="x2")
  for (idx in 1:dim(sparse_spanning_tree)[1]) {
    i <- sparse_spanning_tree[idx,1]
    j <- sparse_spanning_tree[idx,2]
    
    # So as to only plot each line once
    if (i < j) {
      lc <- betweenColor
      if (labels[i] == labels[j]) {
        lc <- setColors[labels[i]]
      }
      lines(x=y[c(i,j),1],y=y[c(i,j),2],col=lc)
    }
  }
  points(y,pch=16,col=cols)
  if ( !is.null(p_value) ) {
    title(paste0("Permutation test p-value: ",round(p_value,3)))
  }
  
  return(NULL)
  
}

#' Helper function used to allow multiple input options for holmesPlot and holmesTest
#' 
#' @param x Either a list of multiPhylo objects or a distance matrix for distances between trees
#' @param dist.fn Function for computing tree to tree distances, must be capable of being called on a multiPhylo.
#' @param labels If x is a distance matrix, this argument assigns each column and row to a bootstrap replicate
#' @return Distance matrix
#' @keywords internal
prepForHolmes <- function(x,dist.fn,labels=NULL) {
  ntrees <- sapply(unique(labels),function(x){sum(labels == x)})
  if ( "dist" %in% class(x) ) {
    x <- as.matrix(x)
  } else if ( "list" %in% class(x) ) {
    if (!("function" %in% class(dist.fn))) {
      stop("If x is a list of trees, dist.fn must be a function for computing the distances between them.")
    }
    ntrees <- lengths(x)
    labels <- NULL
    for (i in 1:length(x)) {
      if ( !("multiPhylo" %in% class(x[[i]])) ) {
        stop("Invalid input.")
      }
      labels <- c(labels,rep(i,length(x[[i]])))
    }
    trees <- do.call(c,x)
    x <- as.matrix(dist.fn(trees))
  } else if ( !("matrix" %in% class(x)) ) {
    stop("Invalid input format, must be trees or distance matrix between trees.")
  }
  
  return(list(x=x, labels=labels, ntrees=ntrees))
}

#' Computes test statistic for Holme's tree test from a pre-computed sparse spanning tree and a set of labels
#' 
#' Allows a list of MSTs, for situations where the MST is not unique, in which case it averages the test statistic over these.
#' 
#' @param sparseSpanningTree A list of spanning trees (each provided as a list of edges)
#' @param labels Which tree goes in which set.
#' @return The test statistic (number of edges connecting trees within a set)
#' @keywords internal
holmesTestStat <- function(sparseSpanningTrees,labels) {
  
  if (class(sparseSpanningTrees) == "matrix") {
    sparseSpanningTrees <- list(sparseSpanningTrees)
  }
  
  stats <- unlist(lapply(sparseSpanningTrees,function(spanning_tree){
    # we've double-counted every edge by not restricting the count to i > j or j > i
    sum(labels[spanning_tree[,1]] == labels[spanning_tree[,2]])/2
  }))
  
  return(mean(stats))
}