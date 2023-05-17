# library(energy)
# 
# uniformlySubsampled.dCov.sq <- function(x,y,dist.fn,n.per.row,fast=TRUE) {
#   n <- length(x)
#   if (length(y) != n) {
#     stop("X and Y are not same length")
#   }
#   
#   # recover()
#   
#   # Sample n.per.row distances in each row, for a total of n.per.row * n
#   col.idx <- lapply(1:n,function(i){
#     sample.int(n,n.per.row)
#   })
#   
#   x <- lapply(1:n,function(i){
#     getSubsampledDistanceMatrix(x,dist.fn,cbind(rep(i,n.per.row),col.idx[[i]]))
#   })
#   y <- lapply(1:n,function(i){
#     getSubsampledDistanceMatrix(y,dist.fn,cbind(rep(i,n.per.row),col.idx[[i]]))
#   })
#   
#   # #TODO this misses samples which have the same _column_ and thus should contribute
#   x_row_means <- x_col_means <- NULL
#   y_row_means <- y_col_means <- NULL
#   if (fast) {
#     x_row_means <- x_col_means <- unlist(lapply(x,mean))
#     y_row_means <- y_col_means <- unlist(lapply(y,mean))
#   } else {
#     x_row_sums <- numeric(n)
#     x_row_counts <- rep(n.per.row,n)
#     y_row_sums <- numeric(n)
#     y_row_counts <- rep(n.per.row,n)
#     for (i in 1:n) {
#       x_row_sums[i] <- x_row_sums[i] + sum(x[[i]])
#       y_row_sums[i] <- y_row_sums[i] + sum(y[[i]])
#       for (idx in 1:n.per.row) {
#         j <- col.idx[[i]][idx]
#         if (j != i) {
#           x_row_sums[j] <- x_row_sums[j] + x[[i]][idx]
#           y_row_sums[j] <- y_row_sums[j] + y[[i]][idx]
#           x_row_counts[j] <- x_row_counts[j] + 1
#           y_row_counts[j] <- y_row_counts[j] + 1
#         }
#       }
#     }
#     x_row_means <- x_col_means <- x_row_sums/x_row_counts
#     y_row_means <- y_col_means <- y_row_sums/y_row_counts
#   }
#   
#   x_grand_mean <- mean(x_row_means)
#   y_grand_mean <- mean(y_row_means)
#   
#   summand <- sapply(1:n,function(i){
#     A <- x[[i]] - x_row_means[i] - x_col_means[col.idx[[i]]] + x_grand_mean
#     B <- y[[i]] - y_row_means[i] - y_col_means[col.idx[[i]]] + y_grand_mean
#     n * mean(A*B)
#   })
#   
#   return(sum(summand)/(n^2))
# 
# }
# 
# diagonallySubsampled.dCov.sq <- function(x,y,dist.fn,n.per.diag) {
#   n <- length(x)
#   if (length(y) != n) {
#     stop("X and Y are not same length")
#   }
#   
#   # recover()
#   
#   x_subsamp <- computeDiagonallySubsampledDistances(x,dist.fn,n.per.diag,FALSE)
#   y_subsamp <- computeDiagonallySubsampledDistances(y,dist.fn,n.per.diag,FALSE)
#   
#   x_diag_means <- sapply(1:(n-1),function(i){
#     mean(x_subsamp[[i]][,3])
#   })
#   y_diag_means <- sapply(1:(n-1),function(i){
#     mean(y_subsamp[[i]][,3])
#   })
#   
#   
#   # We are implicitly here handling the fact that the diagonal has all 0s, and it contributes a constant weight of 1
#   wt <- rep(1,n-1)
#   x_row_means <- numeric(n)
#   y_row_means <- numeric(n)
#   for (i in 1:ceiling(n/2)) {
#     row_first <- i
#     row_last <- n - i + 1
#     
#     x_row_means[row_first] <- x_row_means[row_last] <- 1/n * (sum(wt * x_diag_means))
#     y_row_means[row_first] <- y_row_means[row_last] <- 1/n * (sum(wt * y_diag_means))
#     wt[i] <- wt[i] + 1
#     wt[n - i] <- 0
#   }
#   
#   x_col_means <- x_row_means
#   y_col_means <- y_row_means
#   
#   x_grand_mean <- mean(x_col_means)
#   y_grand_mean <- mean(y_col_means)
#   
#   summand <- sapply(1:(n-1),function(i){
#     A_diag <- x_subsamp[[i]][,3] - x_row_means[x_subsamp[[i]][,1]] - x_col_means[x_subsamp[[i]][,2]] + x_grand_mean
#     B_diag <- y_subsamp[[i]][,3] - y_row_means[y_subsamp[[i]][,1]] - y_col_means[y_subsamp[[i]][,2]] + y_grand_mean
#     return((n - i) * mean(A_diag * B_diag))
#   })
#   
#   return(sum(summand)/(n^2))
# }
# 
# 
# # 
# # # x, y are distance matrices for samples of x, y
# # randsampled.dcov.sq <- function(x,y,n.per.diag=100) {
# #   # recover()
# # 
# #   n <- dim(x)[1]
# #   if (dim(y)[1] != n) {
# #     stop("X and Y are not same size")
# #   }
# # 
# #   matrix.idx <- sample.int(n^2,n.per.diag*n)
# #   row.idx <- sort(rep(1:n,n))[matrix.idx]
# #   col.idx <- (rep(1:n,n))[matrix.idx]
# # 
# #   x_row_means <- x_col_means <- sapply(1:n,function(i){
# #     mean(x[row.idx == i | col.idx == i])
# #   })
# # 
# #   y_row_means <- y_col_means <- sapply(1:n,function(i){
# #     mean(y[row.idx == i | col.idx == i])
# #   })
# # 
# #   x_grand_mean <- sum(sapply(1:length(matrix.idx),function(i){
# #     x[row.idx[i],col.idx[i]]
# #   }))/length(matrix.idx)
# # 
# #   y_grand_mean <- sum(sapply(1:length(matrix.idx),function(i){
# #     y[row.idx[i],col.idx[i]]
# #   }))/length(matrix.idx)
# # 
# #   res <- sum(sapply(1:length(matrix.idx),function(i){
# #     row <- row.idx[i]
# #     col <- col.idx[i]
# #     (x[row,col] - x_row_means[row] - x_col_means[col] + x_grand_mean) *
# #       (y[row,col] - y_row_means[row] - y_col_means[col] + y_grand_mean)
# #   }))/length(matrix.idx)/2
# # 
# #   return(res)
# # }
# 
# 
# 
# eucliddist <- function(x,y) {sqrt(sum((x - y)^2))}
# 
# 
# res <- sapply(1:100,function(i){
#   phi <- runif(1)
#   mcmc <- simAR1(1001,phi,1)
#   
#   x <- mcmc[1:1000]
#   y <- mcmc[2:1001]
#   
#   res <- c(
#     energy::dcov(dist(x),dist(y))^2,
#     diagonallySubsampled.dCov.sq(x,y,eucliddist,200),
#     uniformlySubsampled.dCov.sq(x,y,eucliddist,200,fast=TRUE),
#     uniformlySubsampled.dCov.sq(x,y,eucliddist,200,fast=FALSE)
#   )
#   names(res) <- c("full","diagonal","rowwise.fast","rowwise.slow")
#   return(res)
# })
# 
# plot(res[1,],res[3,],col="blue",log="xy")
# points(res[1,],res[2,],col="green")
# points(res[1,],res[4,],col="blue",pch=2)
# abline(a=0,b=1)
# 
# 
# mean(res[2,] - res[1,])
# mean(res[3,] - res[1,])
# mean(res[4,] - res[1,])
# 
# mean((res[2,] - res[1,])/res[1,])
# mean((res[3,] - res[1,])/res[1,])
# mean((res[4,] - res[1,])/res[1,])
# 
# mean(abs(res[2,] - res[1,]))
# mean(abs(res[3,] - res[1,]))
# mean(abs(res[4,] - res[1,]))
# 
# mean(abs(res[2,] - res[1,])/res[1,])
# mean(abs(res[3,] - res[1,])/res[1,])
# mean(abs(res[4,] - res[1,])/res[1,])
# 
# 
# 
# 
# 
# plot(res[1,],(res[2,] - res[1,])/res[1,],col="green",log="")
# points(res[1,],(res[3,] - res[1,])/res[1,],col="blue")
# points(res[1,],(res[4,] - res[1,])/res[1,],col="blue",pch=2)
# abline(h=0)
