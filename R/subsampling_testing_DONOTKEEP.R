# set.seed(47)
# 
# library(treess)
# library(phangorn)
# 
# source("~/git_repos/tree_convergence_code/real_data_examples/src/utils.R")
# 
# #######
# # Get the logfiles
# #######
# 
# cop <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/cophyline.tar.gz")
# gep <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/gephyromantis.tar.gz")
# het <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/heterixalus.tar.gz")
# par <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/paroedura.tar.gz")
# phe <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/phelsuma.tar.gz")
# uro <- read.gz.nexus.trees("~/git_repos/tree_convergence_code/real_data_examples/data/uroplatus.tar.gz")
# 
# malagasy <- list(cop,gep,het,par,phe,uro)
# names(malagasy) <- c("Cophyline","Gephyromantis","Heterixalus","Paroedura","Phelsuma","Uroplatus")
# 
# # Logfiles are said to contain 4 runs, so presumably they are concatenated sequentially
# chain.idx <- list(
#   1:1001,
#   1002:2002,
#   2003:3003,
#   3004:4004
# )
# 
# malagasy.chains <- lapply(malagasy,function(herps){
#   lapply(chain.idx,function(idx){
#     herps[idx]
#   })
# })
# malagasy.chains <- unlist(malagasy.chains,recursive=F)
# 
# est.s0.lanfear <- function(dmat,max.approximateESS.timelag=100,alpha=0.05) {
#   dmat <- dmat^2
#   
#   N <- n <- dim(dmat)[1]
#   
#   # make sure we don't over shoot
#   max_off_diag <- max.approximateESS.timelag
#   if ( max_off_diag > dim(dmat)[1]-1 ) {
#     max_off_diag <- dim(dmat)[1]-1
#   }
#   
#   # get empirical curve of lag time t vs distance
#   # see also topological.autocorr
#   t <- 1:max_off_diag
#   d_t <- sapply(t,function(t_){
#     mean(dmat[row(dmat) == col(dmat) + t_])
#   })
#   
#   m <- .computeApproximateESSThreshold(d_t, t, alpha)
#   
#   return(m)
# }
# 
# # stolen from jumpDistanceEquivalentESS
# est.s0.full <- function(dmat,alpha=0.05,nsim=1000,bootstrap=FALSE,central_tendency=mean,min.nsamples=5) {
#   dmat <- dmat^2
#   
#   # recover()
#   
#   n <- dim(dmat)[1]
#   # permute matrices to get appropriate autocorrelation
#   off_diag <- row(dmat) == col(dmat)+1
#   null_dist <- sapply(1:nsim,function(i){
#     # The new order of the rows and columns
#     permute <- sample.int(n,replace=bootstrap)
#     # Do not the entire distance matrix just to take an off diagonal
#     # Instead, draw indices of new rows and columns, and extract the diagonal from those indices
#     central_tendency(sapply(1:(n-1),function(i){dmat[permute[i],permute[i+1]]}))
#   })
#   threshold <- quantile(null_dist,probs=alpha)
#   
#   # guarantees if we don't find a better thinning value we set ESS = 1
#   thin <- n 
#   # get distances at time lags until we hit the threshold
#   G_s <- rep(NA,n-min.nsamples+1)
#   G_s[1] <- central_tendency(dmat[row(dmat) == col(dmat)+1])
#   if ( G_s[1] > threshold ) {
#     # Check for ESS = n
#     thin <- 1
#   } else {
#     for (i in 2:(n-min.nsamples-1)) {
#       # distance at this time lag
#       g_s <- central_tendency(dmat[row(dmat) == col(dmat)+i])
#       G_s[i] <- max(g_s,G_s[i-1])
#       # early terimination to avoid unneeded computation
#       if (G_s[i] > threshold) {
#         thin <- i
#         break
#       }
#     }
#   }
#   return(thin)
# }
# 
# est.s0.subsamp <- function(dmat,alpha=0.05,nsim=1000,n.per.diag=200,
#                            bootstrap=FALSE,central_tendency=mean,min.nsamples=5) {
#   dmat <- dmat^2
#   
#   # recover()
#   
#   # guarantees if we don't find a better thinning value we set ESS = 1
#   n <- dim(dmat)[1]
#   # permute matrices to get appropriate autocorrelation
#   off_diag <- row(dmat) == col(dmat)+1
#   idx <- seq(1,n-1,length.out=n.per.diag)
#   null_dist <- sapply(1:nsim,function(i){
#     # The new order of the rows and columns
#     permute <- sample.int(n,replace=bootstrap)
#     # Do not the entire distance matrix just to take an off diagonal
#     # Instead, draw indices of new rows and columns, and extract the diagonal from those indices
#     central_tendency(sapply(1:n.per.diag,function(i){dmat[permute[idx[i]],permute[idx[i]+1]]}))
#   })
#   threshold <- quantile(null_dist,probs=alpha)
#   
#   n_compute_per_diag <- pmin(n.per.diag,(n-1):1)
#   
#   thin <- n 
#   # get distances at time lags until we hit the threshold
#   G_s <- rep(NA,n-min.nsamples+1)
#   for (i in 1:(n - min.nsamples + 1)) {
#     # distance at this time lag
#     n_on_diagonal <- n - i
#     d <- i + 1
#     idx <- round(seq(1,n_on_diagonal,length.out=n_compute_per_diag[i]))
#     ij <- dk2ij(rep(d,n_compute_per_diag[i]),idx)
#     g_s <- central_tendency(apply(ij,1,function(ij_){
#       dmat[ij_[1],ij_[2]]
#     }))
#     G_s[i] <- max(g_s,G_s[i-1])
#     # early terimination to avoid unneeded computation
#     if (G_s[i] > threshold) {
#       thin <- i
#       break
#     }
#   }
#   return(thin)
# }
#  
# est.s0.t.asymototic <- function(dmat,n.per.diag=100,alpha=0.05,nsim=1000,bootstrap=FALSE,central_tendency=mean,min.nsamples=5) {
#   # recover()
#   dmat <- dmat^2
#   
#   n <- dim(dmat)[1]
#   vals <- dmat[upper.tri(dmat)]
#   m <- mean(vals)
#   s <- sd(vals)
#   
#   threshold <- m + qt(alpha,df=n-1) * s/sqrt(n)
#   # thresholds <- m + qt(alpha,df=n-1) * s/sqrt(n) * (sqrt(n)/sqrt(n:1))
#   
#   # guarantees if we don't find a better thinning value we set ESS = 1
#   thin <- n 
#   # get distances at time lags until we hit the threshold
#   G_s <- rep(NA,n-min.nsamples+1)
#   G_s[1] <- central_tendency(dmat[row(dmat) == col(dmat)+1])
#   if ( G_s[1] > threshold ) {
#     # Check for ESS = n
#     thin <- 1
#   } else {
#     for (i in 2:(n-min.nsamples-1)) {
#       # distance at this time lag
#       g_s <- central_tendency(dmat[row(dmat) == col(dmat)+i])
#       G_s[i] <- max(g_s,G_s[i-1])
#       # early terimination to avoid unneeded computation
#       if (G_s[i] > threshold) {
#         thin <- i
#         break
#       }
#     }
#   }
#   return(thin)
# }
# 
# est.s0.t.asymototic.subsamp <- function(x,dist.fn,n.per.diag=100,alpha=0.05,nsim=1000,bootstrap=FALSE,central_tendency=mean,min.nsamples=5) {
#   n <- length(x)
#   
#   i_j_dist <- computeDiagonallySubsampledDistances(x,dist.fn,n.per.diag,FALSE)
#   
#   vals <- unlist(lapply(i_j_dist,function(ijd){ijd[,3]^2}))
#   m <- mean(vals)
#   s <- sd(vals)
#   
#   # threshold <- confint(lm(vals ~ 1),level=(1 - 2 * alpha))
#   threshold <- m + qt(alpha,df=n-1) * s/sqrt(n)
#   
#   # guarantees if we don't find a better thinning value we set ESS = 1
#   thin <- n 
#   # get distances at time lags until we hit the threshold
#   G_s <- rep(NA,n-min.nsamples+1)
#   G_s[1] <- central_tendency(i_j_dist[[1]][,3]^2)
#   if ( G_s[1] > threshold ) {
#     # Check for ESS = n
#     thin <- 1
#   } else {
#     for (i in 2:(n-min.nsamples-1)) {
#       # distance at this time lag
#       g_s <- central_tendency(i_j_dist[[i]][,3]^2)
#       G_s[i] <- max(g_s,G_s[i-1])
#       # early terimination to avoid unneeded computation
#       if (G_s[i] > threshold) {
#         thin <- i
#         break
#       }
#     }
#   }
#   return(thin)
# }
# 
# est.s0.t.subsamp <- function(x,dist.fn,n.per.diag=100,alpha=0.05,nsim=1000,bootstrap=FALSE,central_tendency=mean,min.nsamples=5) {
#   n <- length(x)
#   
#   # recover()
#   
#   null_dist <- sapply(1:1000,function(idx){
#     permute <- sample.int(n,replace=bootstrap)
#     dists <- getSubsampledDistanceMatrix(x,dist.fn,cbind(permute[1:(n.per.diag-1)],permute[2:n.per.diag]))
#     return(mean(dists^2))
#   })
#   
#   threshold <- quantile(null_dist,probs=alpha)
#   
#   i_j_dist <- computeDiagonallySubsampledDistances(x,dist.fn,n.per.diag,FALSE)
#   
#   # guarantees if we don't find a better thinning value we set ESS = 1
#   thin <- n 
#   # get distances at time lags until we hit the threshold
#   G_s <- rep(NA,n-min.nsamples+1)
#   G_s[1] <- central_tendency(i_j_dist[[1]][,3]^2)
#   if ( G_s[1] > threshold ) {
#     # Check for ESS = n
#     thin <- 1
#   } else {
#     for (i in 2:(n-min.nsamples-1)) {
#       # distance at this time lag
#       g_s <- central_tendency(i_j_dist[[i]][,3]^2)
#       G_s[i] <- max(g_s,G_s[i-1])
#       # early terimination to avoid unneeded computation
#       if (G_s[i] > threshold) {
#         thin <- i
#         break
#       }
#     }
#   }
#   return(thin)
# }
# 
# # spanningTreeRatio <- function(dmat) {
# #   sequential_tree_length <- sum(dmat[row(dmat) == col(dmat) + 1])
# #   mst <- getSparseSpanningTree(dmat,TRUE)
# #   mst_length <- sum(apply(mst,1,function(ij){
# #     dmat[ij[1],ij[2]]
# #   }))
# #   return(mst_length/sequential_tree_length)
# # }
# 
# # # For some reason this is never as sensitive to autocorrelation as Holmes' test is
# # est.s0.mst <- function(dmat,nsim=2000,alpha=0.05) {
# #   # recover()
# # 
# #   n <- dim(dmat)[1]
# # 
# #   mst <- getSparseSpanningTree(dmat,TRUE)
# # 
# #   summary_fn_index_dist <- function(idx,mst) {
# #     mean((idx[mst[,1]] - idx[mst[,2]])^2)
# #     # nsteps <- abs(idx[mst[,1]] - idx[mst[,2]])
# #     # sum(nsteps == 1)/length(idx)
# #   }
# # 
# #   obs <- summary_fn_index_dist(1:n,mst)
# #   nd <- sapply(1:nsim,function(i){
# #     summary_fn_index_dist(sample.int(n),mst)
# #   })
# # 
# #   p <- sum(nd > obs)/nsim
# # 
# #   if ( p >= alpha ) {
# #     return(1)
# #   } else {
# #     thin <- n
# #     for (i in 2:n) {
# #       samps <- seq(1,n,i)
# #       dm <- dmat[samps,samps]
# # 
# #       mst <- getSparseSpanningTree(dm,TRUE)
# # 
# #       obs <- summary_fn_index_dist(1:length(samps),mst)
# #       nd <- sapply(1:nsim,function(i){
# #         summary_fn_index_dist(sample.int(length(samps)),mst)
# #       })
# # 
# #       p <- sum(nd > obs)/nsim
# # 
# #       if ( p >= alpha ) {
# #         thin <- i
# #         break
# #       }
# # 
# #       i <- i + 1
# #     }
# #     return(thin)
# #   }
# # 
# # }
# 
# 
# est.s0.dcov.sequential <- function(dmat,nsim=100,alpha=0.05) {
#   # recover()
#   nsamp <- dim(dmat)[1]
# 
#   p <- 0
#   thin <- 0
#   while (p <= alpha) {
#     thin <- thin + 1
#     dm <- dmat[seq(1,nsamp,thin),seq(1,nsamp,thin)]
#     n <- dim(dm)[1]
#     p <- energy::dcov.test(dm[-1,-1],dm[-n,-n],R=nsim)$p.value
#   }
#   return(thin)
# }
# 
# est.s0.dcov <- function(dmat,nsim=100,alpha=0.05) {
#   # recover()
#   min.nsamples <- 5
#   n <- dim(dmat)[1]
#   
#   # permute matrices to get appropriate autocorrelation
#   null_dist <- sapply(1:nsim,function(i){
#     idx <- sample.int(n)
#     permute <- dmat[idx,idx]
#     energy::dcov(permute[1:(n-1),1:(n-1)],permute[2:n,2:n],index=1.0)
#   })
#   threshold <- quantile(null_dist,probs=1 - alpha)
#   
#   thin <- n-min.nsamples+1
#   dc <- rep(NA,n-min.nsamples+1)
#   dc[1] <- energy::dcov(dmat[1:(n-1),1:(n-1)],dmat[2:n,2:n],index=1.0)
#   for (i in 2:(n-min.nsamples+1)) {
#     dc[i] <- min(dc[i-1],energy::dcov(dmat[-c(1:i),-c(1:i)],dmat[-c((n-i+1):n),-c((n-i+1):n)],index=1.0))
#     if ( dc[i] < threshold ) {
#       thin <- i
#       break
#     }
#   }
#   return(thin)
# }
# 
# source("~/git_repos/tree_convergence_code/Normal_random_walk/src/simpleMVRW.R")
# 
# eucliddist <- function(x,y) {
#   sqrt(sum(x - y)^2)
# }
# 
# 
# x <- simAR1(1000,0.9,1)
# dm <- as.matrix(dist(cbind(x)))
# # x <- rnorm(1000)
# # dm <- as.matrix(dist(cbind(x)))
# # midx <- 1
# # dm <- as.matrix(RF.dist(malagasy.chains[[midx]],rooted=T))
# # x <- simpleMVRW(ngen=1000,prop.sd=0.95,dimension=1)
# # dm <- as.matrix(dist(x))
# coda::effectiveSize(x)
# 
# est.s0.full(dm,nsim=5e2)
# est.s0.subsamp(dm)
# est.s0.t.asymototic(dm)
# est.s0.t.asymototic.subsamp(x,eucliddist)
# # est.s0.t.randsamp(dm)
# est.s0.t.subsamp(x,eucliddist)
# # est.s0.t.subsamp(treeSplits(malagasy.chains[[midx]],TRUE),rfDistFromSplits)
# est.s0.lanfear(dm)
# # est.s0.mst(dm)
# est.s0.dcov.sequential(dm,100)
# est.s0.dcov(dm,100,alpha=0.05)
# est.s0.dcov(dm,100,alpha=0.5)
# 
# 
# 
# res <- sapply(1:length(malagasy.chains),function(idx){
#   dm <- as.matrix(RF.dist(malagasy.chains[[idx]],rooted=TRUE))
#   # thin <- sample.int(100,1)
#   # mcmc <- simpleMVRW(ngen=1000*thin,prop.sd=0.45,dimension=1)[seq(1,1000*thin,length.out=1000)]
#   # dm <- as.matrix(dist(mcmc))
#   res <- c(est.s0.full(dm,nsim=1e3),
#            est.s0.t.asymototic(dm),
#            est.s0.t.asymototic.subsamp(treeSplits(malagasy.chains[[idx]],TRUE),rfDistFromSplits,100),
#            est.s0.t.subsamp(treeSplits(malagasy.chains[[idx]],TRUE),rfDistFromSplits,100),
#            est.s0.mst(dm),
#            est.s0.lanfear(dm)
#   )
#   names(res) <- c("bootstrap",
#                   "asymptotic",
#                   "asymptotic.subsamp",
#                   "bootstrap.subsamp",
#                   "mst",
#                   "lanfear")
#   res
# })
# 
# # res
# 
# plot(res[1,],res[2,],col="blue",log="")
# points(res[1,],res[3,],col="green")
# points(res[1,],res[4,],col="purple")
# points(res[1,],res[5,],col="grey")
# points(res[1,],res[6,],col="red")
# abline(a=0,b=1)
# 
# abline(lm(res[2,] ~ res[1,]),col="blue")
# abline(lm(res[3,] ~ res[1,]),col="green")
# abline(lm(res[4,] ~ res[1,]),col="purple")
# abline(lm(res[5,] ~ res[1,]),col="grey")
# abline(lm(res[6,] ~ res[1,]),col="red")
# 
# 
# sqrt(mean((res[1,]-res[2,])^2))
# sqrt(mean((res[1,]-res[3,])^2))
# sqrt(mean((res[1,]-res[4,])^2))
# sqrt(mean((res[1,]-res[5,])^2))
# sqrt(mean((res[1,]-res[6,])^2))
# 
# 
# mean(abs((res[1,]-res[2,])))
# mean(abs((res[1,]-res[3,])))
# mean(abs((res[1,]-res[4,])))
# mean(abs((res[1,]-res[5,])))
# mean(abs((res[1,]-res[6,])))
# 
# 
# plot(res[2,],res[3,])
# abline(a=0,b=1)
