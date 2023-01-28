# getNearestNeighborList <- function(x) {
#   lapply(1:dim(x)[1],function(i){
#     min_d <- min(x[i,-i])
#     return(which(x[i,] == min_d))
#   })
# }
# 
# nearestNeighborTestStat <- function(neighbor.list,labels) {
#   per_tree_props <- sapply(1:length(neighbor.list),function(i){
#     sum(labels[neighbor.list[[i]]] == labels[i])/length(neighbor.list[[i]])
#   })
#   return(mean(per_tree_props))
# }
# 
# nearestNeighborTest <- function(x,
#                                 dist.fn=NULL,
#                                 labels=NULL,
#                                 B=1000,
#                                 returnNullDistribution=FALSE) {
#   
#   tmp <- treess:::prepForHolmes(x=x,dist.fn=dist.fn,labels=labels)
#   x <- tmp$x
#   labels <- tmp$labels
#   
#   # recover()
#   
#   # This will make things much faster for computing test statistics
#   neighbor_list <- getNearestNeighborList(x)
#   
#   T_obs <- nearestNeighborTestStat(neighbor_list,labels)
#   
#   T_star <- sapply(1:B,function(i){
#     nearestNeighborTestStat(neighbor_list,sample(labels))
#   })
#   
#   p <- sum(T_star > T_obs)/B
#   res <- list(
#     p.value = p,
#     T_obs = T_obs
#   )
#   if (returnNullDistribution) {
#     res$nullDistribution = T_star
#   }
#   
#   return(res)
#   
# }