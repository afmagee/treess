#' Reads a MrBayes .trprobs file
#'
#' @param filepath Filepath to MrBayes .trprobs file
#' @param CI.width All trees within the posterior CI.width*100 percentile will be returned (1.0 for all trees).
#' @return A list. $trees contains all unique topologies, and $probs their posterior probabilities
#' @export
readTreeProbs <- function(filepath,CI.width=1.0) {
  # recover()
  trees <- ape::read.nexus(filepath)
  txt <- scan(filepath,what=character(),sep="\n")
  prob_txt <- txt[grepl("[&W",txt,fixed=TRUE)]
  prob_txt <- do.call(rbind,strsplit(prob_txt,"[&W",fixed=TRUE))[,2]
  prob_txt <- do.call(rbind,strsplit(prob_txt,"]",fixed=TRUE))[,1]
  probs <- as.numeric(prob_txt)
  if ( CI.width < sum(probs) ) {
    cumprobs <- cumsum(probs)
    last_index <- min(which(cumprobs >= CI.width))
    trees <- trees[1:last_index]
    probs <- probs[1:last_index]
  }
  probs <- probs/sum(probs)
  return(list(trees=trees,probs=probs))
}