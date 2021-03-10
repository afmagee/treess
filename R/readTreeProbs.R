#' Reads a MrBayes .trprobs file
#' 
#' Reads in trees and their associated posterior probabilities, returns a vector of (re-normalized) probabilities and a multiPhylo object.
#' 
#' @param filepath Filepath to MrBayes .trprobs file
#' @param CI.width All trees within the posterior CI.width*100 percentile will be returned (1.0 for all trees).
#' @param type "MrBayes" for MrBayes-style .trprobs files, "simple" for 2-column style, see details.
#' @param sep If type = "simple", the separator between the two columns, passed to read.table.
#' @return A list. $trees contains all unique topologies, and $probs their posterior probabilities (re-normalized to sum to 1).
#' @details Allowed input is either a MrBayes output file or a simple 2-column table, with probabilities in the first column and trees in the second.
#' @export
readTreeProbs <- function(filepath,CI.width=1.0,type="MrBayes",sep="") {
  # recover()
  probs <- NULL
  trees <- NULL
  if ( tolower(type) == "mrbayes" ) {
    trees <- ape::read.nexus(filepath)
    txt <- scan(filepath,what=character(),sep="\n")
    prob_txt <- txt[grepl("[&W",txt,fixed=TRUE)]
    prob_txt <- do.call(rbind,strsplit(prob_txt,"[&W",fixed=TRUE))[,2]
    prob_txt <- do.call(rbind,strsplit(prob_txt,"]",fixed=TRUE))[,1]
    probs <- as.numeric(prob_txt)
  } else if ( tolower(type) == "simple" ) {
    tmp <- read.table(filepath,header=FALSE,stringsAsFactors=FALSE,sep=sep)
    trees <- lapply(tmp[[2]],function(tree_txt){
      if ( !grepl(";",tree_txt) ) {
        tree_txt <- paste0(tree_txt,";")
      }
      return(ape::read.tree(text=tree_txt))
    })
    class(trees) <- "multiPhylo"
    probs <- as.numeric(tmp[[1]])
  } else {
    stop("Invalid input for argument \"type\"")
  }
  
  if ( CI.width < sum(probs) ) {
    cumprobs <- cumsum(probs)
    last_index <- min(which(cumprobs >= CI.width))
    trees <- trees[1:last_index]
    probs <- probs[1:last_index]
  }
  probs <- probs/sum(probs)
  
  return(list(trees=trees,probs=probs))
}
