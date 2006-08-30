genBPGraph <- function(bpMat, directed=TRUE, bp=TRUE){

  bpMat1 <- bpMat
  b <- rownames(bpMats)
  p <- colnames(bpMats)

  if(!bp){
    if(b != p){
      stop("The rownames and the colnames must be identical.")
    }
  }

  else{
    baits <- union(rownames(bpMat), colnames(bpMat))
    preys <- baits

    bpMat1 <- matrix(0, length(baits), length(preys))
    dimnames(bpMat1) <- list(baits, preys)
    bpMat1[b,p] <- bpMat
  }

  bpGraph <- as(bpMat1, "graphNEL")
  bpGraph
  
}







