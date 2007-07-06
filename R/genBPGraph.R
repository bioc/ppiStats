genBPGraph <- function(bpMat, directed=TRUE, bp=TRUE){

  bpMat1 <- bpMat
  b <- rownames(bpMat)
  p <- colnames(bpMat)

  if(!bp){
    if(sum(b != p) != 0){
      stop("The rownames and the colnames must be identical.")
    }
  }

  else{
    baits <- union(rownames(bpMat), colnames(bpMat))
    preys <- baits

    bpMat1 <- matrix(0, length(baits), length(preys))
    dimnames(bpMat1) <- list(baits, preys)
    bpMat1[b,p] <- bpMat
    if(!directed) {
      bpMat1 <- bpMat1 + t(bpMat1)
      mode(bpMat1) <- "logical"
      mode(bpMat1) <- "numeric"
    }
  }

  bpGraph <- as(bpMat1, "graphNEL")
  
  bpGraph
  
}







