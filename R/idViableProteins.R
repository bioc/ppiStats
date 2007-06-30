idViableProteins <- function(bpGraph, homomer=TRUE){

  bpMat <- as(bpGraph, "matrix")
  if (!homomer){
    diag(bpMat) <- 0
  }
    #deg <- degree(bpGraph)
    #vb <- names(which(deg$outDegree > 0))
    #vp <- names(which(deg$inDegree > 0))
    
  vb <- names(which(rowSums(bpMat)>0))
  vp <- names(which(colSums(bpMat)>0))
  vbp <- names(which(rowSums(bpMat)>0 & colSums(bpMat)>0))
  vbp1 <- intersect(vb, vp)
  if(length(setdiff(vbp, vbp1))>0 || length(setdiff(vbp1, vbp))){
    stop("There is bug")
  }
  vProt <- list()
  vProt$VB <- vb
  vProt$VP <- vp
  vProt$VBP <- vbp
  vProt
}
