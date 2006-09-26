bpMatrix <- function(y2h, symMat = TRUE, homodimer=FALSE, baitAsPrey = FALSE,
                     unWeighted = TRUE, onlyRecip = FALSE, baitsOnly=FALSE){

  if(!baitAsPrey && baitsOnly){
    stop("If baitsOnly is TRUE, then baitAsPrey must also be TRUE")
    ##This is because if we want a bait by bait matrix, the bait proteins
    ##must also be prey proteins (only valid if the experimental procedure
    ##is genome-wide)
  }

  if(symMat && baitsOnly){
    stop("If baitsOnly is TRUE, symMats must be set to FALSE")
    ##symMats will give a matrix where the union of the baits and preys
    ##will index both the rows and columns. This is necessary if we
    ##want to generate an object of class graph. 
  }

  
  
  baits <- unique(names(y2h))
  preys <- unique(unlist(y2h))
  
  if(symMat){

    allProt <- union(baits, preys)
    baits <- allProt
    preys <- allProt
    
  }
  
  if(baitAsPrey){
    preys <- union(baits, preys)
  }
  
  bpMat <- matrix(0, nrow <- length(baits), ncol <- length(preys))
  dimnames(bpMat) <- list(baits, preys)

  nonTrivialIndex <- which(sapply(y2h,length) != 0)
  btmp <- baits[nonTrivialIndex]
  
  if(length(btmp)!=0 && length(length(y2h[nonTrivialIndex])!=0)){
    
    preyCount <- lapply(y2h, table)
    
    preyCount <- preyCount[nonTrivialIndex]
    
    for(i in 1:length(preyCount)){
      bpMat[btmp[i], names(preyCount[[i]])] = as.vector(preyCount[[i]])
    }
    
  }


  if(onlyRecip){
      bpMat <- (bpMat * t(bpMat))
      ##This will give a symmetric matrix, i.e. only those interactions
      ##where reciprocity was seen is given by the matrix.
  }
  

  
  if(!homodimer){diag(bpMat)=0}
  if(unWeighted){
    mode(bpMat) = "logical"
    mode(bpMat) = "numeric"
  }


  
  if(baitsOnly){
    bpMat <- bpMat[,rownames(bpMat), drop=FALSE]
  }
    
  bpMat
  
}

