idStochastic <- function(bpMat, bpGraph=FALSE, pThresh=0.01,
                         pLevels = 1e-4, prob=0.5){

    if(bpGraph)
      bpMat <- as(bpMat, "matrix")

    f <- assessSymmetry(bpMat, pLevels = pLevels, prob=prob)
    sel <- (f$p>pThresh)
    return(names(sel)[sel])
    
}

  
