idStochastic <- function(bpMat, bpGraph=FALSE, pThresh=0.01,
                         pLevels = c(1e-06, 1e-04, 0.01)){

    if(bpGraph)
      bpMat <- as(bpMat, "matrix")

    f <- assessSymmetry(bpMat, pLevels = pLevels)
    sel <- (f$p>pThresh)
    return(names(sel)[sel])
    
}

  
