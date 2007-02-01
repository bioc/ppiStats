idStochastic <- function(bpMat, bpGraph=FALSE, pThresh=0.01,
                         pLevels = c(1e-06, 1e-04, 0.01), prob=0.5){

    if(bpGraph)
      bpMat <- as(bpMat, "matrix")

    f <- assessSymmetry(bpMat, pLevels = pLevels, prob=prob)
    sel <- (f$p>pThresh)
    return(names(sel)[sel])
    
}

  
