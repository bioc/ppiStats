idSystematic <- function(bpMat, viable, bpGraph=FALSE, pThresh=0.01,
                         pLevels = c(1e-06, 1e-04, 0.01)){

    if(bpGraph)
      bpMat <- as(bpMat, "matrix")

    stochastic <- idStochastic(bpMat = bpMat, pThresh = pThresh,
                               pLevels = pLevels)

    systematic <- setdiff(viable, stochastic)
    systematic

    
    
}
