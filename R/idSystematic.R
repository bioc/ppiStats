idSystematic <- function(bpMat, viable, bpGraph=FALSE, pThresh=0.01,
                         pLevels = 1e-4, prob=0.5){

    if(bpGraph)
      bpMat <- as(bpMat, "matrix")

    stochastic <- idStochastic(bpMat = bpMat, pThresh = pThresh,
                               pLevels = pLevels)

    systematic <- setdiff(viable, stochastic)
    systematic

    
    
}
