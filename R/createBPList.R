createBPList <- function(indexSet, baitsSystematic, preysSystematic){

    baits <- lapply(indexSet, function(x) {sapply(x, function(y) y[1])})
    preys <- lapply(indexSet, function(x) {sapply(x, function(y) y[2])})

    
    bSys <- lapply(baits, function(x) {sapply(x, function(y) intAct2Sys(y, baitsSystematic))})
    
    

    pSys <- lapply(preys, function(x) {sapply(x, function(y) intAct2Sys(y, preysSystematic))})
        
    result <- list()

    
    for(i in 1:length(bSys)){
        
        result[[i]] <- split(pSys[[i]], bSys[[i]])
        
    }

    names(result) <- names(indexSet)
    result
    
}
