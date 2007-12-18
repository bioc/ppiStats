estimatePPIErrorRates <- function(matList, GS){

    pfnList <- sapply(matList, function(x) {
        dnames1 <- dimnames(x[[1]])
        dnames2 <- dimnames(x[[2]])
        if(dnames1 != dnames2) {stop("The data matrix and tested matrix do not have the
                                     same names")}
        dimnamesGS <- dimnames(GS)
        ##first we need to be in the same universe to get anything of interest
        ##so we restrict the Gold Standard and the data/tested matrices to the
        ##common proteins 
        restrictedGS <- GS[intersect(dimnamesGS[[1]], dnames1[[1]]),
                           intersect(dimnamesGS[[2]], dnames1[[2]])]
        restrictedTested <- x[[2]][intersect(dimnamesGS[[1]], dnames1[[1]]),
                           intersect(dimnamesGS[[2]], dnames1[[2]])]

        restrictedData <- x[[1]][intersect(dimnamesGS[[1]], dnames1[[1]]),
                           intersect(dimnamesGS[[2]], dnames1[[2]])]

        ##in order to estimate the false negative rate, we need to know the number
        ##of postive interactions in the GS which been tested in the experiment but
        ##not found..this is computed and recorded in testedInGSNotFound. then we
        ##quotient this by all the positive interactions in the GS which have been
        ##tested in the experiment
        testedInGSNotFound <- sum(!restrictedData * restrictedTested * restrictedGS)
        testedInGS <- sum(restrictedTested * restrictedGS)
        
        pfnEst <- testedInGSNotFound/testedInGS
        pfnEst
    })

    ##x123 computes the statistics needed to estimate the pfp rate given the pfn rate
    x123 <- lapply(matList, getx123)

    
    pfpList <- mapply(x123,pfnList,getFPrateFromFNrate)
    unlist(pfpList)

    errorMat <- cbind(pfpList, pfnList)
    colnames(errorMat) <- c("pFP","pFN")
    rownames(errorMat) <- names(matList)
    return(errorMat)
    
}

##      input: a list of two matrices,pos and test, in ppMatrix format.
##      output:
##          x1: the number of reciprocally observed edges.
##          x2: the number of reciprocally unobserved edges.
##          x3: the number of unreciprocally observed edges.

 getx123 = function(ppMatrix) {
      x1<-sum(ppMatrix[["test"]]==2 & ppMatrix[["pos"]]==2)
      x2<-sum(ppMatrix[["test"]]==2 & ppMatrix[["pos"]]==0)
      x3<-sum(ppMatrix[["test"]]==2 & ppMatrix[["pos"]]==1)
      c(x1=x1,x2=x2,x3=x3)
  }


## this function calculate FPrate from FNrate for an experiment.
## input: x1,x2,x3 are from the output of getx123()
##        FNrate is the FNrate for that experiment.
## output: the FPrate.

 getFPrateFromFNrate <- function(y, FNrate) {
     x1 <- y[[1]]
     x2 <- y[[2]]
     x3 <- y[[3]]
     tot <- sum(y) 
      u<-x1-tot*(1-FNrate)^2
      v<-x2-tot*(FNrate)^2
      a<-u-v
      b<-(-2)*u
      c<-u-u*FNrate^2+v*(1-FNrate)^2
      d<-b^2-4*a*c
      if(d>=0)
        {
            FP1<-(-b+sqrt(d))/(2*a)
            FP2<-(-b-sqrt(d))/(2*a)
        } else{
          FP1<-NA
          FP2<-NA
      }
      if(!is.na(FP1) && FP1!=(1-FNrate)){ return(FP1)}
         else{ return(FP2)}
  }
