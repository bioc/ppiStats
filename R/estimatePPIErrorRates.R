estimatePPIErrorRates <- function(matList, GSPos=NULL, GSNeg=NULL){
  printStatement <- "There was no Positive nor Negative Gold Standard Given, so
                     the error rates could not be estimated."

  ##If there is neither a positive nor negative GS, then this function
  ##is irrelevant
  if(is.null(GSPos) && is.null(GSNeg)){
    return(printStatement)
  }

  pfnList = NULL
  pfpList = NULL

  ##if we are given a positive gold standard, we will use it
  ##to estimate the pfn rates for each experiment
  if(!(is.null(GSPos))){
    pfnList <- getFNRateFromGSPos(matList=matList, GSPos=GSPos)
  }

  ##if we are give a negative gold standard, we will use it
  ##to estimate the pfp rates for each experiment
  if(!(is.null(GSNeg))){
    pfpList <- getFPRateFromGSNeg(matList=matList, GSNeg=GSNeg)
  }

  ##if both GS's are given, we will use the pfp and pfn estimates
  ##derived above
  if(!(pfnList) && !(pfpList)){
    errorMat <- cbind(pfnList, pfpList)
    colnames(errorMat) <- c("pfn","pfp")
    rownames(errorMat) <- names(matList)
    return(errorMat)
  }

  ##if the positive GS is given but not the negative, then we use
  ##the methods of moments to estimate the pfp rate
  if(is.null(GSNeg) && !(is.null(GSPos))){
    errorMat <- useMethodOfMoments(matList, pfnList, estFP=TRUE)
    return(errorMat)
  }

  ##if the negative GS is given but not the negative, then we use
  ##the methods of moments to estimate the pfn rate
  if(!(is.null(GSNeg)) && is.null(GSPos)){
    errorMat <- useMethodOfMoments(matList, pfnList, estFP=FALSE)
    return(errorMat)
  }
}

getFNRateFromGSPos <- function(matList, GSPos){

  pfnList <- sapply(matList, function(x) {
      dnames1 <- lapply(dimnames(x[[1]]), sort)
      dnames2 <- lapply(dimnames(x[[2]]), sort)
      if(!identical(dnames1,dnames2))
        {stop("The data matrix and tested matrix do not have the
                 same dimension names")}
      dnGS = dimnames(GSPos)
      ##first we need to be in the same universe to get anything of interest
      ##so we restrict the Gold Standard and the data/tested matrices to the
      ##common proteins 
      restrictedGS <- GSPos[intersect(dnGS[[1]], dnames1[[1]]),
                         intersect(dnGS[[2]], dnames1[[2]])]
      restrictedTested <- x[[2]][intersect(dnGS[[1]], dnames1[[1]]),
                                 intersect(dnGS[[2]], dnames1[[2]])]
      
      restrictedData <- x[[1]][intersect(dnGS[[1]], dnames1[[1]]),
                               intersect(dnGS[[2]], dnames1[[2]])]
      
      ##in order to estimate the false negative rate, we need to know the number
      ##of postive interactions in the GS which been tested in the experiment but
      ##not found..this is computed and recorded in testedInGSNotFound. then we
      ##quotient this by all the positive interactions in the GS which have been
      ##tested in the experiment
      testedInGSNotFound <- sum((!restrictedData) * restrictedTested * restrictedGS)
                      
      testedInGS <- sum(restrictedTested * restrictedGS)
                      
      pfnEst <- testedInGSNotFound/testedInGS
      pfnEst
    })

  return(pfnList)

}

getFPRateFromGSNeg <- function(matList, GSNeg){

  pfpList <- sapply(matList, function(x) {
    dnames1 <- lapply(dimnames(x[[1]]), sort)
    dnames2 <- lapply(dimnames(x[[2]]), sort)
    if(!identical(dnames1,dnames2))
      {stop("The data matrix and tested matrix do not have the
                 same dimension names")}
    dnGS = dimnames(GSNeg)
    ##first we need to be in the same universe to get anything of interest
    ##so we restrict the Gold Standard and the data/tested matrices to the
    ##common proteins 
    restrictedGS <- GSNeg[intersect(dnGS[[1]], dnames1[[1]]),
                          intersect(dnGS[[2]], dnames1[[2]])]
    restrictedTested <- x[[2]][intersect(dnGS[[1]], dnames1[[1]]),
                               intersect(dnGS[[2]], dnames1[[2]])]
    restrictedData <- x[[1]][intersect(dnGS[[1]], dnames1[[1]]),
                             intersect(dnGS[[2]], dnames1[[2]])]
    ##in order to estimate the false positive rate, we need to know the number
    ##of negative interactions in the GS which been tested in the experiment but
    ##were found..this is computed and recorded in notInGSFound. then we
    ##quotient this by all the positive interactions in the GS which have been
    ##tested in the experiment
    notInGSFound <- restrictedData*restrictedTested*restrictedGS
    notInGSTested <- restrictedTested*restrictedGS
    
    pfpEst <- notInGSFound/notInGSTested
    pfpEst
  })

  return(pfpList)
  
}


useMethodOfMoments <- function(matList, errorList, estFP=TRUE){

  ##x123 computes the statistics needed to estimate the pfp rate given the pfn rate
  x123 <- lapply(matList, getx123)
  
  if(estFP){
    errorList2 <- mapply(getFPrateFromFNrate, x123, errorList)
    errorMat <- cbind(errorList2, errorList)
    colnames(errorMat) <- c("pFP","pFN")
  }
  else{
    errorList2 <- mapply(getFNrateFromFPrate, x123, errorList)
    errorMat <- cbind(errorList2, errorList)
    colnames(errorMat) <- c("pFN","pFP")
    errorMat <- errorMat[,2:1]
    
  }

  unlist(errorList2)

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
  ##TC needs to check this...
  ##the roles of FP and FN are interchangable by symmetry, so if
  ##FNrate is really FPrate, the function can also double as
  ##getFNrateFromFPrate
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

getFNrateFromFPrate <- function(y, FPrate){
  ##fix me...need to deconvolute Li's code
  ##symmetry does not work because the number of
  ##truly interacting protein pairs (n) is
  ##not the same as the number of non-interacting
  ##pairs

  
  return("not working yet")

}
