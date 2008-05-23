 

estimateCCMErrorRates <- function(m,GS,filterSystematic=TRUE,
                                  obsPropThresh=1,SystematicpThresh=.01){

  #options(error=recover)
  ##########################################################
  ##                                                      ##
  ##  Create two matrices mFilt and mVBPFilt.             ##
  ##  mFilt has rows for VBPs and columns for VBPs+VPOs.  ##
  ##  mVBPFilt is square with rows and columns for VBPs.  ##
  ##  Systemetic baits removed if indicated (default).    ##
  ##                                                      ##
  ##########################################################

  stopifnot(all(m %in% 0:1))
  
  allBaits <- rownames(m)
  allPrey <- setdiff(colnames(m),allBaits)
  ##TC writes that the allPrey seems to be prey Only and not all prey perse
  ##allPrey <- colnames(m)

  vb <- names(rowSums(m)>0)
  vp <- names(colSums(m)>0)
  vbp <- intersect(vb,vp)
  
  #VBPs <- allBaits[(rowSums(m)>0) & 
  #                 (colSums(m[,allBaits])>0)]
  ##TC writes that the above code fails ... subset out of bounds for
  ##m[,allBaits] since a bait need not be a prey
  VBPs <- vbp

  VPOs <- allPrey[(colSums(m[,allPrey])>0)]
  mVBP <- m[VBPs,c(VBPs,VPOs)]

  
  ## identify baits subject to systematic error
  

  if(filterSystematic){
    mb <- mVBP[VBPs,VBPs]
    mbp <- assessSymmetry(mb)$p
    lose <- which(mbp<SystematicpThresh)    
    goodBaits <- setdiff(VBPs,names(lose))
  }  else goodBaits <- VBPs
  
  
  ##remove baits prone to systematic error


  mFilt <- m[goodBaits,c(goodBaits,VPOs)]
  mVBPFilt <- mFilt[,goodBaits]

  
  ##########################################################
  ##                                                      ##
  ##  Now, look at the gold standard complexes.           ##
  ##  Keep only those for which the proportion specified  ## 
  ##  by obsPropThresh of the constituent members         ##
  ##  are found to be either non-systematic VBPs or VPOs. ##
  ##  NB we don't use VBs for the filtering criteria.     ##
  ##  Put results in GSviable and incidence matrix mGS.   ##                    
  ##                                                      ##
  ##########################################################

  keep1 <- apply(GS,MARGIN=2,FUN=obsProp,y=colnames(mFilt),obsPropThresh=obsPropThresh)


  ##make sure edges are stochastically distributed
  ##within a complex, e.g. we don't want all missing edges
  ##pointing at just one protein  


  keep <- rep(FALSE,length(keep1))
  for (i in 1:length(keep1)){
    if(keep1[i]){
      prots <- names(which(GS[,i]==1))
      keep[i] = checkEdgePattern(prots=prots,m=mFilt,pThresh=SystematicpThresh)
    }
  } 
  
   
  ##make GSviable and mGS
    
  GSviable = GS[,keep]
  mGS <- as.matrix(1*(GSviable %*% t(GSviable) >0))
  diag(mGS) <- 0



  
  ##########################################################
  ##                                                      ##
  ##  Restrict attention to doubly tested edges           ##
  ##  between VBPs.  Calculate numbers of                 ## 
  ##  reciprocated and unreciprocated observations on     ##
  ##  these 'true' edges.                                 ##
  ##                                                      ##
  ##########################################################

    
  #find VBPs in both data sets
  commonNodes <- intersect(rownames(mVBPFilt),rownames(mGS))
  
  mgs <- mGS[commonNodes,commonNodes]
  ms <- mVBPFilt[commonNodes,commonNodes]
  
  #find subset of edges in ms that are present in mgs
  #these are observed 'true positives'
  msTP <- ms * mgs
  
  
  msTPd <- getDegrees(msTP)
  
  mgsd <- getDegrees(mgs)[,"nr"]
  
  totalTPR <- sum(msTPd[,"nr"])/2
  totalTPU <- sum(msTPd[,"ni"] + msTPd[,"no"])/2	  	      
  totalTP <- sum(mgsd)/2
  

  ##########################################################
  ##                                                      ##
  ##  Calculate globalpTP and its standard error.         ##
  ##                                                      ##
  ##########################################################

  globalpTP <- estimatepTP(n=totalTP,r=totalTPR,
                           u=totalTPU)  
  globalpTPse <- sqrt(globalpTP*(1-globalpTP)/(2*totalTP))
  pTP95U = globalpTP+1.96*globalpTPse
  pTP95L = globalpTP-1.96*globalpTPse

  pTP95CI = c(pTP95L,pTP95U)
  names(pTP95CI) = c("95CIlb","95CIub")
  #print(globalpTP)
 
  ##########################################################
  ##                                                      ##  
  ##  Find global pFP estimate using the method of        ##
  ##  moments approach.                                   ##
  ##                                                      ##  
  ##########################################################

  degObs <- getDegrees(mVBPFilt) 
  totalObsR <- sum(degObs[,"nr"])/2  
  totalObsU <- sum(degObs[,"ni"]+degObs[,"no"])/2
  ntot = dim(mVBPFilt)[1]
  
  #find point estimate and CI for pFP
  pFPans = estimatepFP(pTPest=globalpTP,totalObsR=totalObsR,totalObsU=totalObsU,ntot=ntot)
  globalpFP = pFPans$pFPest
  probPairs = pFPans$probPairs

  pFP95U =  
estimatepFP(pTPest=pTP95U,totalObsR=totalObsR,totalObsU=totalObsU,ntot=ntot)$pFPest
  pFP95L =  
estimatepFP(pTPest=pTP95L,totalObsR=totalObsR,totalObsU=totalObsU,ntot=ntot)$pFPest
  
  pFP95CI = c(pFP95L,pFP95U)
  names(pFP95CI) = c("95CIlb","95CIub")


  ##########################################################
  ##                                                      ##  
  ##  Report summary statistics of portion of gold        ##
  ##  standard used.                                      ##
  ##                                                      ##  
  ##########################################################

  nEligComplexes <- dim(GSviable)[2]
  nEligBaits <- sum(rowSums(mgs)>0)
  nEligEdges <- sum(mgs)/2
  nBaitsInComplexes <- colSums(GSviable[commonNodes,])
  complexSizes <- colSums(GSviable)



  
  ans <- list(globalpTP=globalpTP,
       globalpTPse=globalpTPse,
       globalpFP=globalpFP,
       pTP95CI = pTP95CI,
       pFP95CI = pFP95CI,
       probPairs = probPairs,
       nEligComplexes=nEligComplexes,
       nEligBaits=nEligBaits,
       nEligEdges=nEligEdges,
       nBaitsInComplexes=nBaitsInComplexes,
       complexSizes=complexSizes)
  
  return(ans)
}






#####################################################
## other functions used in estimateCCMEErrorRates  ##
#####################################################

#filter function to find complexes meeting tresholding criteria
obsProp <- function(x,y, obsPropThresh=1){
  xn <- names(which(x==1))
  ans <-  length(intersect(xn,y))/sum(x) >= obsPropThresh
  ans
}

#get reciprocated and unreciprocated degrees
getDegrees = function(m){
  stopifnot(all(m %in% 0:1))
  m = m>0
  cbind(nr=rowSums(m & t(m)),no=rowSums(m & (!t(m))),
        ni=rowSums((!m) & t(m)))
}


#MLE estimate of pTP                                        
estimatepTP <- function(n,r,u){
  
  stopifnot(n>0)
  pTPhat <- (2*r+u)/(2*n)
  pTPhat

}


#MOM estimate for pFP corresponding to given pTP
estimatepFP = function(pTPest,totalObsR,totalObsU,ntot){
	
 nintEst = totalObsR/(pTPest^2)
 nint <- (max(floor(nintEst)-500,100)):(min(floor(nintEst)+500,ntot*(ntot-1)/2))

 MOMrange <-
    estErrProbMethodOfMoments(nint=nint,nrec=totalObsR,nunr=totalObsU,ntot=ntot)

 pFPs = c(MOMrange[,"pfp1"],MOMrange[,"pfp2"])
 pFNs = c(MOMrange[,"pfn1"],MOMrange[,"pfn2"])
 probPairs = cbind(pFPs,pFNs)
 probPairs = probPairs[order(pFPs),]
 probPairs = probPairs[probPairs[,"pFPs"]>0 & probPairs[,"pFPs"]<1,]
 probPairs = probPairs[probPairs[,"pFNs"]>0 & probPairs[,"pFNs"]<1,]

 closestpTP <-
    which.min(abs((1-pTPest)-probPairs[,"pFNs"]))
 pFPest = probPairs[,"pFPs"][closestpTP]	    

 ans = list(pFPest=pFPest,probPairs=probPairs)
 ans
}



#check edge pattern within candidate complex to ensure random 
#distribution of missing edges
checkEdgePattern = function(prots,m,pThresh){
      protsub <- intersect(prots,colnames(m))
      baitsub <- intersect(prots,rownames(m))
      preysub <- setdiff(protsub,baitsub)

      tempMat <- m[baitsub,c(baitsub,preysub),drop=FALSE]
      diag(tempMat) <- 0
      
      ##check to make sure edges are stochastically distributed with a complex
      if(dim(tempMat)[1]>1 & dim(tempMat)[2]<15){
        a1 <- rowSums(tempMat)
        b1 <- rowSums(1-tempMat)
        b1[baitsub] <- b1[baitsub]-1
        a2 <- colSums(tempMat)
        b2 <- colSums(1-tempMat)
        b2[baitsub] <- b2[baitsub]-1
        p1 <- fisher.test(rbind(a1,b1),workspace=2e7)$p.value
        p2 <- fisher.test(rbind(a2,b2),workspace=2e7)$p.value
        ans <- p1>0.01 & p2>0.01
      } else ans <- TRUE
      ans
    }
