estimateCCMErrorRates <- function(m,GS,filterSystematic=TRUE,
         obsPropThresh=.7,SystematicpThresh=.01, datamat=NULL){

  ##there is apart of the code when you use datamat, but you never
  ##define it as a local variable. you have defined datamat as a global
  ##variable in your script...it is what was passed into this function
  ##as the parameter m
  if(is.null(datamat)){
    datamat <- m
  }

  stopifnot(all(m %in% 0:1))
  
  mOrig <- m
  allBaits <- rownames(mOrig)
  allPrey <- setdiff(colnames(mOrig),allBaits)
  
  ##get VBP graph
  VBPs <- allBaits[(rowSums(m)>0) & 
                   (colSums(m[,allBaits])>0)]
  VPOs <- allPrey[(colSums(m[,allPrey])>0)]
  
  m <- m[VBPs,c(VBPs,VPOs)]

  ##filter nodes subject to systematic error
  if(filterSystematic){
    mb <- m[VBPs,VBPs]
    mbp <- assessSymmetry(mb)$p
    lose <- which(mbp<SystematicpThresh)
    
    goodBaits <- setdiff(VBPs,names(lose))
  }
  else goodBaits <- VBPs
  
  
  ##remove baits prone to systematic error
  ##mFilt is the matrix with no systematic proteins but has VP
  ##m now has no systematic proteins nor any only VPs
  mFilt <- mOrig[goodBaits,c(goodBaits,VPOs)]
  m <- mFilt[,goodBaits]

  ##now, we look at the gold standard complexes and keep only those
  ##complexes for which a large proportion of the constitute members
  ##are found to be either non-systematic proteins or VPs...NB we
  ##don't use VBs for the filtering criteria
  keep1 <- apply(GS,MARGIN=2,FUN=obsProp,y=colnames(mFilt))
  nEligComplexes <- sum(keep1)
  
  
  keep <- rep(FALSE,length(keep1))
  for (i in 1:length(keep1)){
    if(keep1[i]){
      prots <- names(which(GS[,i]==1))
      protsub <- intersect(prots,c(goodBaits,VPOs))
      baitsub <- intersect(prots,goodBaits)
      preysub <- setdiff(protsub,baitsub)
      ##datamat should be passed into the function as well...i don't know
      ##where datamat came from.
      ##TC says (a day later)...it seems that datamat is the same as mOrig
      tempMat <- datamat[baitsub,c(baitsub,preysub),drop=FALSE]
		diag(tempMat) <- 0
      
      ##don't know what is going on here in this if statement
      if(dim(tempMat)[1]>1 & dim(tempMat)[2]<15){
        a1 <- rowSums(tempMat)
        b1 <- rowSums(1-tempMat)
        b1[baitsub] <- b1[baitsub]-1
        a2 <- colSums(tempMat)
        b2 <- colSums(1-tempMat)
        b2[baitsub] <- b2[baitsub]-1
        p1 <- fisher.test(rbind(a1,b1),workspace=2e7)$p.value
        p2 <- fisher.test(rbind(a2,b2),workspace=2e7)$p.value
        keep[i] <- p1>0.01 & p2>0.01
      } else keep[i] <- TRUE
    }
  } 
  
  
  complexSizes <- colSums(GS[,keep])
  
  
  mGS <- as.matrix(1*(GS[,keep] %*% t(GS[,keep]) >0))
  diag(mGS) <- 0

  complexSizes <- colSums(GS[,keep])
  
  
  mGS <- as.matrix(1*(GS[,keep] %*% t(GS[,keep]) >0))
  diag(mGS) <- 0
  
  
  
  ##tony asks...are you only finding the baits that was also
  ##in the GS...not all nodes...
                                        #find nodes in both data sets
  commonNodes <- intersect(rownames(m),rownames(mGS))
  mgs <- mGS[commonNodes,commonNodes]
  
  nEligBaits <- sum(rowSums(mgs)>0)
  nEligEdges <- sum(mgs)/2
  nBaitsInComplexes <- colSums(GS[commonNodes,keep])
  
  ##tony commented out this print statement...
                                        #print(summary(nBaitsInComplexes))
  
  ##tony asks...why use m and not mFilt?
  ms <- m[commonNodes,commonNodes]
  
                                        #find subset of edges in ms that are present in mgs
                                        #these are observed 'true positives'
  msTP <- ms * mgs
  
  msTPd <- getDegrees(msTP)
  mgsd <- getDegrees(mgs)[,"nr"]
  
                                        #find global pTP estimate
  
  totalTPR <- sum(msTPd[,"nr"])/2
  totalTPU <- sum(msTPd[,"ni"] + msTPd[,"no"])/2	  	      
  totalTP <- sum(mgsd)/2
  
  globalpTP <- estimatepTP(n=totalTP,r=totalTPR,
                           u=totalTPU)[["pTPhat"]]
  
                                        #do bootstrap sampling
  pTPBS <- rep(0,nBS)
  datavec <- c(rep("r",totalTPR),rep("u",totalTPU),
               rep("n",(totalTP-totalTPR-totalTPU)))
  
  for (bss in 1:nBS){
    samp <- sample(x=datavec,size=totalTP,replace=TRUE)
    sampR <- sum(samp=="r")
    sampU <- sum(samp=="u")
    pTPBS[bss] <- estimatepTP(n=totalTP,r=sampR,u=sampU)[["pTPhat"]] 
  }
  
  pTPBSse <- sqrt(sum((pTPBS-mean(pTPBS))^2)/(nBS-1))
  pTPBS95L <- quantile(pTPBS,probs=0.025)
  pTPBS95U <- quantile(pTPBS,probs=0.975)
  
  
  globalpTPse <- sqrt(globalpTP*(1-globalpTP)/(2*sum(mgsd)/2))
  
  
                                        #find global pFP estimate
  
  degObs <- getDegrees(m)
  
  totalObsR <- sum(degObs[,"nr"])/2
  
  totalObsU <- sum(degObs[,"ni"]+degObs[,"no"])/2
  
                                        #find point estimate for pFP
  nintEst <- totalObsR/(globalpTP^2)
                                        #print(nintEst)
  nint <- (max(floor(nintEst)-500,100)):(min(floor(nintEst)+500,dim(m)[1]*(dim(m)[1]-1)/2))
  MOMrange <-
    estErrProbMethodOfMoments(nint=nint,nrec=totalObsR,nunr=totalObsU,ntot=dim(m)[1])
  closestpTP <-
    which.min(abs((1-globalpTP)-c(MOMrange[,3],MOMrange[,5])))
  globalpFP <- c(MOMrange[,2],MOMrange[,4])[closestpTP]
  
                                        #find upper bound for pFP
  
  pTP95U <- globalpTP+1.96*globalpTPse
  
  nintEst <- totalObsR/(pTP95U^2)
                                        #print(nintEst)
  nint <- (max(floor(nintEst)-500,100)):(min(floor(nintEst)+500,dim(m)[1]*(dim(m)[1]-1)/2))
  MOMrange <-
    estErrProbMethodOfMoments(nint=nint,nrec=totalObsR,nunr=totalObsU,ntot=dim(m)[1])
  closestpTP <-
    which.min(abs((1-pTP95U)-c(MOMrange[,3],MOMrange[,5])))
  pFP95U <- c(MOMrange[,2],MOMrange[,4])[closestpTP]
  
  
                                        #find lower bound for pFP
  
  pTP95L <- globalpTP-1.96*globalpTPse
  
  nintEst <- totalObsR/(pTP95L^2)
                                        #print(nintEst)
  nint <- (max(floor(nintEst)-500,100)):(min(floor(nintEst)+500,dim(m)[1]*(dim(m)[1]-1)/2))
  MOMrange <-
    estErrProbMethodOfMoments(nint=nint,nrec=totalObsR,nunr=totalObsU,ntot=dim(m)[1])
  closestpTP <-
    which.min(abs((1-pTP95L)-c(MOMrange[,3],MOMrange[,5])))
  pFP95L <- c(MOMrange[,2],MOMrange[,4])[closestpTP]
  
  
  list(globalpTP=globalpTP,
       globalpTPse=globalpTPse,
       globalpFP=globalpFP,
       pTP95L = pTP95L,
       pTP95U = pTP95U,
       pTPBSse = pTPBSse,
       pTPBS95L = pTPBS95L,
       pTPBS95U = pTPBS95U,
       pFP95L = pFP95L,
       pFP95U = pFP95U,
       nEligComplexes=nEligComplexes,
       nEligBaits=nEligBaits,
       nEligEdges=nEligEdges,
       nBaitsInComplexes=nBaitsInComplexes,
       complexSizes=complexSizes)
  
  
}


                                        #Wolfgang and Tony's functions

getDegrees = function(m){
  stopifnot(all(m %in% 0:1))
  m = m>0
  cbind(nr=rowSums(m & t(m)),no=rowSums(m & (!t(m))),
        ni=rowSums((!m) & t(m)))
}

estimatepTP <- function(n,r,u){
  
  stopifnot(n>0)
  pTPhat <- (2*r+u)/(2*n)
  list(pTPhat=pTPhat)

}

estimatepFP <- function(tE,r,u,pTP){
  
  Nnum <- (2*r+u)^2-4*r*tE
  Ndenom <- 4*pTP*(2*r+u)-4*r-4*tE*(pTP^2)
  N <- Nnum/Ndenom
  
  print(N)
  pFPhat <- ((2*r+u)-2*pTP*N)/(2*(tE-N))
  list(pFPhat=pFPhat)
}

obsProp <- function(x,y, obsPropThresh=.7){
  xn <- names(which(x==1))
  ans <-  length(intersect(xn,y))/sum(x) >= obsPropThresh
  ans
}
