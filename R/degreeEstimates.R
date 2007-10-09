##functions to maximize truncated multinomial pmf for estimating node degree

#degreePMF calculates the joint pmf for 
#reciprocated and unreciprocated degree 
#using multinomial models for TP and FP observations


degreePMF <- function(deltahat,rd,ud,pTP,pFP,N){

       tprobs <- c(pTP^2,2*pTP*(1-pTP),(1-pTP)^2)
       fprobs <- c(pFP^2,2*pFP*(1-pFP),(1-pFP)^2)

       numNonEdges <- N-deltahat
       

       lpmf <- matrix(NA,nrow=(min(rd,deltahat)+1),ncol=(min(ud,deltahat)+1))

       for (i in 0:min(rd,deltahat)){
       for (j in 0:min(ud,deltahat-i)){


           lpmf[i+1,j+1] <- 
		dmultinom(c(i,j,deltahat-i-j),prob=tprobs,log=TRUE) + 
		dmultinom(c(rd-i,ud-j,numNonEdges-(rd-i)-(ud-j)),
						prob=fprobs,log=TRUE) 
    
       }
       }

       lpmf <- lpmf - log(1-(1-pFP)^(2*numNonEdges)*(1-pTP)^(2*deltahat))

       pmf <- sum(exp(lpmf),na.rm=TRUE)
       pmf
}


#rd = observed reciprocated degree
#ud = observed unreciprocated degree
#pTP = true positive probability
#pFP = false positive probability
#N = total number of nodes in the graph

findDegree <- function(rd,ud,pTP,pFP,N){

	   stopifnot(rd+ud>0)

           deltahat <- rd
           diffLPMF <- 1

           while(diffLPMF>0){
           
           deltahat1 <- deltahat+1

           lpmf <- degreePMF(deltahat,rd,ud,pTP,pFP,N) 
           lpmf1 <- degreePMF(deltahat1,rd,ud,pTP,pFP,N) 
           diffLPMF <- lpmf1-lpmf
           if(diffLPMF>0) deltahat <- deltahat1    
           }

           deltahat
}


#degreeEstimates
#finds estimated degree for all nodes in bait-prey matrix


#m is a bait-prey matrix for the VBP graph (a square matrix)
#pTP and pFP are true and false positive probabilities
		
degreeEstimates <- function(m,pTP,pFP){

	nodeNames <- rownames(m)
	N <- length(nodeNames)
	degEst <- rep(0,length(nodeNames))
	names(degEst) <- nodeNames

	stopifnot(all(m %in% 0:1))
	m = m>0
	mRio <- cbind(nr=rowSums(m & t(m)),no=rowSums(m & (!t(m))),
			      ni=rowSums((!m) & t(m)))

 
	mR <- mRio[,"nr"]
	mU <- mRio[,"ni"] + mRio[,"no"]

	for (i in nodeNames){
	    if (mR[i]+mU[i]==0){ 
	       degEst[i] <- NA
	       } else
	    degEst[i] <- findDegree(rd=mR[i],ud=mU[i],pTP=pTP,pFP=pFP,N=N)	
	    }
	degEst
}
