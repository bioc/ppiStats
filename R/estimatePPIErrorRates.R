estimatePPIErrorRates <- function(intactPPIfiles=NULL,
                                  dataSets2BRemoved=c("EBI-1005374", "EBI-895356"),
                                  GS="default"){

  ##A function to estimate local error rates based on the
  ##methods of moments (cf Chiang et al, Genome Biology 2007)

  ##If we are going to use PSI 2.5 XML files (presently only IntAct),
  ##then we can parse these files with the psi25interaction function
  ##from Rintact
  if(!is.null(intactPPIfiles)){
    interactionEntries<-lapply(intactPPIfiles,
                              function(x) psi25interaction(x))
  }

  ##Next we filter out interactions that are not obtained via Y2H
  intactY2HList <- lapply(interactionEntries, getbpPairsFromInteractionEntry2)
  intactY2H <- do.call(rbind, intactY2HList)
  colnames(intactY2H) <- c("bait","prey","PMedID","IntActCode")

  ##Now we map the IntAct accession codes to the gene locus name if
  ##available
  map = sapply(interactionEntrys, function(x) x@interactors[, "locusName"])
  map2 = unlist(map)
  map2 = map2[!duplicated(names(map2))]
  
  baits = map2[intacty2h[,1]]
  prey = map2[intacty2h[,2]]
  intacty2h[,1] = baits
  intacty2h[,2] = prey

  ## remove those interactions without locusID
  ## remove interactions between the same protein

  intacty2h = intacty2h[!is.na(intacty2h[,1])&!is.na(intacty2h[,2]),]
  intacty2h = intacty2h[intacty2h[,1]!=intacty2h[,2],]
  
  ##print(dim(intacty2h))  ## 8289 by 4, Nov 25
  ## split protein pairs according to their expIntact ID
  
  intacty2hbpPairs = split(data.frame(intacty2h[,c("bait","prey")],
    stringsAsFactors=FALSE),intacty2h[,"IntactCode"])
  
  ## bpPairs2bpMatrix() is a function that translates protein pairs into
  ## two bpmatrices.  A bpmatrix is a matrix with rows representing
  ## baits and columns representing preys.  Each experiment corresponds
  ## to two bpmatrices, namely, test bpmatrix and pos bpmatrix. The former
  ## documents the tested bait prey pairs, and the latter documents the
  ## observed bait prey pairs.  Currently, we assume all the edges
  ## between baits and preys have been tested, we will modify this
  ## later.
  
  ##   it seems the output is a list of length = # of experiments, each
  ##   component is a list of length two, containing the two matrices
  ##   described above
  
  intacty2hbpMatrix = lapply(intacty2hbpPairs, bpPairs2bpMatrix)
  
  ## bpMatrix2ppmatrix() is a function that translates a bpmatrix to
  ## ppmatrix, There are also two ppmatrices for each experiment,
  ## namely, test ppmatrix and pos ppmatrix. The difference between
  ## bpmatrix and ppmatrix is that the latter doesn't care about the
  ## bait prey orders.  thus it is not directional, and only upper
  ## triangle of ppmatrix is filled with numbers.
  
  intacty2hppMatrix = lapply(intacty2hbpMatrix, bpMatrix2ppMatrix)

  ## intacty2hhtp documents all the experiments that more than 10 interactions 
  ## are observed
  
  expNum = table(intacty2h[,"IntactCode"])
  intacty2hhtp = names(expNum)[expNum>=10]
  intacty2hhtp2 <- setdiff(intacty2hhtp, data2BRemoved)

  ## expIntact2expPubMed documents the mapping between expIntact ID and PubMed ID
  expIntact2expPubMed = unique(intacty2h[,c("PMedID","IntactCode")])
  rownames(expIntact2expPubMed) = expIntact2expPubMed[,"IntactCode"]
  expIntact2expPubMed = expIntact2expPubMed[,"PMedID"]

  ## for all the small experiments and plus the above two datasets,
  ## we simply set tested edges equal observed edges for ppmatrix

  for(i in 1:length(intacty2hppMatrix)) {
    if(!names(intacty2hppMatrix)[i]%in% b2pall)
      {
        intacty2hppMatrix[[i]][["test"]] <- intacty2hppMatrix[[i]][["pos"]]
      }
  }
  
  ## we do the same thing for bpmatrix
  for(i in 1:length(intacty2hbpMatrix))
    {
      if(!names(intacty2hbpMatrix)[i]%in% b2pall )
        {
          intacty2hbpMatrix[[i]][["test"]]<-intacty2hbpMatrix[[i]][["pos"]]
        }
    }

  if(GS == "default"){
    mpactTable <- system.file("extdata","PPI_180052006", package="ppiStats")
    mpact = as.matrix(read.table(file=mpactTable,sep="|",
      stringsAsFactors=FALSE))
    colnames(mpact) = c("orf1","gene1","orf2","gene2","descr","ref","evi")
    evidence = strsplit(mpact[,"evi"],split=",")
    htp = sapply(evidence,function(x){any(x %in% "902.01.09.02")})
    physical = sapply(evidence,function(x){length(grep("902.01.01.02.01",x))>0;})
    mpactphysical = mpact[physical,]
    mpactGS = mpact[physical&!htp,]
    refIntacty2h = unique(intacty2h[,"expPubMed"])
    toRemain = sapply(strsplit(mpactGS[,"ref"],split=","),
      function(x)any(!x%in%refIntacty2h))
    mpactGSrl = mpactGS[toRemain,]
    refmpactsize = table(unlist(strsplit(mpactGSrl[,"ref"],split=",")))
    refmpactlargesize = names(refmpactsize)[refmpactsize>10]
    toRemain = sapply(strsplit(mpactGSrl[,"ref"],split=","),
      function(x)any(!x%in%refmpactlargesize))
    mpactGSrll = mpactGSrl[toRemain,]
    GS = toupper(mpactGSrll[,c("orf1","orf2")])
    GS = unique(t(apply(GS,1,sort)))
  }

  ## the function getFNrateFromGS tries to get the FNrates for each
  ## experiment without clustering neighboring experiments.
  
  ppMatrix = intacty2hppMatrix
  result = lapply(ppMatrix, function(x) getFNrateFromGS(x, GS))
  GStest = apply(do.call(cbind,lapply(result,function(x)x$GStest)),1,sum)
  GSpos = apply(do.call(cbind,lapply(result,function(x)x$GSpos)),1,sum)
  stats = do.call(rbind,lapply(result,function(x)x$stats))
  
  ##this seems to compute, for each experiment the number of
  ##interactions tested, and the number that were reported as positive
  
  exptestAndpos = do.call(rbind, lapply(ppMatrix, function(x) c(sum(x[["test"]]),
    sum(x[["pos"]]))))
  colnames(exptestAndpos) = c("testNum","posNum")
  
  ## the function getFNrateFromNeighbor() makes use of the above estimates and 
  ## also combining neighboring experiments when the datasets are too small.
  
  FNfromN = getFNrateFromNeighbor(stats[,"testNum"],stats[,"posNum"],
    exptestAndpos[,1], 50, exptestAndpos[,1])
  
  FN = FNfromN[,"FNrate"]
  
  ## get FP rate, given FN rate
  ##getx123() calculates the x1,x2,x3 values for each experiment.
  
  x123 = do.call(rbind, lapply(ppMatrix, getx123))
  
  ## getFPrateFromNeighbor() make use of x123 and the estimated FN rates to 
  ## estimate FP rate,
  ## small experiments are clustered together before estimation
  
  FPfromN = getFPrateFromNeighbor(x123[,1], x123[,2], x123[,3],
    exptestAndpos[,1], 200, exptestAndpos[,1], FN)
  FP = FPfromN[,"FPrate"]
  
  result = cbind(exptestAndpos, FNfromN, FPfromN)
  FPFN = result[order(result[,"testNum"],decreasing=TRUE),]

  return(FPFN)
}

## the function aims to extract edges that are detected using Y2H.
## input: an interactionEntry from the output of psi25interaction()
## output: a baitprey table with four columns, 
## namely,(x@bait,x@prey,x@expPubMed,x@expIntAct).
## Each row represents an interaction.
## FIXME: RG says the any that is used below looks suspicious -
##   I would be surprised if interactionType was of length more than 1?

getbpPairsFromInteractionEntry2 <- function(interactionEntry)
{
    interactionTypeWanted = c("two hybrid","two hybrid array", 
                               "2h fragment pooling", "2 hybrid")
    interactions = interactionEntry@interactions
    bpPairs = lapply(interactions,function(x){
                if(any(x@interactionType %in% interactionTypeWanted)){
                      cbind(x@bait,x@prey,x@expPubMed,x@expIntAct);}
              })
    do.call(rbind,bpPairs)
}

## The function aims to change baitprey table into a list of two matrix in 
## the bpmatrix format,where rows are baits and columns are preys.
## input: a baitprey table for a single experiment, from the output of 
##        getbpPairsFrominteractionEntry2()
## output: a list of two matrix,test and pos, for documenting how many times 
##         a directed edge is tested and observed respectivly.
## Both matrices are in the format of bpmatrix, where rows are baits and 
## columns are prey.

bpPairs2bpMatrix <- function(bpPairs)
  {
      if(any(colnames(bpPairs)!=c("bait","prey")))stop("bpPairs format error")
      vb<-sort(unique(bpPairs[,1]))
      vp<-sort(unique(bpPairs[,2]))
      pos<-matrix(0,length(vb),length(vp),dimnames=list(vb,vp))
      for(i in 1:nrow(bpPairs))
        {
            pos[bpPairs[i,1],bpPairs[i,2]]<-pos[bpPairs[i,1],bpPairs[i,2]]+1
        }
      test<-pos
      test[test==0]<-1
      ##remove the interactions between the same proteins from test matrix
      vbp<-intersect(vb,vp)
      if (length(vbp)>1) diag(test[vbp,vbp])<-0
      if (length(vbp)==1) test[vbp,vbp]<-0

      list(pos=pos,test=test)
  }


## the function aims transform the bpmatrix format  into ppmatrix format.
##  input: a list of two matrix,test and pos,which is the output of 
##         bppairs2bpMatrix().
##  output: a list of two matrix,test and pos, for documenting how many times 
##          an undirected edge is tested and observed respectivly.
##  Both matrix are in the format of ppmatrix,where only upper triangle 
##  has numbers.

 bpMatrix2ppMatrix <- function(bpMatrix)
  {
      if(any(names(bpMatrix)!=c("pos","test"))) stop("bpMatrix has wrong format");
      vb<-rownames(bpMatrix[[1]]);
      vp<-colnames(bpMatrix[[1]]);
      vbp<-sort(union(vb,vp));
      pos<-matrix(0,length(vbp),length(vbp),dimnames=list(vbp,vbp));
      test<-pos;
      pos[vb,vp]<-bpMatrix[["pos"]];
      pos[vp,vb]<-pos[vp,vb]+t(bpMatrix[["pos"]]);
      test[vb,vp]<-bpMatrix[["test"]];
      test[vp,vb]<-test[vp,vb]+t(bpMatrix[["test"]]);
      ## homoprotein interactions are not considered here
      diag(pos)<-0;
      diag(test)<-0;
      ## lower triangle doesn't contain information
      pos[lower.tri(pos)]<-0;
      test[lower.tri(test)]<-0;
      list(pos=pos,test=test);

  }

## this function estimates the FNrate for an experiment from 
## Goldstandard set.
##   input:
##     ppMatrix: a list of two matrix,pos and test, in the ppmatrix format.
##     GS: a table with two columns, each row represent an edge, 
##   output:
##      GStest: a vector documenting if an edge in GS has been tested.
##      GSpos: a vector documenting if an edge in GS has been observed.
##      Stats: a vector documenting the number of tested edges within 
##             GS(testNum),the number of observed edges within GS(posNum), 
##             and the FN rate(FN).

 getFNrateFromGS <- function(ppMatrix,GS) {
      vbp<-rownames(ppMatrix[[1]])
      GSpos<-apply(GS,1,function(x){if(all(x%in%vbp)){ 
               pos<-ppMatrix[["pos"]][x[1],x[2]]} else {pos<-0};pos;})
      GStest<-apply(GS,1,function(x){if(all(x%in%vbp)){ 
               test<-ppMatrix[["test"]][x[1],x[2]]} else {test<-0};test;})
      testNum=sum(GStest)
      posNum=sum(GSpos)
      stats<-c(testNum=testNum,posNum=posNum,FN=1-posNum/testNum)
      list(GStest=GStest,GSpos=GSpos,stats=stats)
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

 getFPrateFromFNrate <- function(x1,x2,x3,FNrate) {
      tot<-x1+x2+x3
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

## this function aims to estimate FNrate for each experiment by combining 
## neighboring expreriments when needed.
## input:
##   GStestNum: a vector documenting the number of tested edges within GS 
##              for each experiment.
##   GSposNum: a vector documenting the number of observed edges within GS 
##             for each experiment.
##   ExptestNum: a vector documenting the number of tested edges for each 
##               experiment.
##   cutoff: the cutoff of the number of tested edges in an experiment, 
##            below which neighboring experiments are combined.
##   newExptestNum: here it is the same as ExptestNum.
## output:
##    a table with two columns. Each row represents an experiment. The first 
##    column documents  FNrate and the second column documents how many edges 
##    are tested within GS in the combined datasets for that experiment.

 getFNrateFromNeighbor <- function(GStestNum, GSposNum, ExptestNum, cutoff, 
           newExptestNum) {
      FNrate<-rep(NA,length(newExptestNum))
      names(FNrate)<-names(newExptestNum)
      n<-FNrate
      for(j in 1:length(newExptestNum)) {
            distance<-abs(ExptestNum-newExptestNum[j])
            distanceOrder<-order(distance,decreasing=FALSE)
            totalGStestNum<-0
            totalGSposNum<-0
            for(i in 1:length(distanceOrder)) {
                  totalGStestNum<-totalGStestNum+GStestNum[distanceOrder[i]]
                  totalGSposNum<-totalGSposNum+GSposNum[distanceOrder[i]]
                  if(totalGStestNum>=cutoff) {
                       FNrate[j]<-1-totalGSposNum/totalGStestNum
                       if(FNrate[j]>0)break
                  }
            }
            n[j]<-totalGStestNum
      }
      cbind(FNrate,n);
  }


## this function aims to estimate FPrate for each experiment by combining 
## neighboring expreriments when needed
## input:
##   x1,x2,x3 are three vectors of values from getx123() for each experiment.
##   ExptestNum: a vector documenting the number of tested edges for each 
##               experiment.
##   cutoff: the cutoff of the sum of x1+x2+x3 for an experiment, below which 
##           neighboring experiments are combined.
##   newExptestNum: here it is the same as ExptestNum.
##   FNrate: a vector of FNrates for each experiment
## output:
##     a table with two columns. Each row represents an experiment. The first 
##     column documents  FPrate and the second column documents the sum of 
##     x1+x2+x3 in the combined datasets for that experiment.

 getFPrateFromNeighbor <- function(x1, x2, x3, ExptestNum, cutoff,
        newExptestNum, FNrate)
  {
      FPrate<-rep(NA,length(newExptestNum))
      names(FPrate)<-names(newExptestNum)
      n<-FPrate
      for(j in 1:length(newExptestNum)) {
            distance<-abs(ExptestNum-newExptestNum[j])
            distanceOrder<-order(distance,decreasing=FALSE)
            totalx1<-0
            totalx2<-0
            totalx3<-0
            for(i in 1:length(distanceOrder)) {
                  totalx1<-totalx1+x1[distanceOrder[i]]
                  totalx2<-totalx2+x2[distanceOrder[i]]
                  totalx3<-totalx3+x3[distanceOrder[i]]
                  if(totalx1+totalx2+totalx3>=cutoff) {
                        FPrate[j]<-getFPrateFromFNrate(totalx1, totalx2, 
                                         totalx3,FNrate[j])
                        if( FPrate[j]>0 ) break
                    }
              }
            n[j]<-totalx1+totalx2+totalx3
        }
      cbind(FPrate,n)
  }

## this function computes the Bayes factor for each edge in an experiment
## input:
##    testMatrix: a matrix documenting how many times an edge is tested for 
##         an experiment. Both bpmatrix format and ppmatrix format are OK here.
##    posMatrix: a matrix documenting how many times an edge is observed for 
##         an experiment.
##    FNrate: FN rate for that experiment.
##    FPrate: FP rate for that expeirment.

 getBayesFactor <- function(testMatrix, posMatrix, FNrate, FPrate) {
      if( !all(dim(testMatrix) == dim(posMatrix)) )
        stop("non-conformable matrices")
      BFpos = ((1-FNrate)/FPrate)
      BFneg = (FNrate/(1-FPrate))
      BFpos^posMatrix*BFneg^(testMatrix-posMatrix)
  }


## FIXME: this seems odd - need to check its use, as there will be 1s
##   where nothing was tested
## this function multiplies the entries in matrices, it expands them to
## so that the final result is of dimentions equal to the union of the row
## and column names of all inputs.
##
## input: any number of matrices, the matrices can be of various sizes, and 
##    with different row and column names.
## output: a matrix that is the product of input.

 combineMatrixByMultiply <- function(...) {
      matrixList = list(...)
      rows = sort(unique(unlist(lapply(matrixList,rownames))))
      columns = sort(unique(unlist(lapply(matrixList,colnames))))
      cMatrix = matrix(1, nrow = length(rows), ncol = length(columns),
                  dimnames=list(rows,columns))
      for(i in 1:length(matrixList)) {
            rowsubset = rownames(matrixList[[i]])
            columnsubset = colnames(matrixList[[i]])
            cMatrix[rowsubset,columnsubset] = cMatrix[rowsubset,columnsubset]*matrixList[[i]]
        }
      cMatrix
  }

## this function aims to add matrices
##  input: any number of matrix, the matrix can  be of various size, and with 
##         different row and column names.
##  output: a matrix that is the sum of input.

 combineMatrixBySum <- function(...)
  {
      matrixList<-list(...)
      rows<-sort(unique(unlist(lapply(matrixList,rownames))))
      columns<-sort(unique(unlist(lapply(matrixList,colnames))))
      cMatrix<-matrix(0,length(rows),length(columns),dimnames=list(rows,columns))
      for(i in 1:length(matrixList)) {
            rowsubset<-rownames(matrixList[[i]])
            columnsubset<-colnames(matrixList[[i]])
            cMatrix[rowsubset,columnsubset] = cMatrix[rowsubset,columnsubset]+
                                                matrixList[[i]]
        }
      cMatrix
  }

## this fuction aims to output if an undirected edge has been reciprocally 
## tested in an experiment, and how the observations.
## input: a list of two matrix,pos and test, in the bpmatrix format.
## output: a list of three matrix in ppmatrix format.
##    reciptest:a matrix documenting if an undircted edge has been 
##         reciprocally tested.
##    recippos: a matrix documenting if a reciprocally tested edge has been 
##              reciprocally observed.
##    unrecippos: a matrix documenting if a reciprocally tested edges has 
##                been unreciprocally observed.

 getRecipResults <- function(bpMatrix) {
      if( any(names(bpMatrix)!=c("pos","test")) ) 
                stop("bpMatrix has wrong format")
      ##whether bp pairs have been tested or positive.
      bpMatrix[["test"]]<-bpMatrix[["test"]]>0
      bpMatrix[["pos"]]<-bpMatrix[["pos"]]>0
      ##
      vb<-rownames(bpMatrix[[1]])
      vp<-colnames(bpMatrix[[1]])
      vbp<-sort(union(vb,vp))
      pos<-matrix(0,length(vbp),length(vbp),dimnames=list(vbp,vbp))
      test<-pos
      pos[vb,vp]<-bpMatrix[["pos"]]
      pos[vp,vb]<-pos[vp,vb]+t(bpMatrix[["pos"]])
      test[vb,vp]<-bpMatrix[["test"]]
      test[vp,vb]<-test[vp,vb]+t(bpMatrix[["test"]])
      ##homoprotein interactions are not considered here
      diag(pos)<-0
      diag(test)<-0
      ##lower triangle doesn't contain information
      pos[lower.tri(pos)]<-0
      test[lower.tri(test)]<-0
      #####
      reciptest<-test
      reciptest[,]<-0
      reciptest[test==2]<-1
      recippos<-pos
      recippos[,]<-0
      recippos[reciptest & pos==2]<-1
      unrecippos<-pos
      unrecippos[,]<-0
      unrecippos[reciptest & pos==1]<-1

      list(reciptest=reciptest,recippos=recippos,unrecippos=unrecippos)
  }


## this function makes a upper triangle matrix into a symmetrix matrix
 getMatrixsym <- function(x){
         x[lower.tri(x)]<-t(x)[lower.tri(x)]
         x
 }

