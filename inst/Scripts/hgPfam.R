library(GO)
library(Category)
library(ppiStats)
library(ppiData)
library("org.Sc.sgd.db")
library("GOstats")


yeastGenome <- names(as.list(org.Sc.sgdALIAS))

  ##------------------------------------------------------------------------------------------
  cat("\n\nThis first code chunk will calculate the hypergeometric tests on the PFam Domains
      to determine over-representation wrt viable baits/prey agains the genome. The resulting
      data is saved in the R objects overVBGen and overVPGen", sep="" )
  ##------------------------------------------------------------------------------------------

parameterVB <- lapply(viableBaits, function(x,y){new("PFAMHyperGParams",
                                      geneIds=x, 
                                      universeGeneIds=yeastGenome,
                                      annotation="org.Sc.sgd", 
                                      testDirection="over",
                                      pvalueCutoff=0.01)})


hgTestVB <- lapply(parameterVB, hyperGTest)
overVBGen <- lapply(hgTestVB, function(x) {names(which(pvalues(x)<0.01))})



parameterVP <- lapply(viablePrey, function(x,y){new("PFAMHyperGParams",
                                      geneIds=x, 
                                      universeGeneIds=yeastGenome,
                                      annotation="org.Sc.sgd", 
                                      testDirection="over",
                                      pvalueCutoff=0.01)})


hgTestVP <- lapply(parameterVP, hyperGTest)
overVPGen <- lapply(hgTestVP, function(x) {names(which(pvalues(x)<0.01))})

  ##------------------------------------------------------------------------------------------
  cat("\n\nNow we will set up the various parameters to conduct the hypergeometric tests
      for biased proteins against the viable bait/prey populations", sep="")
  ##------------------------------------------------------------------------------------------


bpGraphs <- lapply(bpExperimentNames, get) ## could use 'mget'
vbp <- mapply(function(x,y){intersect(x,y)},viableBaits,viablePrey)
vbpGraph <- mapply(function(x,y){subGraph(x, get(y))},vbp,bpExperimentNames)

vbpStochastic <- lapply(vbpGraph, function(x) {if(length(nodes(x))!=0) idStochastic(x, bpGraph=TRUE)})
vbpSystematic <- mapply(function(x,y){setdiff(nodes(x),y)},vbpGraph,vbpStochastic)

vbpSysNonTriv <- which(sapply(vbpSystematic,length) != 0)
vbpSys <- vbpSystematic[vbpSysNonTriv]
vbp1 <- vbp[vbpSysNonTriv]

viableBaits <- viableBaits[vbpSysNonTriv]
viablePrey <- viablePrey[vbpSysNonTriv]



  ##------------------------------------------------------------------------------------------
  cat("\n\nThis last code chunk will calculate the hypergeometric tests on the PFam Domains
      to determine over-representation wrt biased proteins agains the VB and VP. The resulting
      data is saved in the R objects overSysVB and overSysVP", sep="" )
  ##------------------------------------------------------------------------------------------



parameterSB <- mapply(function(x,y){new("PFAMHyperGParams",
                                      geneIds=x, 
                                      universeGeneIds=y,
                                      annotation="org.Sc.sgd", 
                                      testDirection="over",
                                      pvalueCutoff=0.01)},
              vbpSys,viableBaits)

hgTestSB <- lapply(parameterSB, hyperGTest)
overSysVB <- lapply(hgTestSB, function(x) {names(which(pvalues(x)<0.01))})


parameterSP <- mapply(function(x,y){new("PFAMHyperGParams",
                                      geneIds=x, 
                                      universeGeneIds=y,
                                      annotation="org.Sc.sgd", 
                                      testDirection="over",
                                      pvalueCutoff=0.01)},
              vbpSys,viablePrey)
hgTestSP <- lapply(parameterSP, hyperGTest)
overSysVP <- lapply(hgTestSP, function(x) {names(which(pvalues(x)<0.01))})

