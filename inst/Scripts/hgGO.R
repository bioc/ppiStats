#attached base packages:
#[1] "splines"   "tools"     "methods"   "stats"     "graphics"  "grDevices"
#[7] "utils"     "datasets"  "base"     

#other attached packages:
#     xtable     ppiData   Rgraphviz geneplotter    ppiStats       YEAST 
#    "1.4-1"     "0.0.8"    "1.13.0"    "1.13.0"     "1.1.1"   "1.11.14" 
#    GOstats    Category  genefilter    survival        KEGG        RBGL 
#    "2.1.4"     "2.1.6"    "1.13.1"      "2.29"    "1.12.0"    "1.11.0" 
#   annotate     Biobase       graph          GO 
#   "1.13.2"    "1.13.4"    "1.13.1"    "1.13.0" 


##Loading the libraries
library("GO")
library("ppiStats")
library("ppiData")
library("YEAST")
library("GOstats")
library("Category")
yeastGenome <- names(as.list(YEASTALIAS))
##

##Set up to conduct the HyperGeometric tests
bpGraphs <- lapply(bpExperimentNames, get) 
vbp <- mapply(function(x,y){intersect(x,y)},viableBaits,viablePrey)
vbpGraph <- mapply(function(x,y){subGraph(x, get(y))},vbp,bpExperimentNames)

##partition proteins in stochastic and systematic
vbpStochastic <- lapply(vbpGraph, function(x) {if(length(nodes(x))!=0) idStochastic(x, bpGraph=TRUE)})
vbpSystematic <- mapply(function(x,y){setdiff(nodes(x),y)},vbpGraph,vbpStochastic)

vbpSysNonTriv <- which(sapply(vbpSystematic,length) != 0)
vbpSys <- vbpSystematic[vbpSysNonTriv]
vbp1 <- vbp[vbpSysNonTriv]

viableBaits <- viableBaits[-5]
viablePrey <- viablePrey[-5]

GOontologies <- c("CC","BP","MF")
interestingProteins <- list()
interestingProteins[[1]] <- viableBaits
interestingProteins[[2]] <- viablePrey
names(interestingProteins) <- c("vBaits","vPrey")

##A function to build the parameter class:

buildParams4GO <- function(geneSet, universe, direction="over",
                           ontology = "CC", cond=TRUE, pThresh = 0.01){
  parameter <- new("GOHyperGParams", geneIds = geneSet,
                   universeGeneIds = universe, annotation="YEAST",
                   ontology = ontology, conditional=cond,
                   testDirection = direction, pvalueCutoff = pThresh)

  return(parameter)
}

hgTest4GO <- function(parameter, filename, append=TRUE,
                      label = "Experiment name here",
                      typeGeneSet = "Describe the gene set here",
                      cs=50){

  results <- hyperGTest(parameter)
  keep <- any(pvalues(results)<results@pvalueCutoff)
  htmlReport(results, file = paste(filename, ".html", sep=""),
             append = append, label = paste(label,": ", typeGeneSet, sep=""),
             categorySize = cs)
}

filenames <- dNames <- sub("-$", "", sub("BPGraph", "-", bpExperimentNames))
filenames <- filenames[-5]
direction <- c("over", "under")

for(dir in direction){

  for(ont in GOontologies){
    for (i in length(interestingProteins)){
      for(j in length(interestingProteins[[i]])){
        par <- buildParams4GO(interestingProteins[[i]][[j]], universe=yeastGenome,
                       direction=dir, ontology=ont)
        hgTest4GO(par, filename = filenames[j], label = filenames[j],
                  typeGeneSet = names(interestingProteins)[i], append=TRUE)
      }
    }
      
  }
  
}






for(direction in c("over", "under")) {
  ##------------------------------------------------------------------------------------------
  cat("\n\nThis code chunk will build the hypergeometric parameters for both the vbaits and vprey ",
      "wrt the genome over all ontologies searching for ", direction, "-represented categories.\n", sep="")
  ##------------------------------------------------------------------------------------------
  for(i in seq_along(GOontologies)){
    cat(GOontologies[i], "")
    for(j in seq_along(interestingProteins)){
      cat(names(interestingProteins)[j], "")
      parameter <- lapply(interestingProteins[[j]], function(x) {
        new("GOHyperGParams", geneIds=x, 
            universeGeneIds=yeastGenome,
            annotation="YEAST", 
            ontology=GOontologies[i],
            conditional=TRUE,
            testDirection=direction,
            pvalueCutoff=0.01)})
      hgTest <- lapply(parameter, hyperGTest)
      ind2Keep <- sapply(hgTest, function(x) any(pvalues(x)<x@pvalueCutoff))
      hgTest <- hgTest[ind2Keep]
      browser()
      ##print(hgTest)
      ##hgTest <- hgTest[-5]
      Label <- names(hgTest)
      Label <- sub("-$", "", sub("BPGraph", "-", Label))
      mapply(function(x,y){htmlReport(x,file=sprintf("ppi%sRep.html", direction),
                                      append=!(i==1&j==1),
                                      label=paste(y,names(interestingProteins)[j]),
                                      categorySize=50)},hgTest,Label)
    } ## for j
  } ## for i
} ## for direction

ppiHG <- function(proteins, universe, ontology=NULL,
                  conditional=TRUE,
                  annotation="YEAST",pval=0.01,
                  allAnnotated=TRUE){
  
  goOntology <- c("GO:0005575","GO:0008150","GO:0003674")
  names(goOntology) <- c("CC","BP","MF")
  options(error=recover)  
  if(!allAnnotated){
   goNode <- goOntology[ontology]
   
   descendents <- mget(goNode, get(paste("GO",ontology,"OFFSPRING",sep="")))
   #print(descendents)
   descendents <- descendents[[1]]
   index <- descendents%in%names(as.list(YEASTGO2ALLPROBES))
   descendents <- descendents[index]
   realUniverse <- unique(unlist(mget(descendents,
                                      YEASTGO2ALLPROBES)))
   universe <- intersect(universe, realUniverse)
  }
  
  parameter <- new("GOHyperGParams", geneIds=proteins, 
        universeGeneIds=universe,
        annotation=annotation, 
        ontology=ontology,
        conditional=conditional,
        testDirection="over",
        pvalueCutoff=pval)

  parameter2 <- parameter
  testDirection(parameter2) <- "under"
  hg <- hyperGTest(parameter)
  hg2 <- hyperGTest(parameter2)
  hgList <- list(over=hg, under=hg2)
  hgList

}

ppiHG(proteins=viableBaits[[1]], universe=yeastGenome,
      ontology="CC",
      conditional=TRUE,
      annotation="YEAST",pval=0.01,
      allAnnotated=FALSE)

stop("Hollarodullioeh")

##This code chunk will build hypergeometric parameters for the systematic baits
##wrt viable baits over all ontologies searching for over-represented categories
for(i in 1:length(GOontologies)){

   parameter1 <- mapply(function(x,y){new("GOHyperGParams",
                                         geneIds=x,
                                         universeGeneIds=y,
                                         annotation="YEAST",
                                         ontology=GOontologies[i],
                                         conditional=TRUE,
                                         testDirection="over",
                                         pvalueCutoff=0.01)},vbpSys,vbp1)


    
    hgTest <- lapply(parameter1, hyperGTest)
    ind2Keep <- sapply(hgTest, function(x) sum(pvalues(x)<x@pvalueCutoff))
    hgTest <- hgTest[which(ind2Keep != 0)]
    Label <- names(hgTest)
    Label <- sub("-$", "", sub("BPGraph", "-", Label))
    mapply(function(x,y){htmlReport(x,file="ppiSystematicOverRep.html",
                                    label=paste(y,"SystematicBaits"),append=TRUE,
                                    categorySize=50)},hgTest,Label)
}


##This code chunk will build hypergeometric parameters for the systematic baits
##wrt viable baits over all ontologies searching for under-represented categories
for(i in 1:length(GOontologies)){

  parameter1 <- mapply(function(x,y){new("GOHyperGParams",
                                         geneIds=x,
                                         universeGeneIds=y,
                                         annotation="YEAST",
                                         ontology=GOontologies[i],
                                         conditional=TRUE,
                                         testDirection="under",
                                         pvalueCutoff=0.01)},vbpSys,vbp1)


    hgTest <- lapply(parameter1, hyperGTest)
    ind2Keep <- sapply(hgTest, function(x) sum(pvalues(x)<x@pvalueCutoff))
    hgTest <- hgTest[which(ind2Keep != 0)]
    Label <- names(hgTest)
    Label <- sub("-$", "", sub("BPGraph", "-", Label))
    mapply(function(x,y){htmlReport(x,file="ppiSystematicUnderRep.html",
                                    label=paste(y,"SystematicBaits"),
                                    categorySize=50)},hgTest,Label)
}
