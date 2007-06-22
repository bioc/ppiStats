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
names(interestingProteins) <- c("Viable Baits","Viable Prey")

filenames <- dNames <- sub("-$", "", sub("BPGraph", "-", bpExperimentNames))
filenames <- filenames[-5]
direction <- c("over", "under")

for(dir in direction){

  for(ont in GOontologies){
    for (i in 1:length(interestingProteins)){
      for(j in 1:length(interestingProteins[[i]])){
        print(paste("Building the parameter instance for the ", names(interestingProteins)[i],
                    "of the dataset from", filenames[j], sep=" "))
        par <- ppiBuildParams4GO(interestingProteins[[i]][[j]], universe=yeastGenome,
                       direction=dir, ontology=ont)
        ###print(paste("Conducting the HG test for the", names(interestingProteins)[i],
        #            "of the dataset from", filenames[j], "on the", ont, "ontology",sep=" "))
        res <- ppiHGTest4GO(par, filename = filenames[j], label = filenames[j],
                  typeGeneSet = names(interestingProteins)[i], append=TRUE)
      }
    }
      
  }
  
}

