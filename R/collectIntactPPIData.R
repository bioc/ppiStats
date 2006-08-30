collectIntactPPIData <- function(intactID = c("EBI-375746", "EBI-531419", "EBI-295760", "EBI-698096",
                                      "EBI-592695", "EBI-476385", "EBI-476699", "EBI-493706", "EBI-620118",
                                      "EBI-619785", "EBI-492533", "EBI-492535", "EBI-697014", "EBI-538313",
                                      "EBI-538324", "EBI-531492", "EBI-526219", "EBI-491968", "EBI-707901",
                                      "EBI-74070",  "EBI-603898", "EBI-597864", "EBI-455820", "EBI-457380",
                                      "EBI-457455", "EBI-79957",  "EBI-603955", "EBI-600030", "EBI-600812",
                                      "EBI-601450", "EBI-75443",  "EBI-727915", "EBI-607722", "EBI-525791",
                                      "EBI-476132", "EBI-49802",  "EBI-491628", "EBI-491642",
                                      "EBI-783101", "EBI-762635",
                                      "EBI-389903", "EBI-392769")){


    options(error=recover)
    #fileToRead <- system.file("data", "tableList.rda", package = "y2hStat")
    #load(fileToRead)
    data(tableList)
    data(sWAC2Sys)
    
    sLabel <- which(tableList[["acInfo"]][, "ac"] %in% intactID)
    names(sLabel) <- tableList[["acInfo"]][sLabel, "ac"]
    shortL <- unique(tableList[["acInfo"]][sLabel[intactID], "shortLabel"])
    shortLabel <- shortL
        
    n <- length(intactID)


    ##We chose all the id's that were either labelled bait or prey
    isBait <- tableList[["interaction2interactor"]][, "role"] == "bait"
    isPrey <- tableList[["interaction2interactor"]][, "role"] == "prey"
    
    ##For each experiment, we wanted to collect every interaction involved;
    ##So for each experiment, there will be k interactions...we want them all
    ##Here we simple get the row indices...

    isExp <- lapply(intactID, function(x){tableList[["experiment2interaction"]][, "experiment"] == x})

    ##From the indices collected above, we subset the table so all the interactions of
    ##a particular experiment are collected:

    interaction <- lapply(isExp, function(x)
      {tableList[["experiment2interaction"]][x,"interaction"]})

    ##Now we go to another table. We want to find all the players in the interactions we
    ##have collected above. The following returns a logical...we will get the row index
    ##if the interaction is the one in which we are interested: 

    interactor <- lapply(interaction, function(x)
      {tableList[["interaction2interactor"]][,"interaction"]%in% x})
    
    ##Now we go into the same table as above; for each interaction, we want only those
    ##players corresponding to that interaction...since this is y2h, everyone element
    ##should only have two items...here we merely get the row indices

    b2p <- lapply(interaction,function(x){
        lapply(x, function(y){which(tableList[["interaction2interactor"]][, "interaction"] == y)})})
    
    ##In this list, we systematically build our interactions. We take the double row numbers
    ##collected above, and create a list of interactions.
    
    b2pList <- lapply(b2p, function(x){
        lapply(x, function(y){tableList[["interaction2interactor"]][y,]})})
                      
    ##This checks to see which element in b2pList is the bait

    check1 <- lapply(b2pList, function(x) {lapply(x, function(y){which(y[,"role"] == "bait")})})
                     

    ##This checks for the prey

    check2 <- lapply(b2pList, function(x){lapply(x, function(y){which(y[,"role"] == "prey")})})

    ##Later we will use the bait-prey pairs to build either a data.frame, an adjacency
    ##matrix, or a list of matrices. From the b2pList, we simple extract the bait element
    ##and then the prey element.
    indexSetAll <- vector("list", length = n)
    for (i in 1:n){

        indexSet <- vector("list", length = length(b2pList[[i]]))
        for(j in 1:length(indexSet)){
            indexSet[[j]] = c(b2pList[[i]][[j]][check1[[i]][[j]], "interactor"],
                              b2pList[[i]][[j]][check2[[i]][[j]], "interactor"])
            #print(indexSet[[j]])
        }

        indexSetAll[[i]] <- indexSet
    }

    names(indexSetAll) <- shortL

    ##We create a list of all the baits for each experiment.

    baits <- lapply(lapply(interactor, function(x)
                           {tableList[["interaction2interactor"]][(x & isBait),"interactor"]}), unique)

    ##We create a list of all the prey for each experiment

    preys <- lapply(lapply(interactor, function(x)
                           {tableList[["interaction2interactor"]][(x & isPrey),"interactor"]}), unique)


    ##We take the union of all the baits

    allBaits <- unique(unlist(baits))

    ##We take the union of all the preys

    allPreys <- unique(unlist(preys))


    numBaits = length(allBaits)
    numPreys = length(allPreys)

    ##This method will still leave the user with certain proteins that need
    ##to be mapped by hand!

    
    baitsSystematic <- map2Systematic(allProt = allBaits, tableList = tableList, sWAC = sWAC2Sys)
    preysSystematic <- map2Systematic(allProt = allPreys, tableList = tableList, sWAC = sWAC2Sys)


    dataList <- list()
    
    dataList$allBaits <- allBaits
    dataList$allPreys <- allPreys
    dataList$indexSetAll <- indexSetAll
    dataList$baitsSystematic <- baitsSystematic
    dataList$preysSystematic <- preysSystematic
    dataList$shortLabel <- shortLabel
        
    dataList

}

