map2Systematic2 <- function(allProt, tableList, sWAC){

    
    ##The baitSystematic vector will be a mapping of the Intact codes to Yeast Systematic Names
    #protSystematic <- vector(length = length(allProt))
    #yeast2Sys <- as.list(YEASTCOMMON2SYSTEMATIC)
    #yeast2Sys <- yeast2Sys[!is.na(yeast2Sys)]
    #yeastAlias <- names(unlist(as.list(YEASTALIAS)))
    #yro <- unlist(as.list(YEASTREJECTORF))
    #notfound <- vector()
    #notfoundSGD <- vector()
    #notfound2 <- vector()
    #notfound2SGD <- vector()
    
    #for(i in 1:length(allProt)){

        ProtN <- allProt
        whichProt <- tableList[["ac2xref"]][, "ac"] %in% ProtN
        subTable <- tableList[["ac2xref"]][whichProt,1:4]
        ###print(subTable)
        onlyWantAN = which(subTable[,"db"] == "uniprotkb")
        ac2ANCode = split(subTable[onlyWantAN,4], subTable[onlyWantAN,1])
        ac2ANCode <- lapply(ac2ANCode, unique)
        AN <- unlist(ac2ANCode)
        AN = unique(AN)
        if(class(AN) == "character"){
            AN <- strsplit(AN, "_yeast", fixed=TRUE)
            AN <- unlist(AN)
            AN <- toupper(AN)
        }

        
        if(!is.null(AN)){
          
          AN <- sWAC[AN]
          }
        
        if(length(AN)>0){
          return(AN)
        }
        else{
          return(allProt)
        }
    }

