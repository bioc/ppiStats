map2Systematic <- function(allProt, tableList, sWAC){
    #print(allProt)

    ##The baitSystematic vector will be a mapping of the Intact codes to Yeast Systematic Names
    protSystematic <- vector(length = length(allProt))
    yeast2Sys <- as.list(YEASTCOMMON2SYSTEMATIC)
    yeast2Sys <- yeast2Sys[!is.na(yeast2Sys)]
    yeastAlias <- names(unlist(as.list(YEASTALIAS)))
    yro <- unlist(as.list(YEASTREJECTORF))
    #notfound <- vector()
    #notfoundSGD <- vector()
    #notfound2 <- vector()
    #notfound2SGD <- vector()
    one2Many <- vector("list", length = length(allProt))
    
    for(i in 1:length(allProt)){

        ProtN <- allProt[i]
        
        whichProt <- tableList[["ac2xref"]][, "ac"] %in% ProtN
        subTable <- tableList[["ac2xref"]][whichProt,1:4]
        #print(subTable)
        onlyWantSgd = which(subTable[,"db"] == "sgd")
        ac2SgdCode = split(subTable[onlyWantSgd,4], subTable[onlyWantSgd,1])
        ac2SgdCode <- lapply(ac2SgdCode, unique)
        sgdC <- unlist(ac2SgdCode)
        sgdC <- unique(sgdC)
        #print(sgdC)
        #print(ProtN)
        #print(i)
        
        if (length(sgdC) == 1){

            if (!is.null(yeast2Sys[[sgdC]])){
                protSystematic[i] <- yeast2Sys[[sgdC]][1]
                
                one2Many[[i]] <- yeast2Sys[[sgdC]]
                
            }
            
            else {

                
                
                if(sgdC %in% yeastAlias || sgdC %in% yro){
                    protSystematic[i] <- sgdC
                    one2Many[[i]] <- sgdC
                }
                
                else{

                 

                  if(substr(sgdC,1,1)=="Y"){
                        protSystematic[i] <- sgdC
                        one2Many[[i]] <- sgdC
                      }
                  else{

                    aN <- map2Systematic2(ProtN, tableList, sWAC)
                    aN <- unlist(aN)
                      if(!is.na(aN)){
                          one2Many[[i]] <- aN
                      }
                      else{
                          one2Many[[i]] <- "ND"
                      }
                    if(!is.na(aN)){
                        if(substr(aN,1,1)=="Y"){
                            protSystematic[i] = aN
                        }
                        else{
                            protSystematic[i] <- ProtN
                        }
                    }
                    else{
                      protSystematic[i] <- ProtN
                  }
                }
              }
            }
        }
        
        else{
            #print("here")
            if(length(sgdC)==0){
                if(is.null(sgdC)){
                    aN <- map2Systematic2(ProtN, tableList, sWAC)
                    aN <- unlist(aN)
                    if(!is.na(aN)){
                        one2Many[[i]] <- aN
                    }
                    else{
                        one2Many[[i]] <- "ND"
                    }
                    #print(one2Many)
                    if(!is.na(aN)){
                        if(substr(aN,1,1)=="Y"){
                            protSystematic[i] = aN
                        }
                    }


                    
                }

                else{
                    #print("there")
                    protSystematic[i] <- ProtN
                    one2Many[[i]] <- sgdC
                }
            }
            else {
                #print("where")
                one2Many[[i]] <- list()
                
                for(k in 1:length(sgdC)){
                    if(!is.null(yeast2Sys[[sgdC[k]]])){
                        #print("i should get here once")
                        one2Many[[i]][[sgdC[k]]] <- yeast2Sys[[sgdC[k]]]
                
                 
                 
                    }
                    else{

                        if(substr(sgdC[k],1,1)=="Y"){
                            one2Many[[i]][[sgdC[k]]] <- sgdC[k]
                        }
                        else{
                            one2Many[[i]][[sgdC[k]]] <- "ND"
                        }
                    }
                
                }
            }
        }
        #print(one2Many[[ProtN]])
    }
    

    #protSystematic
    names(one2Many) = allProt
    one2Many

}
