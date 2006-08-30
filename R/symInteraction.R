symInteraction <- function(y2h){
    #options(error=recover)

    sym <- vector("list", length=length(y2h))
    
    for(i in 1:length(sym)){
        baits <- names(y2h[[i]])
        sym[[i]] <- list()
        
        for(j in 1:length(baits)){
            
            for(k in 1:length(y2h[[i]][[baits[j]]])){
                
                preyOfInt <- y2h[[i]][[baits[j]]][k]
                
                if(length(preyOfInt) > 0){
                
                    if(!is.null(unlist(y2h[[i]][[preyOfInt]]))){
                        if (baits[j] %in% unlist(y2h[[i]][[preyOfInt]])){
                            
             
                            sym[[i]][[preyOfInt]] <- c(sym[[i]][[preyOfInt]],baits[j])
                        }
                        
                    }
                }
                
            }
        }
        
    }

    names(sym) <- names(y2h)
    sym
}
  
