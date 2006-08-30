intAct2Sys <- function(prot2Sys, bpSysL){
    options(error=recover)
    if(!is.null(names(bpSysL[[prot2Sys]])) && !is.na(names(bpSysL[[prot2Sys]]))){

        mapping <- intAct2Sys(names(bpSysL[[prot2Sys]])[1], bpSysL[[prot2Sys]])
   }

    else{

        if(!is.null(bpSysL[[prot2Sys]]) && bpSysL[[prot2Sys]] != "ND"){

            return(bpSysL[[prot2Sys]][1])

        }

        else{
            return(prot2Sys)
        }
          
        
    }

    mapping
    
}
