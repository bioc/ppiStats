idHomodimers <- function(bpGraph){

    eL <- edges(bpGraph)

    nodes <- names(eL)

    homodimers <- mapply(function(x,y) x%in%y, nodes, eL)
    

    homodimers <- homodimers[homodimers]
    #print(length(homodimers))
    return(names(homodimers))
    
    
}
