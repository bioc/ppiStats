idViableProteins <- function(bpGraph){

    deg <- degree(bpGraph)
    vb <- names(which(deg$outDegree > 0))
    vp <- names(which(deg$inDegree > 0))
    vProt <- list()
    vProt$viableBaits <- vb
    vProt$viablePrey <- vp
    vProt
}
