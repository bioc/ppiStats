assessSymmetry = function(bpMat, bpGraph = FALSE,
  pLevels = c(1e-06, 1e-04, 0.01)) {
    if (bpGraph){
        bpMat <- as(bpMat, "matrix")
    }
    
    stopifnot(all(bpMat %in% 0:1))
    bpMat = bpMat>0
    deg <- cbind(nr=rowSums(bpMat&t(bpMat)), no=rowSums(bpMat&(!t(bpMat))),
                 ni=rowSums((!bpMat)&t(bpMat)))

     
    nunrec = deg[, "no"]+deg[, "ni"]
    nmini  = pmin(deg[, "no"], deg[, "ni"])
    p = pbinom(nmini, size=nunrec, prob=0.5)
    p = pmin(2*p, 1) 
    nmax = 2*max(deg[, "no"], deg[, "ni"])
    contours = sapply(pLevels, qbinom, size=0:nmax, prob=0.5)
    list(deg=deg, p=p, contours=contours)
    
}

