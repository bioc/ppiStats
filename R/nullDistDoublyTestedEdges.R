nullDistDoublyTestedEdges = function(deltaMax, n, pFP, pFN) {
  
  deltaMax = as.integer(deltaMax)
  n = as.integer(n)

  if( any(is.na(n)) || (length(n)!=1 ) || (n<2) )
    stop("'n' must be positive integer scalar >= 2.")
  if( any(is.na(deltaMax)) || (length(deltaMax)!=1 ) || (deltaMax>=n) )
    stop("'deltaMax' must be a positive integer scalar less than 'n'.")
  if( (!is.numeric(pFP)) || (length(pFP)!=1 ) || (pFP<0) || (pFP>1))
    stop("'pFP' must be numeric scalar between 0 and 1.")
  if( (!is.numeric(pFN)) || (length(pFN)!=1 ) || (pFN<0) || (pFN>1))
    stop("'pFN' must be numeric scalar between 0 and 1.")
  
  nMax = deltaMax + as.integer(qbinom(1-(1e-10), size=n, prob=max(pFP*pFP, 2*pFP*(1-pFP))))
  
  probmultinom = c((1-pFN)^2, 2*(1-pFN)*pFN, (pFN^2))
  p = array(as.numeric(0), c(nMax, nMax, deltaMax)+1)

  for(delta in 0:deltaMax) {    
    z = n-delta-1
    du  = dbinom(0:nMax, size=z, prob=2*pFP*(1-pFP))
    dr  = dbinom(0:nMax, size=z, prob=pFP*pFP)
    dfp = outer(dr, du)
    for(nu in 0:delta) {
      idu = 1:(nMax-nu+1)   ## index into du
      ipu = nu+idu          ## index into p (1st dimension)
      for(nr in 0:(delta-nu)) {
        ## dm is the multinomial probability that
        ## nr reciprocated edges and nu unreciprocated edges arise from the true edges
        dm = dmultinom(c(nr, nu, delta-nu-nr), size=delta, prob = probmultinom)

        if(dm > 1e-10) {
          idr = 1:(nMax-nr+1) ## index into dr
          ipr = nr+idr        ## index into p (2nd dimension)
          
          sdfp = dfp[idr, idu]
          ## we can use this hard-coded threshold because we know that probabilities live between 0 and 1
          if(abs(sum(sdfp)-1) > 1e-8) 
            stop(paste("Truncation error: the probabilities in the FP distributions do not sum up to 1.",
                       "Please try increasing 'nMax'.", sep="\n"))
          
          p[ipr, ipu, delta+1] = p[ipr, ipu, delta+1]  + dm*sdfp
        }## if
      } ## for nr
    } ## for nu
  } ## for delta
  return(p)
}
