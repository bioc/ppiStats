estErrProbMethodOfMoments = function(nint, nrec, nunr, ntot) {
  nEdges = ntot*(ntot+1)/2
  stopifnot(length(nrec)==1, length(nunr)==1, length(ntot)==1, 
    nrec+nunr<=nEdges, all(nint<=nEdges))

  nnon = nEdges-nrec-nunr
  n2 = nEdges-nint
  delta = (nnon-n2) - (nrec-nint)
  a = nEdges
  b = delta - 2*nint
  c = nint + delta^2/(4*n2) - nint/n2*nnon
  discr  = b*b-4*a*c

  ## just to double-check:
  discr2 = -delta^2*nint/n2 + 4*nint*(nrec-nint+nint*nnon/n2)
  ## plot(discr, discr2)
  stopifnot(all(abs(discr-discr2)/abs(discr+discr2+1e-10)<1e-10))

  p1 = (-b + sqrt(discr)) / (2*a)
  p2 = (-b - sqrt(discr)) / (2*a)

  q = function(p) ((nnon-n2)-(nrec-nint)+2*n2*p)/(2*nint)
  cbind(nint=nint, pfp1=p1, pfn1=q(p1), pfp2=p2, pfn2=q(p2))
}

