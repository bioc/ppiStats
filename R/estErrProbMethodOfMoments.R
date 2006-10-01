estErrProbMethodOfMoments = function(n1, x1, x2, ntot) {
  nEdges = ntot*(ntot-1)
  stopifnot(length(x1)==1, length(x2)==1, length(ntot)==1, 
    x1+x2<=nEdges, all(n1<=nEdges))

  nnon = nEdges-x1-x2
  n2 = nEdges-n1
  delta = (nnon-n2) - (x1-n1)
  a = nEdges
  b = delta - 2*n1
  c = n1 + delta^2/(4*n2) - n1/n2*nnon
  discr  = b*b-4*a*c

  p1 = (-b + sqrt(discr)) / (2*a)
  p2 = (-b - sqrt(discr)) / (2*a)

  Pfn = function(p) ((nnon-n2)-(x1-n1)+2*n2*p)/(2*n1)
  cbind(n1=n1, pfp1=p1, pfn1=Pfn(p1), pfp2=p2, pfn2=Pfn(p2))
}

