calcInOutDegStats <- function(graphObj, homodimer=FALSE)
{
  #deg <- degree(graphObj)
  #inMinusOut <- deg$inDegree - deg$outDegree
  #outMinusIn <- deg$outDegree - deg$inDegree

  gAM <- as(graphObj, "graphAM")

  adjMat <- gAM@adjMat
  if(!homodimer) diag(adjMat) <- 0
  rownames(adjMat) <- colnames(adjMat)

  inDegree <- colSums(adjMat)
  outDegree <- rowSums(adjMat)
  inMinusOut <- colSums(adjMat) - rowSums(adjMat)
  outMinusIn <- rowSums(adjMat) - colSums(adjMat)
  
  recip <- adjMat*t(adjMat)
  recipIn <- colSums(recip)
  recipOut <- rowSums(recip)
  unrecip <- adjMat - recip
  unrecipIn <- colSums(unrecip)
  unrecipOut <- rowSums(unrecip)
  names(unrecipOut) <- names(unrecipIn)
  totalUnrecip <- sum(unrecip)
  
  degStat <- list()
  degStat$inDegree <- inDegree
  degStat$outDegree <- outDegree
  degStat$inDegreeMinusOutDegree <- inMinusOut
  degStat$outDegreeMinusInDegree <- outMinusIn
  #degStat$unrecip <- unrecip
  degStat$recipIn <- recipIn
  degStat$recipOut <- recipOut
  degStat$totalRecip <- sum(recip)/2
  degStat$unrecipInDegree <- unrecipIn
  degStat$unrecipOutDegree <- unrecipOut
  degStat$totalUnrecipDegree <- totalUnrecip
  degStat

}
