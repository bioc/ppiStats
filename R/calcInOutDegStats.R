calcInOutDegStats <- function(graphObj)
{
  deg <- degree(graphObj)
  inMinusOut <- deg$inDegree - deg$outDegree
  outMinusIn <- deg$outDegree - deg$inDegree

  gAM <- as(graphObj, "graphAM")

  adjMat <- gAM@adjMat

  recip <- adjMat*t(adjMat)
  recipIn <- colSums(recip)
  recipOut <- rowSums(recip)
  unrecip <- adjMat - recip
  unrecipIn <- colSums(unrecip)
  unrecipOut <- rowSums(unrecip)
  names(unrecipOut) <- names(unrecipIn)
  totalUnrecip <- sum(unrecip)
  
  degStat <- list()
  degStat$inDegree <- deg$inDegree
  degStat$outDegree <- deg$outDegree
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
