separateExptBySize <- function(listOfGraphs, bound1, bound2){

  size <- sapply(listOfGraphs, function(x){
    y <- degree(x);
    if(isDirected(x))
      {z <- sum(y$inDegree)}
    else 
      {z <- sum(y)};
    return(z)})

  small <- listOfGraphs[size<bound1]
  medium <- listOfGraphs[size>=bound1 && size<bound2]
  large <- listOfGraphs[size>=bound2]

  partitionList <- list()
  partionList$smallScale <- small
  partitionList$mediumScale <- medium
  partitionList$largeScale <- large

  return(partitionList)

}
