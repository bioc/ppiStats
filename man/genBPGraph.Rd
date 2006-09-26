\name{genBPGraph}
\alias{genBPGraph}

\title{A function to generate the protein-protein interaction (ppi)
	  induced (un)directed graph}
\description{
  This function will take the ppi data and generate instance of the 
  class graph.
}
\usage{
genBPGraph(bpMat, directed=TRUE, bp=TRUE)
}

\arguments{
  \item{bpMat}{An adjacency matrix of PPI. If the matrix is obtained
  by empirical data, or bait to prey, then the  
  rows are indexed by the baits and the columns indexed by the 
  preys. If the rownames are not the same as the column names
  (i.e. a generic bait to prey matrix) the argument bp must be set to
  TRUE.}
  \item{directed}{A logical - if TRUE, the object will be a directed graph
  rather than an ordinary graph. For bait to prey interactions, this 
  parameter must always be set to TRUE.}
  \item{bp}{A logical - if TRUE, it signifies that the adjacency matrix 
  is a bait to prey empirically derived matrix so that the bait population
  (rownames) is usually different from the prey population (colnames).}
}

\value{

	An instance of the class graph.

}

\author{T Chiang}


\examples{
library(ppiData)
data(y2hSysGW)
eg=y2hSysGW[3]
genBPGraph(eg)
}
\keyword{datagen}
