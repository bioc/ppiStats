\name{calcInOutDegStats}
\alias{calcInOutDegStats}

\title{This function calculates the various degree statistic for a
Protein-Protein Interaction (PPI) Graph}
\description{
	This function takes a graphNEL argument and calculates several degree
	statistics for which further analysis can be made.
}
\usage{
calcInOutDegStats(graphObj, homodimer=FALSE)
}

\arguments{
  \item{graphObj}{The argument \code{graphObj} is an instance of the 
    class graphNEL.}
  \item{homodimer}{A logical. If FALSE, the function removes all
    homodimer relationships. It retains them otherwise.}
}
\value{
The return value is a list of various degree statistics of the PPI directed
graph:
  \item{inDegree}{A named numeric vector. The names corresponds to the particular
  yeast Open Reading Frame (ORF), and each entry details how many edges flow towards
  that node, i.e. how many times each node was detected as a prey.}
  \item{outDegree}{A named numeric vector. The names corresponds to the particular
  yeast Open Reading Frame (ORF), and each entry details how many edges flow out from
  that node, i.e. how many prey was detected by each node as a bait.}
  \item{inDegreeMinusOutDegree}{A named numeric vector. This entry is the inDegree
  minus the outDegree entries, or the signed difference for each node as a prey to 
  each node as a bait.}
  \item{outDegreeMinusInDegree}{A named numeric vector. This entry is the outDegree
  minus the inDegree entries, or the signed difference for each node as a bait to 
  each node as a prey.}
  \item{recipInDegree}{A named numeric vector. Again the names represents the ORFs,
  while each entry corresponds to how many times each node was found as a prey in a 
  reciprocated interaction.}
  \item{recipOutDegree}{A named numeric vector. Again the names represents the ORFs,
  while each entry corresponds to how many times each node was found as a bait in a 
  reciprocated interaction.}
  \item{totalRecipDegree}{A numeric. For the PPI graph given, this argument gives the
  total number of reciprocated interactions found.}
  \item{unrecipInDegree}{A named numeric vector. The names corresponds ORF and 
  each entry details how many edges flow towards each node modulo the reciprocated 
  interactions, i.e. for each node $p$, how many baits $b$ detected $p$ as a prey with 
  the condition that $p$ never detects $b$ as a prey.}
  \item{unrecipOutDegree}{A named numeric vector. The names corresponds ORF and 
  each entry details how many edges flow out from each node modulo the reciprocated 
  interactions, i.e. for each node $b$, how many prey $p$ does $b$ detect with the 
  condition that $b$ is never detected by $p$.}
  \item{totalUnrecipDegree}{A numeric. The total number of unreciprocated interactions
  for the PPI graph given.}
  ...
}
\author{T Chiang}
\examples{
	library(ppiData)
	degStat <- calcInOutDegStats(Cagney2001BPGraph)
	degStat$recipOutDeg
}
\keyword{datagen}
