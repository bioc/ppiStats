\name{assessSymmetry}
\alias{assessSymmetry}


\title{A function that calcualtes some directed degree statistics on
  graphs using a binomial error model}
\description{
  This function takes in a bait to prey protein-protein interaction
  graph (undirected) and calculates the reciprocated degree, the
  unreciprocated in and out degrees. Using this information and binomial
  error model, it asses the p-value for the in and out degree of each
  protein. Lastly, it plots the countour curves for these p-values
}
\usage{
assessSymmetry(bpMat, bpGraph = FALSE,
  prob=0.5, pLevels = 1e-4)
}

\arguments{
  \item{bpMat}{Either a bait to prey directed graphNEL or its
    corresponding adjacency matrix.}
  \item{bpGraph}{A logical. If TRUE, than bpMat is passed in by the user
    as a graphNEL.}
  \item{pLevels}{A numeric vector. It gives the levels to calculate the
    countours of the function in p in the (n-in, n-out)-plane}
  \item{prob}{A numeric. The bias of the coin used in the function
    pbinom call.}		
}

\value{
  A list:

  \item{deg}{A 3xn matrix. The rows are indexed by each protein. Column
    one gives the number of reciprocated edges; column two gives the
    number of unreciprocated out edges; colunm three gives the number of
    unreciprocated in-edges}
  \item{p}{The p-value for each protein with experimental in and out degrees}
  \item{countours}{The countours as a function of p}
}

\author{W Huber}

\examples{
library(ppiData)
assessSymmetry(Ito2001BPGraph, bpGraph=TRUE)
}
\keyword{datagen}
