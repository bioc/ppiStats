\name{idProteinErrorType}
\alias{idSystematic}
\alias{idStochastic}

\title{A function to identify those proteins affected by either
stochastic or systematic errors}
\description{
  This function takes in either a bait to prey Graph (matrix) and,
  based on a binomial error model, partitions proteins identified as
  either affected by systematic or stochastic error.
}
\usage{
idSystematic(bpMat, viable, bpGraph = FALSE, pThresh = 0.01, pLevels = 1e-4, prob=0.5)
idStochastic(bpMat, bpGraph = FALSE, pThresh = 0.01, pLevels = 1e-4, prob=0.5)
}

\arguments{
  \item{bpMat}{Either a bait to prey directed graphNEL or its
    corresponding adjacency matrix.}
  \item{viable}{This is a character vector of viable proteins. It is
    only used in the idSystematic function.}
  \item{bpGraph}{A logical. If TRUE, than bpMat is passed in by the user
  as a graphNEL.}
  \item{pThresh}{The p-value threshold for which to partition stochastic
  or systematic errors}
  \item{pLevels}{A numeric vector. It gives the levels to calculate the
    countours of the function in p in the (n-in, n-out)-plane}
  \item{prob}{A numeric. The bias of the coin in the call to pbinom.}
}

\value{
  A character vector of proteins either affected by systematic or
  stochastic errors.
}
\references{ ~put references to the literature/web site here ~ }
\author{T Chiang}

\examples{
library(ppiData)
idSystematic(Ito2001BPGraph, viableBaits[[1]], bpGraph=TRUE)
}
\keyword{datagen}
