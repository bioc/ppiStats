\name{estimatePPIErrorRates}
\alias{estimatePPIErrorRates}
\title{Estimate false positive and false negative error probabilities
  for direct protein interactions}
\description{Estimate false positive and false negative error probabilities
  for direct protein interactions using a protein interaction graph
  gold standard}
\usage{
estimatePPIErrorRates(matList,GSPos=NULL,GSNeg=NULL)
}
\arguments{
  \item{matList}{A named list of list. The names corresponds to a
    particular publication. Each element of the top list is a list of
    two matrices:
           1. The adjacency matrix of the positive interaction
           datapoints.
           2. The adjacency matrix of the tested interaction datapoints.
    The dimension names of these two matrices should be identical (an
    error will be thrown if this is not the case) where the rows are
    indexed by the bait proteins and the columns are indexed by the
    prey. The adjacency matrices are 0,1 matrices where an 1 entry
    signifies either a positive interaction relationship or a positive
    tested relationship. A 0 entry is then well defined.}	   
  \item{GSPos}{A positive gold standard protein interaction graph given by its
  adjacency matrix. The dimension names of this matrix must be the same
  identifying system used for the dimension names for the matrices found
  within the matList parameter. An entry of one in this matrix indicates
  a positive interaction while an entry of zero is inconclusive.}
  \item{GSNeg}{A negative gold standard protein interaction graph given by its
  adjacency matrix. The dimension names of this matrix must be the same
  identifying system used for the dimension names for the matrices found
  within the matList parameter. An entry of one in this matrix indicates
  a negative interaction while an entry of zero is inconclusive.}
}

\value{
  A matrix.

  The return value is an kx2 matrix where k is the length of
  matList. The first colum returns the estimated false positive error
  rates and the second column returns the estimated false negative error
  rates.
}



\details{The model is described in the manuscript
\emph{Estimating node degree in bait-prey graphs.}
by D. Scholtens et al. 
}

\references{Scholtens D, Chiang T, Huber W, Gentleman R.  Estimating node
degree in bait-prey graphs. \emph{Bioinformatics}.}

\author{T. Chiang and L. Wang}

\examples{
load(system.file("extdata", "intacty2hppMatrix.rda", package="ppiStats"))
load(system.file("extdata", "PPIpos", package="ppiStats"))
x = intacty2hppMatrix[c("EBI-389903","EBI-783101")]
estimatePPIErrorRates(x, GSPos=PPIpos)
}


\keyword{manip}
