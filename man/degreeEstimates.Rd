\name{degreeEstimates}
\alias{degreeEstimates}
\alias{degreePMF}
\alias{findDegree}
\title{Estimate per-node degree in viable bait-prey (VBP) graph using maximum
likelihood.}

\description{Estimate per-node degree in viable bait-prey (VBP) graph using
maximum likelihood.}
\usage{
degreeEstimates(m,pTP,pFP)
findDegree(rd,ud,pTP,pFP,N)
degreePMF(deltahat,rd,ud,pTP,pFP,N)
}
\arguments{
  \item{m}{Square incidence matrix for VBP graph.}
  \item{pTP}{True positive probability.}
  \item{pFP}{False positive probability.}
  \item{rd}{Observed number of reciprocated edges.}
  \item{ud}{Observed number of unreciprocated edges.}
  \item{deltahat}{Estimate of node degree.}
  \item{N}{Number of proteins which were tested twice
    (e.g. both as viable bait and as viable prey.)  Should equal the number of
    rows of m.}
}

\value{
  \code{degreeEstimates} returns a named vector of degree estimates for each
  node in the graph.

  \code{findDegree} returns a single degree estimate for one node.

  \code{degreePMF} returns the value of the pmf for the multinomial model at a
  specified estimate of node degree, given the number of observed reciprocated
  and unreciprocated edges incident on the node and pTP and pFP.

}

\details{\code{degreeEstimates} returns per-node degree estimates using the
maximum likelihood method of Scholtens et al. (Submitted).  It takes arguments
\code{m}, which is an incidence matrix of bait-prey relationships, typically a
VBP graph filtered for proteins prone to systematic bias, as well as pTP and
pFP values, globally applicable to the entire graph.

\code{degreeEstimates} calls the function \code{findDegree} which estimates
degree for a single node, given its observed number of reciprocated and
unreciprocated incident edges.  \code{findDegree} takes an argument \code{N}
which is the number of proteins in the graph that were tested twice.  When
\code{degreeEstimates} calls \code{findDegree}, \code{N} is set to the first
dimension of the incidence matrix \code{m}.  

\code{degreePMF} calculates the value of the pmf for an estimated degree,
given observed numbers of reciprocated and unreciprocated edges, as well as
pTP, pFP, and N.  It is not generally called directly by the user.  It is used
to locate the maximum likelihood estimator for degree.

}

\references{Scholtens D, Chiang T, Huber W, Gentleman R.  Estimating node
degree in bait-prey graphs. (Submitted).}

\author{Denise Scholtens, dscholtens@northwestern.edu}

\examples{

findDegree(rd=2,ud=2,pTP=0.75,pFP=0.001,N=1000)

}

\keyword{graphs}
