\name{estErrProbMethodOfMoments}
\alias{estErrProbMethodOfMoments}
\title{Estimate false positive and false negative error probabilities
  by method moments.}
\description{Estimate false positive and false negative error probabilities
  by method moments.}
\usage{
estErrProbMethodOfMoments(nint, nrec, nunr, ntot)
}
\arguments{
  \item{nint}{Integer vector. True number of interactions. Typically,
    the function is called for a range of these, returning all possible
    solutions for that range.}
  \item{nrec}{Integer scalar. Observed number of reciprocated edges.}
  \item{nunr}{Integer scalar. Observed number of unreciprocated edges.}
  \item{ntot}{Integer scalar. Number of proteins which were tested twice
    (e.g. both as viable bait and as viable prey).}
}

\value{
  Matrix with 5 columns \code{nint} (a copy of the input argument),
  \code{pfp1}, \code{pfn1}, \code{pfp2} and \code{pfn2}, and as many
  rows as the length of \code{nint}.
}

\details{The model is described in the vignette
\emph{Stochastic and systematic errors in PPI data, by looking
at unreciprocated in- or out-edges}
by W. Huber, T. Chiang and R. Gentleman. 
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{

est = estErrProbMethodOfMoments(nint=seq(8000, 40000, by=100), nrec=9722, nunr=15856, ntot=2000)
if(interactive()) {
  plot(est[, c("pfp2", "pfn2")], type="l", col="blue", lwd=2,
       xlab=expression(p[FP]), ylab=expression(p[FN]))
  abline(h=0, v=0, lty=2)
}
}

\keyword{manip}
