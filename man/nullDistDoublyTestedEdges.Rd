\name{nullDistDoublyTestedEdges}
\alias{nullDistDoublyTestedEdges}
\title{Null distribution of number of reciprocated,
    unreciprocated and missing edges in stochastic model.}
\description{
  Calculate the null distribution of the number of reciprocated,
  unreciprocated and missing edges in a stochastic model where each
  edge is tested twice.
}
\usage{
nullDistDoublyTestedEdges(deltaMax, n, pFP, pFN)
}
\arguments{
  \item{deltaMax}{Integer. Distributions will be calculated
   for model parameter \emph{delta}\code{=0, 1, 2, ..., deltaMax}.}
  \item{n}{Integer. The parameter \emph{n} of the model.}
  \item{pFP}{Numeric. The parameter \emph{pFP} of the model.}
  \item{pFN}{Numeric. The parameter \emph{pFN} of the model.}
}

\value{
3d array with dimensions \code{nMax+1} x \code{nMax+1} x
\code{deltaMax+1} whose element \code{p[nr+1, nu+1, delta+1]}
is the corresponding joint probability.  \code{nMax+1} is calculated
(probably too conservatively) by the function to make sure that no
probability leaks out of the array.
}

\details{The model is described in the vignette
\emph{Stochastic and systematic errors in PPI data, by looking
at unreciprocated in- or out-edges}
by W. Huber, T. Chiang and R. Gentleman. 

This function can be quite slow, its runtime grows quickly
with \code{deltaMax} (and is roughly independent of \code{n}, \code{pFP},
\code{pFN}). The example below should take only a few seconds on a reasonable
computer, though.
}

\author{Wolfgang Huber \url{http://www.ebi.ac.uk/huber}}

\examples{
p = nullDistDoublyTestedEdges(32, 1000, pFP=0.001, pFN=0.15)

if(interactive() && require("RColorBrewer"))
  for(k in 1:dim(p)[3]) {
    image(sqrt(p[,,k]), xlab=expression(N[rec]), ylab=expression(N[unrec]),
        main = expression(P(N[rec], N[unrec]~";"~ delta^"*", n, p[FP], p[FN])),
        x = 1:dim(p)[1], y = 1:dim(p)[2], 
        col = colorRampPalette(brewer.pal(9, "GnBu"))(256))
    text(35, 35, paste("delta", k, sep="="))
  }
}

\keyword{manip}
