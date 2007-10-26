\name{estimateCCMErrorRates}
\alias{estimateCCMErrorRates}
\title{Estimate false positive and false negative error probabilities
  }
\description{Estimate false positive and false negative error probabilities
  for complex comembership edges using a protein complex interactome
  gold standard}
\usage{
estimateCCMErrorRates(m,GS,filterSystematic=TRUE,	
         obsPropThresh=.7,SystematicpThresh=.01,
         nBS = 100, seed = 123)
}
\arguments{
  \item{m}{The bait to prey data adjacency matrix. Baits index the rows and prey index
           the column.}
  \item{GS}{A gold standard protein complex interaction incidence matrix. Proteins index
            the rows and protein complex index the columns. }
  \item{filterSystematic}{A logical. If true, all proteins with with highly uneven
                          directed degree will be filtered out of the data.}
  \item{obsPropThresh}{A numeric between 0 and 1. The proportion of tested proteins found
                       within a protein complex needed to keep that protein complex within
        	       the gold standard set.}
  \item{SystematicpThresh}{A numeric between 0 and 1. The p-value threshold by which
                           systematic errors are filtered.}
  \item{nBS}{A positive integer. The number of times needed to bootstrap.}
  \item{seed}{A positive integer. The seed for the random number generator. }
}

\value{
 A list:
}

\details{The model is described in the manuscript
\emph{Estimating node degree in bait-prey graphs.}
by D. Scholtens et al. 
}

\references{Scholtens D, Chiang T, Huber W, Gentleman R.  Estimating node
degree in bait-prey graphs. (Submitted).}

\author{T. Chiang and D. Scholtens}

\examples{

}


\keyword{manip}
