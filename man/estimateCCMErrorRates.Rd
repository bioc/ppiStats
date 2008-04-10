\name{estimateCCMErrorRates}
\alias{estimateCCMErrorRates}
\title{Estimate false positive and false negative error probabilities}
\description{Estimate false positive and false negative error probabilities
  for complex comembership edges using a protein complex interactome
  gold standard}
\usage{
estimateCCMErrorRates(m,GS,filterSystematic=TRUE,	
         obsPropThresh=1,SystematicpThresh=.01)
}
\arguments{
  \item{m}{The bait to prey data adjacency matrix. Baits index the rows and prey index
           the columns.}
  \item{GS}{A gold standard protein complex interaction incidence matrix. Proteins index
            the rows and protein complexes index the columns. }
  \item{filterSystematic}{A logical. If TRUE, all baits with with highly uneven
                          directed degree will be filtered out of the data.}
  \item{obsPropThresh}{A numeric between 0 and 1. The proportion of tested proteins found
                       within a protein complex needed to keep that protein complex within
        	       the gold standard set.}
  \item{SystematicpThresh}{A numeric between 0 and 1. The p-value threshold by which
                           systematic errors are filtered.}
}

\value{
 A list:


  \item{globalpTP}{A numeric between 0 and 1.  Estimate of pTP.}
  \item{globalpTPSE}{A numeric.  Estimate of standard error of globalpTP estimate.}
  \item{globalpFP}{A numeric between 0 and 1.  Estimate of pFP.}
  \item{pTP95CI}{A vector of length 2.  95 percent confidence interval upper and lower
  bounds for globalpTP estimate.}
  \item{pFP95CI}{A vector of length 2.  95 percent confidence interval upper and lower
  bounds for globalpFP estimate.}
  \item{nEligComplexes}{A numeric.  Number of complexes from GS that met
  obsPropThresh criteria.}
  \item{nEligBaits}{A numeric. Total number of eligible baits in GS set.}
  \item{nEligEdges}{A numeric.  Total number of eligible edges in GS set.}
  \item{nBaitsInComplexes}{A vector.  Number of baits in each eligible complex.}
  \item{complexSizes}{A vector. Size of each complex in GS set.}

}

\details{The model is described in the manuscript
\emph{Estimating node degree in bait-prey graphs.}
by D. Scholtens et al. 
}

\references{Scholtens D, Chiang T, Huber W, Gentleman R.  Estimating node
degree in bait-prey graphs. \emph{Bioinformatics}.  To appear.}

\author{T. Chiang and D. Scholtens}

\examples{
data(Ho2002BPGraph)
data(ScISIC)
Ho2002mat = as(Ho2002BPGraph,"matrix")
estimateCCMErrorRates(Ho2002mat,ScISIC)$globalpTP
}


\keyword{manip}
