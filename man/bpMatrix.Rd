\name{bpMatrix}
\alias{bpMatrix}

\title{This function generates the Bait-Prey Matrix for an
	    protein-protein interaction (ppi) experiment}
\description{
  This function takes in experimental ppi data and generates the bait to
  prey adjacency matrix. 
}
\usage{
bpMatrix(y2h, symMat = TRUE, homodimer = FALSE, baitAsPrey = FALSE, 
	      unWeighted=TRUE, onlyRecip=FALSE, baitsOnly=FALSE)
}

\arguments{
	  \item{y2h}{A named list of character vectors. 
  	  The names of the list correspond to the baits
  	  used in the experimental technology. The 
	  entries of the character vectors are those 
	  proteins found as prey or a character(0) if 
  	  the bait did not find any prey.}
  	  \item{symMat}{A logical, if TRUE, the matrix will be square with all
    	  the proteins documented in the experiment indexing both the row and
    	  column; if FALSE, only the baits index the rows, preys the columns.}
  	  \item{homodimer}{A logical. If TRUE, the matrix will record the
    	  presence of homodimers; if FALSE, all homodimers data will be
    	  deleted.}
  	  \item{baitAsPrey}{A logical; if TRUE, the columns will be indexed by
    	  both the bait and prey population while the rows will remain indexed
    	  by the baits exclusively.}
  	  \item{unWeighted}{A logical. If TRUE, the entries of the adjacency will be 
    	  binary (0,1) which records the presence of interactions or not. If FALSE,
    	  the entries of the matrix will be a natrual number to record the 
    	  multiplicity of the interaction found by the experiment.}
  	  \item{onlyRecip}{A logical. If TRUE, the adjacency matrix will be restricted 
    	  to only those interactions which are symmetric.}
  	  \item{baitsOnly}{A logical. If TRUE, the matrix will be indexed by the baits
    	  for both the rows and columns. If baitsOnly is TRUE, then baitsAsPrey must
	  also be TRUE and symMat must be FALSE}
}
\details{
  It is important to note that the weight of each directed edge is
  recorded by the number of instances bait b finds prey p.
}
\value{
  An adjacency matrix with weighting. The rows are indexed by those
  proteins sampled as baits (if symMat is true, only those rows with
  non-trivial row sums were sampled as baits) and the colunms are
  indexed by proteins detected as baits (et cetera).
}


\author{T Chiang}

\examples{
library(ppiData)
data(y2hSysGW)
eg <- y2hSysGW[2:4]
lapply(eg, bpMatrix)
}
\keyword{datagen}
