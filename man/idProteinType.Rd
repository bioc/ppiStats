\name{idProteinType}
\alias{idViableProteins}
\alias{idHomodimers}

\title{A function to determine viable baits, viable preys, or
  homodimers within experimental data-sets.}
\description{
  These functions take a bait to prey directed graphNEL and returns
  either a character vector of all proteins which participates in
  homodimer relationships or a list of three character vectors: a vector
  of viable baits, a vector of viable prey, and a vector of the viable
  bait/prey.
}
\usage{
idViableProteins(bpGraph, homomer=TRUE)
idHomodimers(bpGraph)
}

\arguments{
  \item{bpGraph}{A directed graphNEL}
  \item{homomer}{A logical. If True, homomer relationships will also be used to 
                 characterize proteins as viable baits, viable prey, and also as
		 viable bait/prey simultaneously. If False, those proteins whose
		 only intraction is a homomer relationship will not be 
		 characterized as viable in any context.}
}

\value{

  The return value for idHomodimers is a character vector of those
  proteins which participates in homomer relationships.

  The return value for idViableProteins is a list of two character
  vectors:
  
  \item{VB}{A vector of baits that finds at least one prey in
            the experimental graphNEL, i.e. the nodes with out-degree at greater
    	    than 0. If homomer = FALSE, those proteins which participate in 
   	    only homomer relationships will not be returned.}
  \item{VP}{A vector of prey which is found by at least one bait in
            the experimental graphNEL, i.e. the nodes with in-degree at greater
    	    than 0. If homomer = FALSE, those proteins which participate in 
    	    only homomer relationships will not be returned.}
  \item{VBP}{A vector of proteins that are both viable baits and viable
             prey in an experimental dataset as defined above. If 
	     homomer = FALSE, those proteins which participate in 
	     only homomer relationships will not be returned.}	
 }

\author{T Chiang}
\examples{
library(ppiData)
idViableProteins(ItoCore2001BPGraph)
idHomodimers(ItoCore2001BPGraph)
}
\keyword{datagen}
