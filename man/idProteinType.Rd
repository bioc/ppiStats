\name{idProteinType}
\alias{idViableProteins}
\alias{idHomodimers}

\title{A function to determine viable baits, viable preys, or
  homodimers within experimental data-sets.}
\description{
  These functions take a bait to prey directed graphNEL and returns
  either a character vector of all proteins which participates in
  homodimer relationships or a list of two character vectors: a vector
  of viable baits as well as a vector of viable prey.
}
\usage{
idViableProteins(bpGraph)
idHomodimers(bpGraph)
}

\arguments{
  \item{bpGraph}{A direced graphNEL}
}

\value{

  The return value for idHomodimers is a character vector of those
  proteins which participates in homodimer relationships.

  The return value for idViableProteins is a list of two character
  vectors:
  
  \item{viableBaits}{A vector of baits that finds at least one prey in
    the experimental graphNEL, i.e. the nodes with out-degree at greater
    than 0}
  \item{viablePrey}{A vector of prey which is found by at least one bait in
    the experimental graphNEL, i.e. the nodes with in-degree at greater
    than 0}
 }

\author{T Chiang}
\examples{
idViableProteins(Ito2001BPGraph)
idHomodimers(Ito2001BPGraph)
}
\keyword{datagen}
