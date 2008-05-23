\name{findAdjacent}
\alias{findAdjacent}


\title{A function that returns the set of prey detected by a bait in a
bait-prey experiment.}
\description{
  This function returns a list of the unique prey detected by a bait in a
  bait-prey experiment.
}
\usage{
findAdjacent(x,bait)
}

\arguments{
  \item{x}{A graph object containing the set of bait-prey interaction data.}
  \item{bait}{A character specifying tjhe bait node of interest.}
}

\value{
  A vector of prey that are detected by the bait in the graph, or NULL if the
  supplied bait is not a node in the graph.
}

\author{Denise Scholtens}

\examples{

}
\keyword{graphs}
