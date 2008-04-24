\name{separateExptBySize}
\alias{separateExptBySize}
\title{A function to partition experiments in user defined small, medium, large
scale experiments}
\description{This function partitions protein interaction data into one of three classes
defined by boundaries supplied by the users.}
\usage{
separateExptBySize(listOfGraphs, bound1, bound2)
}
\arguments{
  \item{listOfGraphs}{A list containing graphNEL objects to be sorted by the number of edges}
  \item{bound1}{An integer. The strict upper bound for the small scale experiments and the
                         lower bound for the medium scale experiments.}
  \item{bound2}{An integer. The upper bound for the medium scale experiments, and the 
                          strict lower bound for the large scale experiments.}
}
\value{
 A list:
  \item{small}{A list of graph objects with the number of edges in each graph less than bound1.}
  \item{medium}{A list of graph objects with the number of edges in each graph of a value between
                           that of bound1 and bound2.}
  \item{large}{A list of graph objects with the number of edges in each graph greater than bound2.}
}

\details{
A function to partition graphs based on the number of edges.
}

\author{T. Chiang }

\keyword{manip}
