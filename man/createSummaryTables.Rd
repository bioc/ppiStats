\name{createSummaryTables}
\alias{createSummaryTables}

\title{A function to create summary statistic tables for directed graphs}
\description{
  This function takes a list of directed graph objects and creates a summary table
  based on the directed connectivity.
}
\usage{
createSummaryTables(dataGraphs)
}

\arguments{
  \item{dataGraphs}{A named list of directed graphNELs}
}

\value{

	A dataframe with the rows indexed by the names of each directed graphNEL 
	and the columns indexed by viable baits (VB), viable prey (VP), viable bait/prey
	(VBP), the ratio of viable bait/prey by viable baits (VBP/VB), the ratio of viable
	prey by viable baits (VP/VB), the total number of directed interactions (TI), and
	the ratio of total intractions by viable baits (TI/VB).

}

\author{T Chiang}


\examples{
data("bpExperimentNames", package = "ppiData")
graphs <- new.env(parent = emptyenv())
data(list = bpExperimentNames, package = "ppiData", envir = graphs)
graphs <- lapply(bpExperimentNames, function(x) graphs[[x]])
names(graphs) <- bpExperimentNames
createSummaryTables(graphs)
}
\keyword{datagen}
