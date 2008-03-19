\name{viabilityCharts}
\alias{viabilityCharts}

\title{A function to create summary bar charts for directed graphs}
\description{
  This function takes a list of directed graph objects and creates a summary of the 
  relative viable baits, viable prey, and viable bait/prey of each data graph.
}
\usage{
viabilityCharts(dataGraphs, total=6466)
}

\arguments{
  \item{dataGraphs}{A named list of directed graphNELs}
  \item{total}{The total number of potential nodes (i.e. proteins) that could have been queried.}
}

\value{

	A barchart object. Each item in the bar-chart represents one experimental
	data graph (set). 

}

\author{T Chiang}


\examples{
graphs = lapply(bpExperimentNames, function(x) get(x))
names(graphs) = bpExperimentNames
viabilityCharts(graphs)
}
\keyword{datagen}
