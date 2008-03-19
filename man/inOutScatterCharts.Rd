\name{inOutScatterCharts}
\alias{inOutScatterCharts}

\title{A function to create scatter plots of the in-degree and out-degree
	 for each vertex of a directed graph}
\description{
  This function takes a list of directed graph objects and plots the out-degree
  against the in-degree for each vertex. A binomial test is used to determine	
  if there is a highly disproportionate in-degree or out-degree given a p-value 
  threshold (the binomial test presumes p = 0.5).  
}
\usage{
inOutScatterCharts(dataGraphs, pThresh=0.01, pLevels=1e-4)
}

\arguments{
  \item{dataGraphs}{A named list of directed graphNELs}
  \item{pThresh}{The two-sided p-value threshold for the binomial testing
  for potentially biased nodes}
  \item{pLevels}{Figure out what this does}
}

\value{

	For each directed data graph, this function generates 3 pdf scatter-plots and 1 eps 
	scatter plot. The pdf files are:
	1. A scatter plot of each nodes out- vs in-degree where points outside the staircase
	    bundaries are those rejected in the binomial test.
	2. A scatter plot of each nodes out- vs in-degree scaled by sqrt where points outside 
	    the conic region are those rejected in the binomial test.    
	3. A histogram for the distribution of p-values for each node

	The eps file contains the same information as the second pdf file.	    

}

\author{T Chiang}


\examples{
graphs = lapply(bpExperimentNames, function(x) get(x))
names(graphs) = bpExperimentNames
#inOutScatterCharts(graphs)
}
\keyword{datagen}
