\name{buildGOhtml}
\alias{buildGOhtml}

\title{A function to generate a .html file for the
       output of the GO hypergeometric tests.}
\description{
This file takes in a list of the outputs from
the goHyperGeo function and produces a .html
file showing the GO terms which are either 
over or under represented.
}
\usage{
buildGOhtml(hgTestOutPut, direc = "Over", Onto = "CC")
}

\arguments{
  \item{hgTestOutPut}{A name list. Each entry in the list is 
	an output from the function goHyperGeoTest of the 
	GOstats package.}
  \item{direc}{A character. Tells whether the geometric test
	is testing over or under representation.}
  \item{Onto}{A character. Describes which ontology was used
	in the hyper-geometric test.}
}

\value{
An .html file is produced in your working directory.
}

\author{ T Chiang }

\examples{
	
}
\keyword{datagen}
