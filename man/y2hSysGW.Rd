\name{y2hSysGW}
\alias{y2hSysGW}
\docType{data}
\title{A data object documenting the bait to prey associations for 
	the y2h empirical data-sets with two restrictions.}
\description{
  This list is a sublist of the file y2h.rda of bait to prey 
  interactions which have 2 restrictions: 
  Genome Wide prey population and the use of GAL4 transcription factor.
  There are currently 7 experiments that fall within these restrictions.
  The protein names in this list have been mapped to the systematic
  names from the IntAct codes.
}
\usage{data(y2hSysGW)}
\format{
  The format is:
 chr "y2hSysGW"
}
\details{
  The format of this data object is the same as that of y2h.rda. Please 
  see the help page for that object.
  
}

\examples{
data(y2hSysGW)
y2hSysGW[[2]]
}
\keyword{datasets}
