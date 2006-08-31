\name{y2h}
\alias{y2h}
\docType{data}
\title{A data set documenting the bait to prey associations for the
  y2h empirical data-sets as found within the IntAct repository.}

\usage{data(y2h)}
\format{
  The format is:
 chr "y2h"
}
\description{
  This data set is the bait to prey assocation after the function
  createBPList has been called on three entries of intActPPIData:
  indexSetAll, baitsSystematic, preysSystematic.

  The structure of this R data object is a three-tiered list. The 
  first (or top) tier of the list has 42 sub-lists corresponding to 
  the 42 y2h data-sets found within IntAct. Each sub-list is named 
  by the corresponding experimentor and the date of publication. 
  Each one of these sub-lists contains a set of lists; each of these
  lists is named by an individual viable bait of that experiment. 
  Each viable bait list contains a character vector of the viable
  preys dectect by the corresponding viable bait.

  For a more in depth summary, please see the intactPPIData.Rnw file.
}
\examples{
data(y2h)

}
\keyword{datasets}
