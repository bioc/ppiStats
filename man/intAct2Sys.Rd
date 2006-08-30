\name{intAct2Sys}
\alias{intAct2Sys}

\title{A function that maps the Intact ID's to the first YEAST
  systematic name using the YEAST package}
\description{
This function takse the Intact ID's maps to the yeast common names via intact
repository and uses the YEAST package to map the
common names to the YEAST Systematic names. There is problems if more
than one systematic name corresponds to the common name. This function
will simply always return the first name given by each mapping: from
Intact to common names and also from common names to systematic names.
}
\usage{
intAct2Sys(prot2Sys, bpSysL)
}

\arguments{
  \item{prot2Sys}{An Intact accession number for some protein which is
    needed to be mapped to its systematic name}
  \item{bpSysL}{A named list with the IntAct accession numbers as the names and the
    corresponding systematic names as the entries}
}

\value{
  The return value is a mapping of each Intact ID to it's respecting
  systematic name.

}

\author{T Chiang}

\examples{

}
\keyword{datagen}
