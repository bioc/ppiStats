\name{map2Systematic}
\alias{map2Systematic}

\title{A function that maps the Intact ID's to Yeast Systematic Names if
possible}
\description{
This function takse the IntAct accension numbers and maps to the yeast common names via intact
repository and uses the YEAST package to map the
common names to the YEAST Systematic names.
}
\usage{
map2Systematic(allProt, tableList, sWAC)
}

\arguments{
  \item{allProt}{A character vector of proteins to be mapped}
  \item{tableList}{The tableList element of the output from collectIntactData}
  \item{sWAC}{A named character vector; the names are SwissProt Accesion
    Codes and the elemnets are the Yeast Systematic names}
}

\value{
  The return value is a mapping of each Intact ID to it's respecting
  systematic name.

}

\author{T Chiang}

\examples{
data(sWAC2Sys)
dataL = collectIntactData("EBI-375746")
sysN = map2Systematic(dataL$allBaits, dataL$tableList, sWAC2Sys)
}
\keyword{datagen}
