\name{map2Systematic2}
\alias{map2Systematic2}

\title{A secondary function that maps the Intact ID's to Yeast Systematic Names if
possible}
\description{
  This function takse the Intact accession numbers and  maps to the SGD ID's via intact
  repository and truncates the SGD ID's to find the yeast systematic
  names. 
}
\usage{
map2Systematic2(allProt, tableList, sWAC)
}
\arguments{
  \item{allProt}{A character vector of protein ID's}
  \item{tableList}{A List of tables. The parsed XLM files of the IntAct
    database}
  \item{sWAC}{The data file sWAC2Sys}
}
}
\value{
  A character vector - The return value is a mapping of each Intact ID
  to it's respecting systematic name.
}

\author{T Chiang}


\keyword{datagen}
