\name{createBPList}
\alias{createBPList}
\title{A function to create the Bait to Prey association list.}
\description{
  This function takes the indexSetAll, baitsSystematic, and
  preysSystematic entries respectively from the output list of the R
  function collectIntactPPIData, and it generates a three tier-ed list
  which ultimately records the bait to prey association from each
  protein-protein interaction (ppi) experiment obtained from the IntAct 
  repository.

  
}
\usage{
createBPList(indexSet, baitsSystematic, preysSystematic)
}

\arguments{
  \item{indexSet}{A list of named lists. The names correspond to each
    experiment while the entries to the lower level lists are character
    vectors with two elements: the first element is the id for the bait
    protein and the second an id for the prey. This is an entry of the
    output of the collectIntactData function called indexSetAll}
  \item{baitsSystematic}{A named list. The names correspond to the
    IntAct acension codes for any proteins used as baits, and the
    entries of the list are character vectors consisting of Systematic
    gene names with respect to the IntAct ID. This is also an entry from
    the output of the function collectIntactData.}
  \item{preysSystematic}{A named list. The names correspond to the
    IntAct acension codes for any proteins used as preys, and the
    entries of the list are character vectors consisting of Systematic
    gene names with respect to the IntAct ID. This is also an entry from
    the output of the function collectIntactData.}
  
}

\value{
  This function generates a three tiered list. The top tier is the main
  list. The second tier-ed set of lists
  should be named or referenced by experimentors. The third tier is a
  list of bait to prey associations. The names in the third tier of the
  list are those proteins sampled as baits, and the entries are
  character vectors of proteins detected as prey. 
}

\author{T Chiang}

\examples{
data(y2h)
eg = y2h$indexSetAll[3]
eg1 = y2h$baitsSystematic[3]
eg2 = y2h$preysSystematic[3]
createBPList(eg, eg1, eg2)
}
\keyword{datagen}
