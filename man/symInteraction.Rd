\name{symInteraction}
\alias{symInteraction}

\title{A function that returns symmetric bait to prey interactions from
  the y2h experiments.}
\description{
  This function parses through the argument y2h and returns only those
  bait to prey interactions for which the bait b detected prey p and
  vice versa. 
}
\usage{
symInteraction(y2h)
}

\arguments{
  \item{y2h}{A three tiered list. Names of the second tier corresponds to
  ID of experiment. Names of the third tier should be the bait
  proteins. The entries of the third tier are character vectors of
  proteins detected as prey}
}

\value{
  The return value is identical in structure to the argument. It is
  still a three tiered list. Names of the second tier corresponds to
  ID of experiment. Names of the third tier should be the bait
  proteins, but this is now restricted so that these baits must be in
  symmetric interactions with respect to its prey set. The entries of the
  third tier are 
  character vectors of proteins detected as prey, but this is also
  restricted to those prey which have participated in symmetric
  interaction relative the the bait which indexes them. . 
}

\author{T Chiang}
\examples{
data(ppiData)
data(y2hSysGW)
eg = y2hSysGW[3]
symInteraction(eg)
}
\keyword{datagen}
