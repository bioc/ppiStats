\name{twPERM}
\alias{twPERM}
\alias{twowayPERM}
\alias{makeBinVect}


\title{A function to compute a permutation test on a two way table.h}
\description{
  This function takes the entries of a two-way table and creates
  two binary logical vectors upon which it computes a permutation
  test to check for independence.
}
\usage{
twowayPERM( v1, v2, NPERM, stat, seed = 123)
makeBinVect(n11, n12, n21, n22)
}

\arguments{
  \item{v1}{A logical vector.}
  \item{v2}{A logical vector the same length as v1.}
  \item{NPERM}{A positive integer. The number of tests to be conducted.}
  \item{stat}{A statistical function passed into to test the
	permutation on two vectors.}
  \item{seed}{A positive integer. To set the seed for the random number
	generator.}
  \item{n11}{A postive integer. The value of (1,1) in the two way table.}
  \item{n12}{A postive integer. The value of (1,2) in the two way table.}	
  \item{n21}{A postive integer. The value of (2,1) in the two way table.}	
  \item{n22}{A postive integer. The value of (2,2) in the two way table.}	
}	

\value{
	A numeric vector containing the p-values for the test that
	the permutation of the vector v1 is independent from v2.
}

\author{T Chiang}

\examples{
x <- makeBinVect(13, 17, 23, 71)
}
\keyword{datagen}
