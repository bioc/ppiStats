\name{hgTests}
\alias{ppiHGTest4GO}
\alias{ppiHGTest4PFAM}

\title{A wrapper function to implement the Hypergeometric test,  
HyperGTest found withing the Category and GOstats packages.}
\description{
  This function takes the instances of the hyperGParams and conducts
  a (conditional) test for over/under representation of some category
  under the Hypergeometric distribution.
}
\usage{
ppiHGTest4GO(parameter, filename, append=TRUE,
                      label = "Experiment name here",
                      typeGeneSet = "Describe the gene set here",
                      cs=50)

ppiHGTest4PFAM(parameter, filename, append = TRUE,
                           label = "Experiment Name Here",
                           typeGeneSet = "Describe the Gene Set Here",
                           cs = 50)
}

\arguments{
  \item{parameter}{An object of hyperGParams.}
  \item{filename}{A character vector. The name given to the .html file produced.}
  \item{append}{A logical. If multiple tests are conducted using the same filename, then the 
                          .html file will be appended if TRUE or over-written if FALSE.}
  \item{label}{A character. The should give a description of the experiment used to 
                      obtain the gene set.}
  \item{typeGeneSet}{A character: this character vector should adequately describe the gene
                                 set itself.}
  \item{cs}{A numeric. This gives a cut-off for the size of the categories in which we 
                  conduct the test.}
 }

\value{
  A instance of the class hyperGTest as well as an .html file.
}
\author{T Chiang}

\keyword{datagen}
