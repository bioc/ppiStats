\name{hgParams}
\alias{ppiBuildParams4GO}
\alias{ppiBuildParams4PFAM}

\title{A wrapper function to build a parameter class for the 
input of the HyperGTest.}
\description{
  This function takes a gene set and conducts test for either over or under
  representation of some category using the Hypergeometric distribution.

  The two differences when building the parameter classes are that a conditional
  test can be performed on the GO dag but not on PFAM categories and an ontology
  can be assigned to GO but not to PFAM.
}
\usage{
ppiBuildParams4GO(geneSet, universe, direction="over", annot="org.Sc.sgd",
                           ontology = "CC", cond=TRUE, pThresh = 0.01)
ppiBuildParams4PFAM(geneSet, universe, annot ="org.Sc.sgd",
                              direction = "over", pThresh=0.01)
}

\arguments{
  \item{geneSet}{A character vector of genes given by the gene locus name.}
  \item{universe}{The set of genes by which the geneSet is tested against 
                            for over/under representation. The genes are also given 
                            by the gene locus names.}
  \item{direction}{A character. This parameter can be either set to over or under
                            when testing for GO categories.}
  \item{annot}{A character. The annotation package used.}
  \item{ontology}{A character: either CC, MF, or BP to describe the GO ontology.}
  \item{cond}{A logical. To test within the CO dag, a conditional hypergeometric test
                      can be conducted.}
  \item{pThresh}{A numeric. A p-value threshold by which the null hypothesis is rejected.}
}

\value{
  A object of class hyperGParams.
}
\author{T Chiang}

\keyword{datagen}
