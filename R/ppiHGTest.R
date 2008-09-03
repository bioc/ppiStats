ppiBuildParams4GO <- function(geneSet, universe, direction="over",
                           annot = "org.Sc.sgd", ontology = "CC", cond=TRUE, pThresh = 0.01){
  parameter <- new("GOHyperGParams", geneIds = geneSet,
                   universeGeneIds = universe, annotation=annot,
                   ontology = ontology, conditional=cond,
                   testDirection = direction, pvalueCutoff = pThresh)

  return(parameter)
}

ppiHGTest4GO <- function(parameter, filename, append=TRUE,
                      label = "Experiment name here",
                      typeGeneSet = "Describe the gene set here",
                      cs=50){

  print(paste("Processing", filename, "for", parameter@testDirection,
              "representation in the ", parameter@ontology, "GO ontology."))
  results <- hyperGTest(parameter)
  keep <- any(pvalues(results)<results@pvalueCutoff)
  htmlReport(results, file = paste(filename, ".html", sep=""),
             append = append, label = paste(label,": ", typeGeneSet, sep=""))
  return(results)
  
}

ppiBuildParams4PFAM <- function(geneSet, universe, annot ="org.Sc.sgd",
                              direction = "over", pThresh=0.01){
  parameter <- new("PFAMHyperGParams", geneIds = geneSet, universeGeneIds =
                   universe, annotation = annot, testDirection = direction,
                   pvalueCutoff = pThresh)

  return(parameter)
}

ppiHGTest4PFAM <- function(parameter, filename, append = TRUE,
                           label = "Experiment Name Here",
                           typeGeneSet = "Describe the Gene Set Here",
                           cs = 50){
  print(paste("Processing", filename, "for", parameter@testDirection,
              "representation of PFAM."))
  results <- hyperGTest(parameter)
  keep <- any(pvalues(results)<results@pvalueCutoff)
  htmlReport(results, file = paste(filename, ".html", sep=""),
             append = append, label = paste(label,": ", typeGeneSet, sep=""))
  return(results)
}
