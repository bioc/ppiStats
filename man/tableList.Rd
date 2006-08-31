\name{tableList}
\alias{tableList}
\docType{data}
\title{A list of tables for the downloaded and parsed xml files from the
IntAct database}
\usage{data(tableList)}
\format{
  The format is:
 chr "tableList"
}
\description{
  This is a list of tables generated from the a script used to parse the
  xml files from the IntAct data repository.

  tableList is a list that contains six data-frames (tables) that records
  all the various data given by IntAct. The names of these six tables
  are:

  1. "acInfo"
  2. "experimentInfo"
  3. "ac2xref"
  4. "interactionInfo"
  5. "experiment2interaction"
  6. "interaction2inteactor"

  The "acInfo" data frame contains 4 columns: "ac", "type",
  "shortLabel", "fullName". They are IntAct accession numbers, types
  (experiment, interaction, or interactor), and short/long names,
  respectively.  The reason to have the second column is because IntAct
  assigns an accession number to each experiment, interaction, and
  interactor. The second column is an indicator between these three. The
  accession numbers in the first column of this data frame are those
  accession numbers in the yeast xml files from the IntAct repository.

  There are also 4 columns in the "experimentInfo" data frame:
  "ac","hostOrganism","interactionDetection", and
  "participantDetection". The IntAct accession numbers of experiments
  will be in the first column. The "hostOrganism" column is the tax id
  of the host organism (yeast in this case).  The next two columns are
  the detection methods of the interaction and of the participant,
  respectively.  The values in the 3rd and 4th columns are Proteomics
  Standards Initiative (PSI) Molecular Interaction (MI) identification
  codes. IntAct documents both the PSI MI 1.0 as well as PSI MI 2.5
  codes.

  The "ac2xref" data frame contains 4 columns:
  "ac", "db", "id", and "secondary".
  This data frame provides the cross references (among several databases)
  for each IntAct accession number.  The column "ac" is the IntAct
  accession numbers, and the columns "db" and "id" are the reference
  database names and the corresponding identification code for that
  database. The secondary column provides secondary id if provided.
  
  There are also 4 columns in the "interactionInfo" data frame:
  "ac", "interactionType", "confidenceUnit", and "confidenceValue".
  The first column contains all the interaction accession numbers among
  these xml files.  The second column is the type of interaction.  The
  values here are also the PSI MI codes. The third and forth columns are
  the confidence information for this interaction.

  There are 2 columns in the "experiment2interaction" data frame:
  "experiment", and "interaction". This data frame tells us the
  what types of interactions were derived from the experiment.
  Both the experiments and the interactions are given by the IntAct
  accession numbers.

  There are 3 columns in the "interaction2inteactor" data frame:
  "interaction", "interactor", and "role". The first column provide the
  type of interaction while the second column provides the interactors
  within the interaction. Both the "interaction" and the "interactor"
  are given by the IntAct accession numbers. The third column is the
  role (bait or prey) of the interactors.
  
}
\examples{
data(tableList)

}
\keyword{datasets}
