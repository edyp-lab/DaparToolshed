% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agregation.R
\name{GetDetailedNbPeptidesUsed}
\alias{GetDetailedNbPeptidesUsed}
\title{Computes the detailed number of peptides used for aggregating
each protein}
\usage{
GetDetailedNbPeptidesUsed(pepData, X)
}
\arguments{
\item{X}{An adjacency matrix}

\item{qdata.pep}{A data.frame of quantitative data}
}
\value{
A list of two items
}
\description{
Method to compute the detailed number of quantified peptides used for
aggregating each protein
}
\examples{
library(DaparToolshedData)
data(Exp1_R25_pept, package="DaparToolshedData")
obj.pep <- Exp1_R25_pept[seq_len(10)]
last.obj <- obj.pep[[length(obj.pep)]]
X <- BuildAdjacencyMatrix(last.obj)
ll.n <- GetDetailedNbPeptidesUsed(assay(last.obj), X)

}
\author{
Samuel Wieczorek
}
