% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregation.R
\name{GetDetailedNbPeptides}
\alias{GetDetailedNbPeptides}
\title{Computes the detailed number of peptides for each protein}
\usage{
GetDetailedNbPeptides(X)
}
\arguments{
\item{X}{An adjacency matrix}
}
\value{
A data.frame
}
\description{
Method to compute the detailed number of quantified peptides for each
protein
}
\examples{
\dontrun{
data(Exp1_R25_pept, package="DaparToolshedData")
obj.pep <- Exp1_R25_pept[seq_len(10)]
last.obj <- obj.pep[[length(obj.pep)]]
X <- BuildAdjacencyMatrix(last.obj)
n <- GetDetailedNbPeptides(X)
}

}
\author{
Samuel Wieczorek
}
