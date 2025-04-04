% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/agregation.R
\name{aggregateIterParallel}
\alias{aggregateIterParallel}
\title{xxxx}
\usage{
aggregateIterParallel(
  obj.pep,
  i,
  X,
  init.method = "Sum",
  method = "Mean",
  n = NULL
)
}
\arguments{
\item{obj.pep}{xxxxx}

\item{i}{xxx}

\item{X}{xxxx}

\item{init.method}{xxxxx}

\item{method}{xxxxx}

\item{n}{xxxx}
}
\value{
xxxxx
}
\description{
xxxx
}
\examples{
\dontrun{
data(Exp1_R25_pept, package="DaparToolshedData")
obj.pep <- Exp1_R25_pept[seq_len(10)]
last.se <- obj.pep[[length(obj.pep)]]
X <- BuildAdjacencyMatrix(last.se)
conds <- colData(obj.pep)$Condition
obj.agg <- aggregateIterParallel(obj.pep = obj.pep, i = length(obj.pep), X = X, n = 3)
}

}
\author{
Samuel Wieczorek
}
