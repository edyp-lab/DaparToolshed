% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DaparToolshed_aggregate.R
\name{RunAggregation}
\alias{RunAggregation}
\title{xxx}
\usage{
RunAggregation(
  qf = NULL,
  i = 0,
  includeSharedPeptides = "Yes_As_Specific",
  operator = "Sum",
  considerPeptides = "allPeptides",
  n = NULL,
  addRowData = FALSE
)
}
\arguments{
\item{qf}{xxx}

\item{i}{xxx}

\item{includeSharedPeptides}{A boolean}

\item{operator}{xxx}

\item{considerPeptides}{Available values are 'allPeptides' (default) and 'topN'}

\item{n}{In case of "top n peptides', specify then number of peptides}

\item{addRowData}{xxx}
}
\description{
xxx
A short description...
}
\examples{
\dontrun{
data(Exp1_R25_pept, package="DaparToolshedData")
ft <- Exp1_R25_pept[1:100]
obj.agg <- RunAggregation(ft, length(ft), "Yes_As_Specific", 'Sum', 'allPeptides', addRowData = TRUE)
obj.agg <- RunAggregation(ft, length(ft), "Yes_As_Specific", 'Mean', 'allPeptides', addRowData = TRUE)
obj.agg <- RunAggregation(ft, length(ft), "Yes_As_Specific", 'Sum', "topN", n = 4, addRowData = TRUE)
obj.agg <- RunAggregation(ft, length(ft), "Yes_As_Specific", 'Mean', "topN", n = 4, addRowData = TRUE)

# Case E
obj.agg <- RunAggregation(ft, length(ft), "No", 'Sum', 'allPeptides', addRowData = FALSE)

obj.agg <- RunAggregation(ft, length(ft), "No", 'Sum', "topN", n = 4)
obj.agg <- RunAggregation(ft, length(ft), "Yes_Redistribution", 'Sum', 'allPeptides', addRowData = TRUE)
obj.agg <- RunAggregation(ft, length(ft), "Yes_Redistribution", 'Sum', "topN", n = 4, addRowData = TRUE)
}

}
