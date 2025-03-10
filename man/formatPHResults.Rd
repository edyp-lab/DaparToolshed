% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anova_analysis.R
\name{formatPHResults}
\alias{formatPHResults}
\title{Extract logFC and raw pvalues from multiple post-hoc models summaries}
\usage{
formatPHResults(post_hoc_models_summaries)
}
\arguments{
\item{post_hoc_models_summaries}{a list of summaries of post-hoc models.}
}
\value{
a list of 2 dataframes containing the logFC values and pvalues for
each comparison.
}
\description{
Extract logFC and raw pvalues from multiple post-hoc models summaries
}
\examples{
\dontrun{examples/ex_formatPHResults.R}

}
\author{
Helene Borges
}
