% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregation.R
\name{DaparToolshed-aggregate}
\alias{DaparToolshed-aggregate}
\alias{aggregateFeatures4Prostar,QFeatures-method}
\alias{aggregateFeatures4Prostar,SummarizedExperiment-method}
\alias{aggQmetacell}
\alias{aggregateMethods}
\title{Aggregate an assay's quantitative features}
\usage{
\S4method{aggregateFeatures4Prostar}{QFeatures}(
  object,
  i,
  fcol,
  name = "newAssay",
  fun = MsCoreUtils::robustSummary,
  shared = TRUE,
  n = NULL,
  ...
)

\S4method{aggregateFeatures4Prostar}{SummarizedExperiment}(
  object,
  fcol,
  fun = MsCoreUtils::robustSummary,
  conds,
  shared = TRUE,
  n = NULL,
  ...
)

aggQmetacell(qMeta, X, level, conds)

aggregateMethods()
}
\arguments{
\item{object}{An instance of class \code{QFeatures} or \code{SummarizedExperiment}}

\item{i}{The index or name of the assay which features will be aggregated the create the new assay.}

\item{fcol}{A \code{character(1)} naming a rowdata variable (of assay \code{i} in case of a \code{QFeatures})
defining how to aggregate the features of the assay.
This variable is a (possibly sparse) matrix. See below for details.}

\item{name}{A \code{character(1)} naming the new assay. Default is \code{newAssay}.
Note that the function will fail if there's already an assay with \code{name}.}

\item{fun}{A function used for quantitative feature aggregation.
See details for examples.}

\item{shared}{A \code{boolean} indication if shared peptides should be considered. If \code{TRUE}, shared peptides}

\item{n}{A \code{numeric(1)} specifying the number of peptides to use for each protein. If \code{NULL}, all peptides are considered.}

\item{...}{Additional parameters passed the \code{fun}.}

\item{conds}{A \code{character()} vector which is the names of conditions
for each sample in the dataset.}

\item{qMeta}{An object of class 'SummarizedExperiment'}

\item{X}{xxxx}

\item{level}{A \code{character(1)} which is the type of dataset}
}
\value{
A \code{QFeatures} object with an additional assay or a \code{SummarizedExperiment} object (or subclass thereof).

xxxxx
}
\description{
This function aggregates the quantitative features of an assay,
applying a summarization function (\code{fun}) to sets of features.
The \code{fcol} variable name points to a rowData column that defines
how to group the features during aggregate. This variable has to
be an adjacency matrix. This function uses \code{\link[QFeatures:QFeatures-aggregate]{QFeatures::aggregateFeatures()}}
to aggregate quantitative data.

The list of agregation methods can be obtained with the function
\code{\link[=aggregateMethods]{aggregateMethods()}}. This function compiles both methods from the
packages \code{DaparToolshed} and \code{QFeatures}.

xxx
}
\details{
This function uses \code{\link[QFeatures:QFeatures-aggregate]{QFeatures::aggregateFeatures()}} to aggregate quantitative data.
}
\section{Iterative aggregation function}{

xxxxxx
xxxxx
}

\section{Quantitative metadata aggregation}{

The function to aggregate the quantitative metadata is \code{aggQmetadat()}.
}

\examples{

## ---------------------------------------
## An example QFeatures with PSM-level data
## ---------------------------------------
\dontrun{
data(ft, package='DaparToolshed')
ft

## Aggregate peptides into proteins
## using the adjacency matrix
feat1 <- aggregateFeatures4Prostar(object = ft,
i = 1,
name = 'aggregated',
fcol = 'adjacencyMatrix',
fun = 'colSumsMat')
feat1

assay(feat1[[1]])
assay(feat1[[2]])
aggcounts(feat1[[2]])
assay(feat1[[3]])
aggcounts(feat1[[3]])
rowData(feat1[[2]])
}
data(ft, package='DaparToolshed')
qMeta <- qMetacell(ft, 1)
X <- adjacencyMatrix(ft[[1]])
level <- typeDataset(ft[[1]])
conds <- colData(ft)$Condition
aggQmeta <- aggQmetacell(qMeta, X, level, conds)

}
\seealso{
The \emph{QFeatures} vignette provides an extended example and
the \emph{Aggregation} vignette, for a complete quantitative
proteomics data processing pipeline.
}
