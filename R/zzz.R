.onLoad <- function(libname, pkgname) {
    # shiny::addResourcePath(
    #     prefix = "images",
    #     directoryPath = system.file("images",package = "DaparToolshed")
    # )
  
  #library(QFeatures)
}



#' @title List of available pipelines in the package
#'
#' @description
#' Get the list of pipelines available in the package
#'
#' @export
#'
#' @return NA
#'
Pipelines <- function() {
    list(
        Protein = c("protein"),
        Peptide = c("peptide"),
        P2p = c("protein"),
        Peptidomic = c("peptide")
    )
}



#' @title Adds an instance of `SummarizedExperiment` to a `QFeatures` object
#'
#' @description
#' Adds one or more items to the dataset. This function is specific of the
#' type of dataset.
#'
#' @param dataset An instance of `SummarizedExperiment` class
#'
#' @param name A `character()` for the new assay
#'
#' @importFrom QFeatures addAssay
#'
#' @return The dataset minus some items
#'
#' @export
#'
Add_Item_to_Dataset <- function(dataset, name) {
    QFeatures::addAssay(dataset,
        dataset[[length(dataset)]],
        name = name
    )
}

#' @title Removes assay from a `QFeatures` object
#'
#' @description
#' Removes one or more items from the dataset. This function is specific of the
#' type of dataset.
#'
#' @param dataset xxx
#'
#' @param range xxx
#'
#' @return The dataset minus some items
#'
#' @export
#'
Keep_Items_from_Dataset <- function(dataset, range) {
    dataset[, , range]
}
