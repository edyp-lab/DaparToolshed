.onLoad <- function(libname, pkgname) {
    # shiny::addResourcePath(
    #     prefix = "images",
    #     directoryPath = system.file("images",package = "DaparToolshed")
    # )
  
  #library(QFeatures)
}

.onUnload <- function(libname, pkgname) {
    shiny::removeResourcePath("images")
}



#' @title xxx
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



#' @title xxx
#'
#' @description
#' Adds one or more items to the dataset. This function is specific of the
#' type of dataset.
#'
#' @param dataset xxx
#'
#' @param name xxx
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

#' @title xxx
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
