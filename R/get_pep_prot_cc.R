#' @title Build the list of connex composant of the adjacency matrix
#'
#' @param X An adjacency matrix
#'
#' @return A list of CC
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
#' ll <- getPepProtCC(X)
#'
#' @export
#' @importFrom Matrix crossprod Matrix
#' @importFrom igraph components graph_from_adjacency_matrix
#' @importFrom graph graphAM connComp 
#'
getPepProtCC <- function(X) {
  
  if (is.null(X)) {
    warning("The adjacency matrix is empty")
    return()
  }
  
  # Protein-protein adjacency (shared peptides), sparse, no self-loops
  A <- Matrix::crossprod(X, boolArith = TRUE)
  Matrix::diag(A) <- FALSE
  
  # Connected components via igraph (works for single-protein CCs too)
  g <- igraph::graph_from_adjacency_matrix(A, mode = "undirected")
  membership <- igraph::components(g)$membership   # named by protein
  
  prot.cc <- split(names(membership), membership)
  
  # Sub-matrix of X for each CC: keep the CC's proteins, then only the
  # peptides that have at least one non-zero entry in those columns
  res <- lapply(prot.cc, function(prots) {
    sub <- X[, prots, drop = FALSE]
    sub[Matrix::rowSums(sub) != 0, , drop = FALSE]
  })
  
  # Order CCs by their (character-sorted) first protein name, as in your target
  first.prot <- vapply(prot.cc, function(p) sort(p)[1], character(1))
  res <- res[order(first.prot)]
  
  return(unname(res))
}





#' @title Jitter plot of CC
#'
#' @param list.of.cc List of cc such as returned by the function getPepProtCC
#'
#' @return A plot
#'
#' @author Thomas Burger
#'
#' @examples
#' data(subR25pept)
#' X <- BuildAdjacencyMatrix(subR25pept[[1]])
#' ll <- getPepProtCC(X)
#' plotJitter(ll)
#'
#' @export
#'
plotJitter <- function(list.of.cc = NULL) {
  if (is.null(list.of.cc)) {
    return()
  }

  # One column per CC: number of proteins (columns) and peptides (rows)
  cc.summary <- vapply(
    list.of.cc,
    function(x) c(ncol(x), nrow(x)),
    integer(2)
  )

  rownames(cc.summary) <- c("Nb_proteins", "Nb_peptides")
  colnames(cc.summary) <- seq_along(list.of.cc)

  cc.summary <- as.data.frame(t(cc.summary))

  plot(
    jitter(cc.summary[, "Nb_peptides"]),
    jitter(cc.summary[, "Nb_proteins"]),
    type = "p",
    xlab = "#peptides in CC",
    ylab = "#proteins in CC"
  )
}



#' @title Display a CC
#'
#' @param The.CC A cc (a list)
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
#' ll <- getPepProtCC(X)
#' g <- buildGraph(ll[[1]])
#'
#' @export
#' @importFrom Matrix Matrix
#'
buildGraph <- function(The.CC) {
  # The.CC is the sparse sub-matrix of one connected component
  subX <- as.matrix(The.CC)
  
  nb.pep  <- nrow(subX)
  nb.prot <- ncol(subX)
  nb.total <- nb.prot + nb.pep
  
  pep.deg <- rowSums(subX != 0)   # number of proteins per peptide
  
  edge.list <- as.data.frame(which(subX != 0, arr.ind = TRUE))
  
  # Peptide groups: shared if linked to >1 protein, specific if exactly 1
  def.grp <- c(
    ifelse(pep.deg == 1, "spec.peptide", "shared.peptide"),
    rep("protein", nb.prot)
  )
  
  nodes <- data.frame(
    id = seq_len(nb.total),
    group = def.grp,
    label = c(rownames(subX), colnames(subX)),
    title = paste0("<p>", seq_len(nb.total), "<br>Tooltip !</p>"),
    size = c(rep(10, nb.pep), rep(20, nb.prot)),
    stringsAsFactors = FALSE
  )
  
  edges <- data.frame(
    from = edge.list$row,
    to = edge.list$col + nb.pep,
    stringsAsFactors = FALSE
  )
  
  list(nodes = nodes, edges = edges)
}


#' @title Display a CC
#'
#' @param g A cc (a list)
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @examples
#' data(subR25pept)
#' X <- QFeatures::adjacencyMatrix(subR25pept[[1]])
#' ll <- getPepProtCC(X)
#' g <- buildGraph(ll[[1]])
#' displayCCvisNet(g)
#'
#' @export
#'
#'
displayCCvisNet <- function(g) {
  
  pkgsRequire('visNetwork')
  
  col.prot <- "#ECB57C"
  col.spec <- "#5CA3F7"
  col.shared <- "#0EA513"
  
  
  visNetwork::visNetwork(g$nodes, g$edges, width = "100%", 
                         height = "100%") |>
    visNetwork::visNodes(shape = "dot") |> # square for all nodes
    visNetwork::visGroups(groupname = "spec.peptide", 
                          color = col.spec) |> # darkblue for group "A"
    visNetwork::visGroups(groupname = "shared.peptide", 
                          color = col.shared) |> # darkblue for group "A"
    visNetwork::visGroups(groupname = "protein", 
                          color = col.prot, shape = "dot") |>
    visNetwork::visOptions(highlightNearest = FALSE) |>
    visNetwork::visEdges(color = "#A9A9A9", width = 2)
}



#' @title Display a a jitter plot for CC
#'
#' @param df A `data.frame()`
#'
#' @return A plot
#'
#' @author Thomas Burger, Samuel Wieczorek
#'
#' @export
#' 
#' @examples 
#' data(subR25pept)
#' X <- BuildAdjacencyMatrix(subR25pept[[1]])
#' ll <- getPepProtCC(X)[seq_len(4)]
#' n.prot <- vapply(ll, ncol, integer(1))
#' n.pept <- vapply(ll, nrow, integer(1))
#' df <- tibble::tibble(
#' x = jitter(n.pept),
#' y = jitter(n.prot),
#' index = seq_len(length(ll))
#' )
#' plotJitter_rCharts(df)
#'
plotJitter_rCharts <- function(df) {
  
  i_tooltip <- which(startsWith(colnames(df), "tooltip"))
  txt_tooltip <- NULL
  
  if (length(i_tooltip) == 0){
    warning("There is no tooltip in the object.")
  }
  
  txt_tooltip <- apply(df[, i_tooltip, drop = FALSE], 1, function(row) {
    paste(
      paste0("<b>", gsub("tooltip_", "", colnames(df)[i_tooltip]), "</b>: ", row),
      collapse = "<br>"
    )
  })
  
  p <- plot_ly(
    data = df,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "markers",
    text = txt_tooltip,
    hoverinfo = "text",
    showlegend = FALSE
  ) |>
    plotly::layout(
      xaxis = list(title = "Nb of peptides ic CC"),
      yaxis = list(title = "Nb of proteins ic CC"),
      hovermode = "closest",
      margin = list(b = 60) 
    )
  
  return(p)
}
