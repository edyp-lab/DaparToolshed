
#' @title Display a heatmap for data with missing values
#' @param x A numeric matrix
#' @param col See `graphics::image()`
#' @param srtCol See `graphics::text()`
#' @param labCol See `graphics::text()`
#' @param labRow See `graphics::axis()`
#' @param key.title See `graphics::title()`
#' @param key See `graphics::par()`
#' @param main See `graphics::title()`
#' @param ylab See `graphics::image()`
#' @rdname heatmap
#' @export
#' @return A heatmap
#' @import graphics
#' @import grDevices
#'
heatmapForMissingValues <- function(
    x,
  col = NULL,
  srtCol = NULL,
  labCol = NULL,
  labRow = NULL,
  key = TRUE,
  key.title = NULL,
  main = NULL,
  ylab = NULL) {
  if (is.null(col)) {
    col <- grDevices::heat.colors(100)
  }
  
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low) / (high - low)
    x
  }
  
  offsetCol <- 0.5
  offsetRow <- 0.5
  srtRow <- NULL
  colRow <- NULL
  colCol <- NULL
  xlab <- NULL
  key.par <- list()
  margins <- c(5, 5)
  sepcolor <- "white"
  na.color <- "white"
  keysize <- 1.5
  breaks <- NULL
  na.rm <- TRUE
  
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) {
    stop("`x' must be a numeric matrix")
  }
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) {
    stop("`x' must have at least 2 rows and 2 columns")
  }
  x <- x[nr:1, ]
  cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  cexCol <- 0.2 + 1 / log10(nc)
  cexRow <- 0.2 + 1 / log10(nr)
  iy <- seq_len(nr)
  breaks <- length(col) + 1
  breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
    length = breaks
  )
  
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  lhei <- c(keysize, 4)
  lwid <- c(keysize, 4)
  lmat <- rbind(4:3, 2:1)
  lmat[is.na(lmat)] <- 0
  
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  graphics::par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  
  
  graphics::image(seq_len(nc), seq_len(nr), x,
    xlim = 0.5 + c(0, nc), ylim = 0.5 +
      c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col,
    breaks = breaks
  )
  
  
  if (!is.null(labCol)) {
    graphics::axis(1,
      seq_len(nc),
      label = labCol,
      las = 2,
      line = -0.5 + offsetCol,
      tick = 0,
      cex.axis = cexCol,
      hadj = NA,
      padj = 0
    )
  } else {
    adjCol <- c(1, NA)
    xpd.orig <- graphics::par("xpd")
    graphics::par(xpd = NA)
    xpos <- graphics::axis(1,
      seq_len(nc),
      label = rep("", nc),
      las = 2,
      tick = 0
    )
    graphics::text(
      x = xpos,
      y = graphics::par("usr")[3] - (1 + offsetCol) *
        graphics::strheight("M"),
      label = labCol,
      adj = adjCol,
      cex = cexCol,
      srt = srtCol,
      col = colCol
    )
    graphics::par(xpd = xpd.orig)
  }
  
  
  if (!is.null(labRow)) {
    graphics::axis(4,
      iy,
      label = labRow,
      las = 5,
      line = -0.5 + offsetRow,
      tick = 0,
      cex.axis = cexRow,
      hadj = 0,
      padj = NA
    )
  } else {
    xpd.orig <- graphics::par("xpd")
    graphics::par(xpd = NA)
    ypos <- graphics::axis(4,
      iy,
      label = rep("", nr),
      las = 2,
      line = -0.5,
      tick = 0
    )
    
    .strw <- graphics::strwidth("M")
    graphics::text(
      x = graphics::par("usr")[2] + (1 + offsetRow) * .strw,
      y = ypos,
      label = labRow,
      adj = c(0, NA),
      cex = cexRow,
      srt = srtRow,
      col = colRow
    )
    graphics::par(xpd = xpd.orig)
  }
  
  graphics::par(mar = c(margins[1], 0, 0, 0))
  graphics::plot.new()
  graphics::par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  
  graphics::plot.new()
  if (!is.null(main)) {
    graphics::title(main, cex.main = 1.5 * op[["cex.main"]])
  }
  
  
  if (key) {
    mar <- c(5, 4, 2, 1)
    graphics::par(mar = mar, cex = 0.75, mgp = c(2, 1, 0))
    if (length(key.par) > 0) {
      do.call(par, key.par)
    }
    
    tmpbreaks <- breaks
    min.raw <- min.breaks
    max.raw <- max.breaks
    
    z <- seq(min.raw, max.raw, by = min(diff(breaks) / 100))
    graphics::image(
      z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
      xaxt = "n", yaxt = "n"
    )
    graphics::par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    xargs <- list(at = xv, label = lv)
    
    xargs$side <- 1
    do.call(graphics::axis, xargs)
    key.xlab <- "Intensity value"
    
    graphics::mtext(
      side = 1,
      key.xlab,
      line = graphics::par("mgp")[1],
      padj = 0.5,
      cex = graphics::par("cex") * graphics::par("cex.lab")
    )
    
    if (is.null(key.title)) {
      graphics::title("Color Key")
    }
  }
}
