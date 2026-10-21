#' Maps of the variable and of its local feature contributions
#'
#' Draws three maps: the ranks of the variable, the local level feature
#' \code{c.lin|n.lin} and the local scale feature \code{c.quad|n.quad} of
#' \eqn{T_V(4,2)} (see \code{\link{local_features}}), with a blue-red scale.
#' Polygons are drawn when the graph was built by \code{\link{polygon_graph}};
#' otherwise the sites are drawn as points (squares for lattices) at their
#' coordinates. The local contributions are descriptive, not tests.
#'
#' @param x The variable at the sites.
#' @param graph Neighbourhood graph with coordinates.
#' @param seed Optional seed (ties are broken at random).
#' @param main Optional overall title.
#' @return The local feature contributions (a data frame), invisibly.
#' @examples
#' g <- lattice_graph(20, 25)
#' map_local(mercer_hall$grain, g)
#' @export
map_local <- function(x, graph, seed = NULL, main = NULL) {
  check_graph(graph)
  if (is.null(graph$coords)) stop("this graph has no coordinates (edge list): no map can be drawn", call. = FALSE)
  lf <- local_features(x, graph, d = 4, seed = seed)
  rk <- (rank(x, ties.method = "average") - 0.5) / length(x) - 0.5
  clip <- function(v) pmin(pmax(v, -3), 3)
  panels <- list(list("ranks of the variable", rk, c(-0.5, 0.5)),
                 list("local level feature c.lin|n.lin", clip(lf[["c.lin|n.lin"]]), c(-3, 3)),
                 list("local scale feature c.quad|n.quad", clip(lf[["c.quad|n.quad"]]), c(-3, 3)))
  pal <- grDevices::hcl.colors(101, "Blue-Red 3")
  colour <- function(v, lim) pal[1 + round(100 * (v - lim[1]) / diff(lim))]
  op <- graphics::par(mfrow = c(1, 3), mar = c(1, 1, 2.5, 1), oma = c(2.5, 0, if (is.null(main)) 0 else 2, 0))
  on.exit(graphics::par(op))
  xy <- graph$coords
  geo <- graph$geometry
  rings <- NULL
  if (!is.null(geo) && requireNamespace("sf", quietly = TRUE)) {
    cc <- sf::st_coordinates(geo)
    lc <- grep("^L[0-9]+$", colnames(cc))
    rings <- split(seq_len(nrow(cc)), do.call(paste, c(as.data.frame(cc[, lc, drop = FALSE]), sep = "_")))
    owner <- vapply(rings, function(i) cc[i[1], lc[length(lc)]], numeric(1))
  }
  for (p in panels) {
    col <- colour(p[[2]], p[[3]])
    if (!is.null(rings)) {
      graphics::plot(range(cc[, 1]), range(cc[, 2]), type = "n", asp = 1, axes = FALSE, xlab = "", ylab = "")
      for (k in seq_along(rings))
        graphics::polygon(cc[rings[[k]], 1], cc[rings[[k]], 2], col = col[owner[k]], border = "white", lwd = 0.2)
    } else {
      if (grepl("lattice", graph$kind)) {
        graphics::plot(range(xy[, 1]) + c(-0.5, 0.5), range(xy[, 2]) + c(-0.5, 0.5), type = "n", asp = 1,
                       axes = FALSE, xlab = "", ylab = "")
        graphics::rect(xy[, 1] - 0.5, xy[, 2] - 0.5, xy[, 1] + 0.5, xy[, 2] + 0.5, col = col, border = "white",
                       lwd = 0.3)
      } else {
        graphics::plot(xy[, 1], xy[, 2], pch = 21, bg = col, col = "grey40", lwd = 0.3, asp = 1, axes = FALSE,
                       xlab = "", ylab = "", cex = max(0.5, min(2, 25 / sqrt(graph$n))))
      }
    }
    graphics::title(p[[1]], cex.main = 1)
    graphics::box(col = "grey80")
  }
  graphics::mtext("blue: low / negative    red: high / positive    (local contributions: mean 0, variance 1 under independence)",
                  side = 1, outer = TRUE, line = 1, cex = 0.75)
  if (!is.null(main)) graphics::mtext(main, side = 3, outer = TRUE, line = 0.5, font = 2)
  invisible(lf)
}
