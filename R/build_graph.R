#' Build the neighbourhood graph from a data set
#'
#' A single entry point to the graph builders, used by the graphical
#' interface: the graph is built from the columns of a data frame (or from the
#' polygons of an \pkg{sf} data frame).
#'
#' @param data A data frame with one row per site, or an \pkg{sf} data frame
#'   of polygons.
#' @param type \code{"knn"}, \code{"delaunay"} or \code{"band"} (from the
#'   coordinate columns \code{x}, \code{y}); \code{"lattice"} (from the
#'   row and column index columns \code{row}, \code{col}); \code{"grid"}
#'   (\code{nrow} x \code{ncol} grid, rows stored row by row);
#'   \code{"polygons"} (contiguity of the polygons of an \pkg{sf} data
#'   frame); \code{"edges"} (an edge list \code{edges}).
#' @param x,y Names of the coordinate columns.
#' @param k Number of nearest neighbours.
#' @param radius Distance of the distance band.
#' @param row,col Names of the row and column index columns.
#' @param nrow,ncol Dimensions of the grid.
#' @param contiguity \code{"rook"} or \code{"queen"} (default: rook for
#'   lattices and grids, queen for polygons).
#' @param edges A data frame or matrix whose first two columns are site
#'   numbers (1-based), or the path of a CSV file with them.
#' @return An object of class \code{"ivspatial_graph"}.
#' @examples
#' build_graph(baltimore, "knn", x = "X", y = "Y", k = 4)
#' build_graph(mercer_hall, "lattice", row = "row", col = "col")
#' @export
build_graph <- function(data, type = c("knn", "delaunay", "band", "lattice", "grid", "polygons", "edges"),
                        x = "x", y = "y", k = 4, radius = 0, row = "row", col = "col", nrow = 0, ncol = 0,
                        contiguity = NULL, edges = NULL) {
  type <- match.arg(type)
  if (is.null(contiguity)) contiguity <- if (type == "polygons") "queen" else "rook"
  df <- if (inherits(data, "sf")) sf::st_drop_geometry(data) else as.data.frame(data)
  need <- function(v, what) {
    if (!(v %in% names(df))) stop(sprintf("%s column '%s' not found; available: %s", what, v,
                                          paste(names(df), collapse = ", ")), call. = FALSE)
    as.numeric(df[[v]])
  }
  if (type %in% c("knn", "delaunay", "band")) {
    xy <- cbind(need(x, "coordinate"), need(y, "coordinate"))
    if (!all(is.finite(xy))) stop("missing coordinates", call. = FALSE)
    return(switch(type, knn = knn_graph(xy, k), delaunay = delaunay_graph(xy),
                  band = {
                    if (radius <= 0) stop("give a positive distance for the distance band", call. = FALSE)
                    distance_band_graph(xy, radius)
                  }))
  }
  if (type == "lattice") return(lattice_from_indices(need(row, "row index"), need(col, "column index"), contiguity))
  if (type == "grid") {
    if (nrow * ncol != nrow(df)) stop(sprintf("nrow x ncol = %d but the data have %d rows", nrow * ncol, nrow(df)),
                                      call. = FALSE)
    return(lattice_graph(nrow, ncol, contiguity))
  }
  if (type == "polygons") {
    if (!inherits(data, "sf")) stop("polygon contiguity needs polygons (a shapefile, GeoJSON or GeoPackage file)",
                                    call. = FALSE)
    return(polygon_graph(data, contiguity))
  }
  if (is.null(edges)) stop("give the edge list", call. = FALSE)
  e <- if (is.character(edges)) read_sites(edges) else as.data.frame(edges)
  if (ncol(e) < 2) stop("the edge list needs two columns of site numbers", call. = FALSE)
  graph_from_edges(e[[1]], e[[2]], n = nrow(df))
}
