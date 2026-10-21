# Neighbourhood graphs: points, lattices, polygons, edge lists and spatial lags.
#
# An object of class "ivspatial_graph" is a list with
#   neighbors  list of integer vectors (1-based), neighbors[[v]] = N(v) without v,
#   n, degree, kind (a short description), coords (n x 2 matrix or NULL).

new_graph <- function(neighbors, kind = "graph", coords = NULL) {
  n <- length(neighbors)
  nb <- lapply(seq_len(n), function(v) {
    a <- sort(unique(as.integer(neighbors[[v]])))
    a <- a[!is.na(a) & a != 0L]
    if (length(a) && (min(a) < 1L || max(a) > n)) stop("neighbour index out of range at site ", v, call. = FALSE)
    a[a != v]
  })
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    storage.mode(coords) <- "double"
    if (nrow(coords) != n) stop("coords must have one row per site", call. = FALSE)
  }
  structure(list(neighbors = nb, n = n, degree = lengths(nb), kind = kind, coords = coords),
            class = "ivspatial_graph")
}

#' Neighbourhood graphs
#'
#' Functions that build the neighbourhood graph on which the tests are
#' computed. Every site \eqn{v} has a set of neighbours \eqn{N(v)}; the
#' neighbourhood of a centre is \eqn{v} together with \eqn{N(v)}. The graph
#' may be directed (k nearest neighbours) or undirected (contiguity,
#' Delaunay, lattices). The order of the sites is that of the data.
#'
#' \describe{
#'   \item{\code{knn_graph}}{each site and its \code{k} nearest neighbours
#'     (Euclidean distance; ties broken by index);}
#'   \item{\code{delaunay_graph}}{the Delaunay triangulation of the points,
#'     that is the contiguity graph of their Voronoi cells (needs the package
#'     \pkg{deldir});}
#'   \item{\code{distance_band_graph}}{all sites within distance
#'     \code{radius};}
#'   \item{\code{lattice_graph}}{an \code{nrow} x \code{ncol} grid stored row
#'     by row (site index \code{(i - 1) * ncol + j}), rook or queen
#'     contiguity;}
#'   \item{\code{lattice_from_indices}}{a grid given by integer row and
#'     column indices of the sites (cells may be missing);}
#'   \item{\code{graph_from_edges}}{an edge list;}
#'   \item{\code{graph_from_neighbors}}{a list of neighbour indices, for
#'     example an \code{nb} object of \pkg{spdep};}
#'   \item{\code{polygon_graph}}{queen or rook contiguity of polygons (an
#'     \pkg{sf} object or a polygon file; see \code{\link{polygon_graph}}).}
#' }
#'
#' @param coords A matrix or data frame with two columns of coordinates, one
#'   row per site.
#' @param k Number of nearest neighbours.
#' @param radius Distance of the distance band.
#' @param nrow,ncol Numbers of rows and columns of the grid.
#' @param rows,cols Integer row and column indices of the sites.
#' @param contiguity \code{"rook"} (common edge) or \code{"queen"} (common
#'   edge or corner).
#' @param i,j Vectors of site numbers (1-based) of the edges.
#' @param n Number of sites (default: the largest site number).
#' @param symmetric If \code{TRUE} (default), every edge is used in both
#'   directions.
#' @param neighbors A list whose element \code{v} gives the neighbours of site
#'   \code{v} (1-based; 0 means none, as in \pkg{spdep}).
#' @param kind A short description of the graph, used in reports.
#' @return An object of class \code{"ivspatial_graph"}.
#' @examples
#' g <- lattice_graph(20, 25)
#' g
#' w <- mercer_hall
#' lattice_from_indices(w$row, w$col, "queen")
#' b <- baltimore
#' knn_graph(cbind(b$X, b$Y), k = 4)
#' graph_from_edges(c(1, 2, 3), c(2, 3, 1))
#' @name graphs
NULL

check_coords <- function(coords) {
  xy <- as.matrix(coords)
  if (ncol(xy) < 2) stop("coords must have two columns", call. = FALSE)
  xy <- xy[, 1:2, drop = FALSE]
  storage.mode(xy) <- "double"
  if (!all(is.finite(xy))) stop("coordinates must be finite", call. = FALSE)
  unname(xy)
}

#' @rdname graphs
#' @export
knn_graph <- function(coords, k = 4) {
  xy <- check_coords(coords)
  n <- nrow(xy)
  if (k < 1 || k >= n || k != round(k)) stop("k must be a whole number between 1 and n - 1", call. = FALSE)
  if (anyDuplicated(xy)) warning("duplicate coordinates: their order among nearest neighbours is arbitrary", call. = FALSE)
  nb <- vector("list", n)
  step <- max(1L, floor(2e6 / n))
  for (s in seq(1L, n, by = step)) {
    idx <- s:min(n, s + step - 1L)
    D <- outer(xy[idx, 1], xy[, 1], "-")^2 + outer(xy[idx, 2], xy[, 2], "-")^2
    for (r in seq_along(idx)) {
      v <- idx[r]
      D[r, v] <- Inf
      nb[[v]] <- order(D[r, ], seq_len(n))[seq_len(k)]
    }
  }
  new_graph(nb, sprintf("%d-nearest-neighbour graph", as.integer(k)), xy)
}

#' @rdname graphs
#' @export
delaunay_graph <- function(coords) {
  xy <- check_coords(coords)
  n <- nrow(xy)
  if (n < 4) stop("at least 4 points are needed", call. = FALSE)
  if (anyDuplicated(xy))
    stop("duplicate coordinates: Delaunay triangulation needs distinct points ",
         "(use a k-nearest-neighbour graph or jitter the points)", call. = FALSE)
  if (!requireNamespace("deldir", quietly = TRUE))
    stop("the Delaunay graph needs the package deldir: install.packages(\"deldir\")", call. = FALSE)
  dd <- suppressWarnings(deldir::deldir(xy[, 1], xy[, 2]))
  e <- dd$delsgs
  graph_edges(e$ind1, e$ind2, n, TRUE, "Delaunay graph", xy)
}

#' @rdname graphs
#' @export
distance_band_graph <- function(coords, radius) {
  xy <- check_coords(coords)
  if (!is.finite(radius) || radius <= 0) stop("radius must be positive", call. = FALSE)
  n <- nrow(xy)
  nb <- vector("list", n)
  step <- max(1L, floor(2e6 / n))
  for (s in seq(1L, n, by = step)) {
    idx <- s:min(n, s + step - 1L)
    D <- outer(xy[idx, 1], xy[, 1], "-")^2 + outer(xy[idx, 2], xy[, 2], "-")^2
    for (r in seq_along(idx)) nb[[idx[r]]] <- setdiff(which(D[r, ] <= radius^2), idx[r])
  }
  new_graph(nb, sprintf("distance band %g", radius), xy)
}

lattice_offsets <- function(contiguity) {
  contiguity <- match.arg(contiguity, c("rook", "queen"))
  off <- rbind(c(0, 1), c(1, 0), c(0, -1), c(-1, 0))
  if (contiguity == "queen") off <- rbind(off, c(1, 1), c(1, -1), c(-1, 1), c(-1, -1))
  off
}

#' @rdname graphs
#' @export
lattice_graph <- function(nrow, ncol, contiguity = c("rook", "queen")) {
  contiguity <- match.arg(contiguity)
  if (nrow < 1 || ncol < 1) stop("nrow and ncol must be positive", call. = FALSE)
  off <- lattice_offsets(contiguity)
  ii <- rep(seq_len(nrow), each = ncol)
  jj <- rep(seq_len(ncol), times = nrow)
  nb <- lapply(seq_along(ii), function(v) {
    r <- ii[v] + off[, 1]
    c <- jj[v] + off[, 2]
    ok <- r >= 1 & r <= nrow & c >= 1 & c <= ncol
    (r[ok] - 1L) * ncol + c[ok]
  })
  new_graph(nb, sprintf("%d x %d lattice, %s contiguity", as.integer(nrow), as.integer(ncol), contiguity),
            cbind(jj - 1, nrow - ii))
}

#' @rdname graphs
#' @export
lattice_from_indices <- function(rows, cols, contiguity = c("rook", "queen")) {
  contiguity <- match.arg(contiguity)
  rows <- as.integer(round(rows))
  cols <- as.integer(round(cols))
  if (length(rows) != length(cols) || anyNA(rows) || anyNA(cols))
    stop("rows and cols must be complete integer vectors of the same length", call. = FALSE)
  key <- paste(rows, cols)
  if (anyDuplicated(key)) stop("duplicate (row, column) cells", call. = FALSE)
  off <- lattice_offsets(contiguity)
  nb <- lapply(seq_along(rows), function(v) {
    m <- match(paste(rows[v] + off[, 1], cols[v] + off[, 2]), key)
    m[!is.na(m)]
  })
  new_graph(nb, sprintf("lattice (%d x %d cells), %s contiguity", diff(range(rows)) + 1L,
                        diff(range(cols)) + 1L, contiguity), cbind(cols, -rows))
}

graph_edges <- function(i, j, n, symmetric, kind, coords = NULL) {
  i <- as.integer(i)
  j <- as.integer(j)
  if (symmetric) {
    ii <- c(i, j)
    jj <- c(j, i)
  } else {
    ii <- i
    jj <- j
  }
  ok <- ii != jj
  nb <- split(jj[ok], factor(ii[ok], levels = seq_len(n)))
  new_graph(unname(nb), kind, coords)
}

#' @rdname graphs
#' @export
graph_from_edges <- function(i, j, n = NULL, symmetric = TRUE) {
  i <- as.integer(round(i))
  j <- as.integer(round(j))
  if (length(i) != length(j) || anyNA(i) || anyNA(j)) stop("i and j must be complete vectors of the same length",
                                                          call. = FALSE)
  if (is.null(n)) n <- max(i, j)
  if (min(i, j) < 1 || max(i, j) > n) stop("site numbers must be between 1 and n", call. = FALSE)
  graph_edges(i, j, n, symmetric, "graph from edge list")
}

#' @rdname graphs
#' @export
graph_from_neighbors <- function(neighbors, kind = "graph", coords = NULL) {
  if (!is.list(neighbors)) stop("neighbors must be a list", call. = FALSE)
  if (inherits(neighbors, "nb") && kind == "graph") kind <- "graph from spdep nb object"
  new_graph(unclass(neighbors), kind, coords)
}

#' Contiguity graph of polygons
#'
#' Queen (common vertex) or rook (common edge) contiguity of polygons, such as
#' counties or census tracts. Vertices are matched after rounding to
#' \code{tol}; this gives the same graph as \code{spdep::poly2nb} for polygons
#' whose common boundaries share their vertices, which is the case for
#' boundary files built from a common topology. Needs the package \pkg{sf}.
#'
#' @param polygons An \pkg{sf} object, an \code{sfc} geometry column, or the
#'   path of a polygon file readable by \code{sf::st_read} (shapefile,
#'   GeoJSON, GeoPackage, ...).
#' @param contiguity \code{"queen"} (default) or \code{"rook"}.
#' @param tol Tolerance for matching vertices.
#' @return An object of class \code{"ivspatial_graph"} with the centroids of
#'   the polygons as coordinates; the polygons are kept in the element
#'   \code{geometry} for maps.
#' @examples
#' if (requireNamespace("sf", quietly = TRUE)) {
#'   nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'   g <- polygon_graph(nc)
#'   g
#' }
#' @export
polygon_graph <- function(polygons, contiguity = c("queen", "rook"), tol = 1e-9) {
  contiguity <- match.arg(contiguity)
  if (!requireNamespace("sf", quietly = TRUE))
    stop("polygon contiguity needs the package sf: install.packages(\"sf\")", call. = FALSE)
  if (is.character(polygons)) polygons <- sf::st_read(polygons, quiet = TRUE)
  geom <- if (inherits(polygons, "sf")) sf::st_geometry(polygons) else polygons
  if (!inherits(geom, "sfc")) stop("polygons must be an sf object, an sfc column or a file name", call. = FALSE)
  n <- length(geom)
  xy <- sf::st_coordinates(geom)
  lc <- grep("^L[0-9]+$", colnames(xy))
  feat <- xy[, lc[length(lc)]]
  ring <- do.call(paste, c(as.data.frame(xy[, lc, drop = FALSE]), sep = "_"))
  kx <- round(xy[, 1] / tol)
  ky <- round(xy[, 2] / tol)
  vkey <- paste(format(kx, scientific = FALSE, trim = TRUE), format(ky, scientific = FALSE, trim = TRUE))
  if (contiguity == "queen") {
    item <- vkey
    own <- feat
  } else {
    nx <- length(vkey)
    same <- ring[-1] == ring[-nx]
    a <- vkey[-nx][same]
    b <- vkey[-1][same]
    keep <- a != b
    item <- ifelse(a < b, paste(a, b), paste(b, a))[keep]
    own <- feat[-nx][same][keep]
  }
  u <- unique(data.frame(item = item, own = own, stringsAsFactors = FALSE))
  grp <- split(u$own, u$item)
  grp <- grp[lengths(grp) > 1]
  ii <- unlist(lapply(grp, function(o) rep(o, each = length(o))), use.names = FALSE)
  jj <- unlist(lapply(grp, function(o) rep(o, times = length(o))), use.names = FALSE)
  cen <- matrix(NA_real_, n, 2)
  cs <- rowsum(xy[, 1:2, drop = FALSE], feat)
  cen[as.integer(rownames(cs)), ] <- cs / as.vector(table(feat)[rownames(cs)])
  g <- graph_edges(ii, jj, n, FALSE, sprintf("polygons, %s contiguity", contiguity), cen[, 1:2])
  g$geometry <- geom
  g
}

#' Graph of a spatial lag
#'
#' The sites at graph distance exactly \code{lag} from each site, in the
#' symmetrized graph. Used by \code{spatial_test(lags = L)}.
#'
#' @param graph An object of class \code{"ivspatial_graph"}.
#' @param lag A whole number of at least 1.
#' @return An object of class \code{"ivspatial_graph"}.
#' @examples
#' lag_graph(lattice_graph(10, 10), 2)
#' @export
lag_graph <- function(graph, lag) {
  check_graph(graph)
  if (lag < 1 || lag != round(lag)) stop("lag must be a whole number of at least 1", call. = FALSE)
  n <- graph$n
  A <- adjacency(graph)
  A <- (A + Matrix::t(A)) > 0
  A <- methods::as(A, "dMatrix")
  reach <- Matrix::Diagonal(n) > 0
  prev <- reach
  ring <- reach
  for (l in seq_len(lag)) {
    new <- (methods::as(prev, "dMatrix") %*% A) > 0
    ring <- new & !reach
    reach <- reach | new
    prev <- ring
  }
  ring <- methods::as(methods::as(ring, "dMatrix"), "TsparseMatrix")
  out <- graph_edges(ring@i + 1L, ring@j + 1L, n, FALSE, sprintf("%s, lag %d", graph$kind, as.integer(lag)),
                     graph$coords)
  out$geometry <- graph$geometry
  out
}

check_graph <- function(graph) {
  if (!inherits(graph, "ivspatial_graph"))
    stop("graph must be an ivspatial graph (see ?graphs, e.g. knn_graph() or lattice_graph())", call. = FALSE)
  invisible(TRUE)
}

# binary adjacency matrix A[v, u] = 1 if u is a neighbour of v
adjacency <- function(graph) {
  Matrix::sparseMatrix(i = rep(seq_len(graph$n), graph$degree), j = unlist(graph$neighbors),
                       x = 1, dims = c(graph$n, graph$n))
}

is_symmetric_graph <- function(graph) {
  A <- adjacency(graph)
  Matrix::nnzero(A - Matrix::t(A)) == 0
}

graph_summary <- function(graph) {
  deg <- graph$degree
  s <- sprintf("%s: n = %d sites, degree mean %.2f (min %d, max %d)%s", graph$kind, graph$n, mean(deg),
               min(deg), max(deg), if (is_symmetric_graph(graph)) "" else ", directed")
  iso <- sum(deg == 0)
  if (iso) s <- paste0(s, sprintf(", %d isolated site(s)", iso))
  s
}

#' @method format ivspatial_graph
#' @export
format.ivspatial_graph <- function(x, ...) graph_summary(x)

#' @method print ivspatial_graph
#' @export
print.ivspatial_graph <- function(x, ...) {
  cat("<ivspatial graph> ", graph_summary(x), "\n", sep = "")
  invisible(x)
}
