#' ivspatial: nonparametric tests of spatial independence
#'
#' Rank-based tests of the hypothesis that the values of a variable observed
#' at spatial units (grid cells, polygons, irregular points or the nodes of a
#' network) are independent and identically distributed, with the procedure of
#' Ashtari Nezhad, \emph{A nonparametric test of spatial independence on neighbourhood graphs}.
#'
#' Main functions:
#' \describe{
#'   \item{\code{\link{spatial_test}}}{the test for one variable;}
#'   \item{\code{\link{graphs}}, \code{\link{polygon_graph}}, \code{\link{build_graph}}}{neighbourhood graphs;}
#'   \item{\code{\link{screen_variables}}}{many variables on one graph, with Benjamini-Hochberg adjustment;}
#'   \item{\code{\link{residual_test}}}{residuals of a regression;}
#'   \item{\code{\link{local_features}}, \code{\link{map_local}}}{where the dependence is;}
#'   \item{\code{\link{null_limit}}}{the weighted chi-square null limit;}
#'   \item{\code{\link{read_sites}}, \code{\link{read_polygons}}}{reading data;}
#'   \item{\code{\link{ivspatial_app}}}{the graphical interface;}
#'   \item{\code{\link{ivspatial_guide}}}{the user guide in the browser.}
#' }
#' The user guide is \code{vignette("ivspatial")}. A Python version with the
#' same functions is available from PyPI (\code{pip install ivspatial}).
#'
#' @keywords internal
#' @import Matrix
#' @importFrom methods as
"_PACKAGE"

#' Mercer and Hall wheat yields
#'
#' Grain and straw yields (pounds) of 500 plots of wheat on a 20 x 25 grid,
#' from the uniformity trial of Mercer and Hall (1911); one of the
#' applications in the paper.
#'
#' The same data are in the file
#' \code{system.file("extdata", "mercer_hall_wheat.csv", package = "ivspatial")}.
#'
#' @format A data frame with 500 rows and 4 columns: \code{row} (1 to 20),
#'   \code{col} (1 to 25), \code{grain} and \code{straw}.
#' @source Mercer, W. B. and Hall, A. D. (1911). The experimental error of
#'   field trials. \emph{Journal of Agricultural Science}, 4, 107--132.
#' @examples
#' g <- lattice_from_indices(mercer_hall$row, mercer_hall$col)
#' spatial_test(mercer_hall$grain, g, B = 499, seed = 1)
"mercer_hall"

#' Baltimore house sales
#'
#' Prices and characteristics of 211 houses sold in Baltimore, Maryland, in
#' 1978, with the coordinates of the houses; one of the applications in the
#' paper.
#'
#' The same data are in the file
#' \code{system.file("extdata", "baltimore_houses.csv", package = "ivspatial")}.
#'
#' @format A data frame with 211 rows and 17 columns: \code{STATION} (id),
#'   \code{PRICE} (thousands of dollars), \code{NROOM}, \code{DWELL},
#'   \code{NBATH}, \code{PATIO}, \code{FIREPL}, \code{AC}, \code{BMENT},
#'   \code{NSTOR}, \code{GAR}, \code{AGE}, \code{CITCOU}, \code{LOTSZ},
#'   \code{SQFT} and the coordinates \code{X}, \code{Y}.
#' @source Dubin, R. A. (1992). Spatial autocorrelation and neighborhood
#'   quality. \emph{Regional Science and Urban Economics}, 22, 433--452; data
#'   distributed with the GeoDa Center sample data.
#' @examples
#' g <- knn_graph(cbind(baltimore$X, baltimore$Y), k = 4)
#' spatial_test(log(baltimore$PRICE), g, B = 499, seed = 1)
"baltimore"
