#' Open the user guide
#'
#' Opens the user guide of the package (the same text as
#' \code{vignette("ivspatial")}) in the web browser. The guide is installed
#' with the package, also when it is installed from GitHub without building
#' the vignettes.
#'
#' @return The path of the guide, invisibly.
#' @examples
#' if (interactive()) ivspatial_guide()
#' @export
ivspatial_guide <- function() {
  path <- system.file("guide", "ivspatial-guide.html", package = "ivspatial")
  if (!nzchar(path)) stop("the guide is not installed; see vignette(\"ivspatial\")", call. = FALSE)
  utils::browseURL(path)
  invisible(path)
}
