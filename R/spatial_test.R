# Nonparametric tests of spatial independence (Ashtari Nezhad, manuscript).

RESOLUTION_FAMILY <- rbind(cbind(d = 2:6, s = 2), cbind(d = 2:6, s = 3))

with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = genv) else assign(".Random.seed", old, envir = genv))
  set.seed(seed)
  expr
}

fmt_p <- function(p) ifelse(p < 1e-4, "< 0.0001", sprintf("%.4f", p))

# ----------------------------------------------------------------------
# components
# ----------------------------------------------------------------------
make_component <- function(name, kind, graph, d = 0, s = 0, lag = 1) {
  comp <- list(name = name, kind = kind, graph = graph, d = d, s = s, lag = lag)
  if (kind == "TV") {
    comp$eng <- tv_engine(graph, d, s)
    comp$dim <- length(comp$eng$feats)
  } else {
    comp$rm <- rank_moran(graph)
    comp$dim <- 1
  }
  comp
}

component_values <- function(comp, R) {
  R <- as.matrix(R)
  n <- nrow(R)
  if (comp$kind == "TV") return(tv_T(comp$eng, (comp$d * R) %/% n))
  rank_moran_stats(comp$rm, R)[[comp$kind]]
}

component_description <- function(comp) {
  lag <- if (comp$lag > 1) sprintf(", lag %d", comp$lag) else ""
  if (comp$kind == "TV") {
    order <- if (comp$s == 2) "2" else sprintf("2..%d", comp$s)
    return(sprintf("T_V(%d,%d)%s: %d features of order %s, %d classes", comp$d, comp$s, lag, comp$dim, order, comp$d))
  }
  paste0(comp$kind, lag, if (comp$kind == "RM") ": rank Moran (level)" else ": Moran of squared ranks (scale)")
}

build_components <- function(graph, spec, lags) {
  graphs <- c(list(graph), lapply(seq_len(max(lags - 1, 0)) + 1, function(l) lag_graph(graph, l)))
  comps <- list()
  for (l in seq_along(graphs)) {
    gr <- graphs[[l]]
    suf <- if (l == 1) "" else sprintf("@lag%d", l)
    if (spec == "C4") {
      comps <- c(comps, list(make_component(paste0("V43", suf), "TV", gr, 4, 3, l),
                             make_component(paste0("V42", suf), "TV", gr, 4, 2, l),
                             make_component(paste0("RM", suf), "RM", gr, lag = l),
                             make_component(paste0("RS", suf), "RS", gr, lag = l)))
    } else {
      for (r in seq_len(nrow(RESOLUTION_FAMILY))) {
        d <- RESOLUTION_FAMILY[r, "d"]
        s <- RESOLUTION_FAMILY[r, "s"]
        comps <- c(comps, list(make_component(sprintf("V%d%d%s", d, s, suf), "TV", gr, d, s, l)))
      }
    }
  }
  comps
}

#' Reference table of random rank arrangements
#'
#' The statistics of \code{B} random rank arrangements on a graph. The table
#' does not depend on the data, so one table serves every variable observed
#' on the same graph; \code{\link{spatial_test}} and
#' \code{\link{screen_variables}} build it when needed, and a table can be
#' passed to them with \code{reference} to save time.
#'
#' @param graph An object of class \code{"ivspatial_graph"}.
#' @param resolution \code{"C4"} (the four components of the combined test)
#'   or \code{"auto"}/\code{"schwarz"} (the family \eqn{T_V(d,s)},
#'   \eqn{d = 2..6}, \eqn{s = 2, 3}).
#' @param lags Number of spatial lags.
#' @param B Number of random rank arrangements.
#' @param seed Optional seed.
#' @param chunk Number of arrangements computed at once.
#' @param progress Optional function called with the fraction done.
#' @return An object of class \code{"ivspatial_reference"}.
#' @examples
#' g <- lattice_graph(20, 25)
#' ref <- reference_table(g, B = 199, seed = 1)
#' spatial_test(mercer_hall$grain, g, reference = ref)
#' @export
reference_table <- function(graph, resolution = c("C4", "auto", "schwarz"), lags = 1, B = 9999, seed = NULL,
                            chunk = 250, progress = NULL) {
  check_graph(graph)
  resolution <- match.arg(resolution)
  spec <- if (resolution == "C4") "C4" else "family"
  B <- as.integer(B)
  if (B < 1) stop("B must be at least 1", call. = FALSE)
  comps <- build_components(graph, spec, lags)
  n <- graph$n
  out <- matrix(0, B, length(comps), dimnames = list(NULL, vapply(comps, `[[`, "", "name")))
  with_seed(seed, {
    done <- 0L
    while (done < B) {
      b <- min(chunk, B - done)
      R <- vapply(seq_len(b), function(i) sample.int(n) - 1L, integer(n))
      if (!is.matrix(R)) R <- matrix(R, n)
      for (j in seq_along(comps)) out[done + seq_len(b), j] <- component_values(comps[[j]], R)
      done <- done + b
      if (is.function(progress)) progress(done / B)
    }
  })
  P <- out
  for (j in seq_len(ncol(out))) {
    s <- sort(out[, j])
    P[, j] <- (B - findInterval(out[, j] - 1e-9 * abs(out[, j]), s, left.open = TRUE)) / B
  }
  structure(list(graph = graph, spec = spec, lags = lags, B = B, comps = comps, ref = out, refP = P),
            class = "ivspatial_reference")
}

#' @method print ivspatial_reference
#' @export
print.ivspatial_reference <- function(x, ...) {
  cat(sprintf("<ivspatial reference table> %d arrangements of %s on %s\n", x$B,
              paste(colnames(x$ref), collapse = ", "), graph_summary(x$graph)))
  invisible(x)
}

ref_p_values <- function(ref, T) {
  vapply(seq_along(T), function(j) (1 + sum(ref$ref[, j] >= T[j] - 1e-9 * abs(T[j]))) / (ref$B + 1), numeric(1))
}

ref_p_min <- function(ref, p, cols) {
  mref <- apply(ref$refP[, cols, drop = FALSE], 1, min)
  (1 + sum(mref <= min(p[cols]) + 1e-12)) / (ref$B + 1)
}

# ----------------------------------------------------------------------
# main function
# ----------------------------------------------------------------------
ASY_SIGMA_LIMIT <- 6e6

standardized_features <- function(graph, r, cache) {
  if (is.null(cache$feat)) {
    eng <- tv_engine(graph, 3, 2)
    assign("feat", list(eng = eng, sd = sqrt(diag(tv_sigma(eng)))), envir = cache)
  }
  f <- cache$feat
  w <- tv_W(f$eng, matrix((3 * r) %/% length(r)))[, 1]
  stats::setNames(w / f$sd, f$eng$names)
}

asymptotic_p <- function(comp, T, cache, asymptotic) {
  if (isFALSE(asymptotic) || comp$kind != "TV") return(NA_real_)
  if (!identical(asymptotic, "all") && n_sets(comp$graph, comp$d, comp$s) > ASY_SIGMA_LIMIT) return(NA_real_)
  key <- paste("lam", comp$name)
  if (is.null(cache[[key]])) {
    lam <- eigen(tv_sigma(comp$eng), symmetric = TRUE, only.values = TRUE)$values
    assign(key, lam[lam > 1e-10], envir = cache)
  }
  weighted_chi2_sf(T, cache[[key]])
}

#' Nonparametric test of spatial independence
#'
#' Tests whether the values of a variable observed at the sites of a
#' neighbourhood graph are independent and identically distributed. The values
#' are replaced by their ranks and by quantile classes; the statistics compare
#' the classes of every site with those of its neighbours.
#'
#' The proposed test \eqn{C_4} is the minimum P-value of four components:
#' \describe{
#'   \item{\eqn{T_V(4,3)}}{features of order 2 and 3 with 4 quantile classes:
#'     dependence in the joint distribution of a site and two neighbours,
#'     including dependence that is invisible in pairs;}
#'   \item{\eqn{T_V(4,2)}}{pairs with resolved tails;}
#'   \item{RM}{the Moran coefficient of the ranks (level);}
#'   \item{RS}{the Moran coefficient of the squared standardized ranks
#'     (scale).}
#' }
#' All statistics depend on the data only through the ranks, so under
#' independence the permutation P-values are exact for every sample size,
#' every graph and every continuous distribution; ties are broken at random.
#' The minimum is calibrated exactly from the same table of random rank
#' arrangements (\code{\link{reference_table}}).
#'
#' With \code{resolution = "auto"} the family \eqn{T_V(d,s)}, \eqn{d = 2..6},
#' \eqn{s = 2, 3}, is combined by the minimum-P rule and the resolution chosen
#' by the data is reported; with \code{"schwarz"} the statistic maximizing
#' \eqn{T_V - p \log n} (\eqn{p} features) is used. In both cases the selection
#' is repeated in every rank arrangement, so the P-value is exact after the
#' selection.
#'
#' @param x Numeric vector: the variable at the \eqn{n} sites of
#'   \code{graph}, in the same order. Missing values are not allowed (remove
#'   the sites and build the graph for the remaining ones).
#' @param graph Neighbourhood graph (see \code{\link{graphs}} and
#'   \code{\link{polygon_graph}}).
#' @param B Number of random rank arrangements (default 9999). The P-values
#'   are exact for every \code{B}.
#' @param seed Optional seed for reproducible P-values. The random numbers of
#'   the session are not changed.
#' @param alpha Level of the test.
#' @param resolution \code{"C4"} (default), \code{"auto"} or
#'   \code{"schwarz"}.
#' @param lags 1 (default) or more: the same components on the graphs of
#'   sites at distance 2, ..., \code{lags} are added and all are combined by
#'   the minimum-P rule.
#' @param label Optional name of the variable, shown in the report.
#' @param features If \code{TRUE}, the standardized features of
#'   \eqn{T_V(3,2)} are reported.
#' @param asymptotic If \code{TRUE}, asymptotic P-values of the \eqn{T_V}
#'   components are reported (skipped for very large designs unless
#'   \code{"all"}).
#' @param reference Optional reference table from
#'   \code{\link{reference_table}}.
#' @param residuals Internal: description of residuals.
#' @return An object of class \code{"ivspatial"}: a list with \code{p_value}
#'   (the combined test), \code{reject}, \code{components} (a data frame with
#'   statistics, permutation and asymptotic P-values), \code{features},
#'   \code{selected} and the settings. \code{print} gives the report and
#'   \code{as.data.frame} one row of results.
#' @references Ashtari Nezhad, E. (2026). A nonparametric test of spatial independence on neighbourhood graphs. Manuscript.
#' @examples
#' g <- lattice_from_indices(mercer_hall$row, mercer_hall$col)
#' spatial_test(mercer_hall$grain, g, B = 499, seed = 1, label = "grain")
#'
#' b <- baltimore
#' gk <- knn_graph(cbind(b$X, b$Y), k = 4)
#' r <- spatial_test(log(b$PRICE), gk, B = 499, seed = 1)
#' r$p_value
#' r$components
#' @export
spatial_test <- function(x, graph, B = 9999, seed = NULL, alpha = 0.05, resolution = c("C4", "auto", "schwarz"),
                         lags = 1, label = "", features = TRUE, asymptotic = TRUE, reference = NULL,
                         residuals = "") {
  spatial_test_impl(x, graph, B, seed, alpha, match.arg(resolution), lags, label, features, asymptotic, reference,
                    new.env(), residuals)
}

spatial_test_impl <- function(x, graph, B, seed, alpha, resolution, lags, label, features, asymptotic, reference,
                              cache, residuals) {
  check_graph(graph)
  x <- as.numeric(x)
  n <- graph$n
  if (length(x) != n) stop(sprintf("the variable has %d values but the graph has %d sites", length(x), n),
                           call. = FALSE)
  if (any(!is.finite(x)))
    stop(sprintf("%d missing or non-finite values: remove those sites and rebuild the graph", sum(!is.finite(x))),
         call. = FALSE)
  if (n < 20) stop("at least 20 sites are needed", call. = FALSE)
  if (length(unique(x)) < 3) stop("the variable is (almost) constant", call. = FALSE)
  if (lags < 1 || lags != round(lags)) stop("lags must be a whole number of at least 1", call. = FALSE)
  if (min(graph$degree) == 0) warning(sprintf("%d isolated site(s) contribute nothing", sum(graph$degree == 0)),
                                      call. = FALSE)
  spec <- if (resolution == "C4") "C4" else "family"
  ref <- reference
  if (is.null(ref) || ref$spec != spec || ref$lags != lags || !identical(ref$graph$neighbors, graph$neighbors))
    ref <- reference_table(graph, resolution, lags, B, seed)
  r <- with_seed(if (is.null(seed)) NULL else seed + 7919, rank0(x))
  T <- vapply(ref$comps, function(cp) component_values(cp, matrix(r)), numeric(1))
  p <- ref_p_values(ref, T)
  pa <- vapply(seq_along(ref$comps), function(j) asymptotic_p(ref$comps[[j]], T[j], cache, asymptotic), numeric(1))
  names(T) <- vapply(ref$comps, `[[`, "", "name")
  comps <- data.frame(name = names(T),
                      description = vapply(ref$comps, component_description, ""),
                      statistic = T, p_value = p, p_asymptotic = pa, stringsAsFactors = FALSE)
  rownames(comps) <- NULL
  selected <- NULL
  lagtxt <- if (lags > 1) sprintf(", lags 1..%d", lags) else ""
  if (resolution %in% c("C4", "auto")) {
    pC <- ref_p_min(ref, p, seq_along(p))
    method <- paste0(if (resolution == "C4") "C4 = min-P(V43, V42, RM, RS)" else "min-P over T_V(d,s), d=2..6, s=2,3",
                     lagtxt)
    if (resolution == "auto") {
      cand <- which(p <= min(p) + 1e-12)
      z <- (T[cand] - colMeans(ref$ref[, cand, drop = FALSE])) /
        apply(ref$ref[, cand, drop = FALSE], 2, function(v) sqrt(mean((v - mean(v))^2)))
      j <- cand[which.max(z)]
      selected <- sub(":.*$", "", comps$description[j])
    }
  } else {
    dims <- vapply(ref$comps, `[[`, 0, "dim")
    pen <- log(n)
    sel <- function(Tv) {
      crit <- Tv - dims * pen
      best <- which(crit >= max(crit) - 1e-12)
      best[which.min(dims[best])]
    }
    refsel <- vapply(seq_len(ref$B), function(b) ref$ref[b, sel(ref$ref[b, ])], numeric(1))
    j <- sel(T)
    pC <- (1 + sum(refsel >= T[j] - 1e-9 * abs(T[j]))) / (ref$B + 1)
    selected <- sub(":.*$", "", comps$description[j])
    method <- paste0("Schwarz-selected T_V", lagtxt)
  }
  feats <- if (isTRUE(features)) standardized_features(graph, r, cache) else numeric()
  structure(list(label = label, n = n, graph = graph_summary(graph), B = ref$B, alpha = alpha, method = method,
                 p_value = pC, reject = pC <= alpha, components = comps, features = feats, selected = selected,
                 residuals = residuals, resolution = resolution, lags = lags),
            class = "ivspatial")
}

FEATURE_EXPLANATION <- c("c.lin|n.lin" = "level: site and neighbours similar",
                         "c.quad|n.quad" = "scale: extremes next to extremes",
                         "n.lin^2" = "level: neighbours similar to each other",
                         "n.quad^2" = "scale: neighbours jointly extreme",
                         "c.lin|n.quad" = "low/high site among extreme neighbours",
                         "c.quad|n.lin" = "extreme site among low/high neighbours",
                         "n.lin*quad" = "mixed level-scale between neighbours")

#' Report of a spatial independence test
#'
#' @param x An object of class \code{"ivspatial"}.
#' @param row.names,optional Ignored.
#' @param ... Ignored.
#' @return \code{format}: the lines of the report; \code{print}: \code{x},
#'   invisibly; \code{as.data.frame}: a data frame with one row.
#' @examples
#' g <- lattice_graph(20, 25)
#' r <- spatial_test(mercer_hall$straw, g, B = 199, seed = 1)
#' print(r)
#' as.data.frame(r)
#' @method format ivspatial
#' @export
format.ivspatial <- function(x, ...) {
  title <- paste0("Nonparametric test of spatial independence", if (nzchar(x$label)) sprintf("  [%s]", x$label) else "")
  L <- c(title, strrep("=", nchar(title)), sprintf("Observations n = %d", x$n))
  if (nzchar(x$residuals)) L <- c(L, paste("Data:", x$residuals))
  L <- c(L, paste("Neighbourhood graph:", x$graph), sprintf("Random rank arrangements B = %d", x$B), "",
         sprintf("%-56s%11s%11s%11s", "Test", "Statistic", "P-value", "asympt."),
         sprintf("%-56s%11s%11s%11s", paste0(x$method, "  [proposed]"), "", fmt_p(x$p_value), ""))
  cp <- x$components
  for (i in seq_len(nrow(cp))) {
    pa <- if (is.na(cp$p_asymptotic[i])) "" else fmt_p(cp$p_asymptotic[i])
    L <- c(L, sprintf("%-56s%11.3f%11s%11s", paste0("  ", substr(cp$description[i], 1, 54)), cp$statistic[i],
                      fmt_p(cp$p_value[i]), pa))
  }
  L <- c(L, "", sprintf("Spatial independence is %s at the %g level (%s).",
                        if (x$reject) "REJECTED" else "not rejected", x$alpha, x$method))
  if (!is.null(x$selected))
    L <- c(L, sprintf("Resolution selected by the data: %s (exactly calibrated after selection).", x$selected))
  if (length(x$features)) {
    L <- c(L, "", "Standardized features of T_V(3,2) (approximately N(0,1) under independence):")
    ex <- FEATURE_EXPLANATION[names(x$features)]
    ex[is.na(ex)] <- ""
    L <- c(L, sprintf("  %-16s%8.2f   %s", names(x$features), x$features, ex),
           "  |value| > 3 is notable; the sign of mixed features depends on the coding direction.")
  }
  L
}

#' @rdname format.ivspatial
#' @method print ivspatial
#' @export
print.ivspatial <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

#' @rdname format.ivspatial
#' @method as.data.frame ivspatial
#' @export
as.data.frame.ivspatial <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- data.frame(variable = x$label, n = x$n, method = x$method, p_value = x$p_value, reject = x$reject,
                    selected = if (is.null(x$selected)) NA_character_ else x$selected, stringsAsFactors = FALSE)
  for (i in seq_len(nrow(x$components))) {
    nm <- x$components$name[i]
    out[[paste0("T_", nm)]] <- x$components$statistic[i]
    out[[paste0("p_", nm)]] <- x$components$p_value[i]
  }
  for (k in names(x$features)) out[[paste0("feature ", k)]] <- x$features[[k]]
  out
}

# ----------------------------------------------------------------------
# several variables, residuals, local features, null limit
# ----------------------------------------------------------------------

#' Test many variables on one graph (screening)
#'
#' Applies \code{\link{spatial_test}} to every numeric column of a data frame
#' (or the columns named in \code{columns}), with one reference table shared
#' by all variables, and adjusts the P-values by the method of Benjamini and
#' Hochberg.
#'
#' @param data A data frame (one row per site, in the order of the graph).
#' @param graph Neighbourhood graph.
#' @param columns Optional names of the columns to test (default: every
#'   numeric column).
#' @param B,seed,alpha,resolution,lags,features,asymptotic Options of
#'   \code{\link{spatial_test}}.
#' @return An object of class \code{"ivspatial_screen"}: a list with
#'   \code{table} (one row per variable with the adjusted P-value and the
#'   decision), \code{results} and \code{skipped}.
#' @examples
#' g <- lattice_graph(20, 25)
#' screen_variables(mercer_hall, g, columns = c("grain", "straw"), B = 199, seed = 1)
#' @export
screen_variables <- function(data, graph, columns = NULL, B = 9999, seed = NULL, alpha = 0.05,
                             resolution = c("C4", "auto", "schwarz"), lags = 1, features = TRUE, asymptotic = TRUE) {
  resolution <- match.arg(resolution)
  check_graph(graph)
  data <- as.data.frame(data)
  num <- names(data)[vapply(data, is.numeric, logical(1))]
  columns <- if (is.null(columns)) num else columns
  miss <- setdiff(columns, num)
  if (length(miss)) stop("column(s) not found or not numeric: ", paste(miss, collapse = ", "), call. = FALSE)
  if (nrow(data) != graph$n) stop(sprintf("the data have %d rows but the graph has %d sites", nrow(data), graph$n),
                                  call. = FALSE)
  ref <- reference_table(graph, resolution, lags, B, seed)
  cache <- new.env()
  results <- list()
  skipped <- character()
  for (nm in columns) {
    r <- tryCatch(suppressWarnings(spatial_test_impl(data[[nm]], graph, B, seed, alpha, resolution, lags, nm,
                                                     features, asymptotic, ref, cache, "")),
                  error = function(e) conditionMessage(e))
    if (is.character(r)) skipped[nm] <- r else results[[nm]] <- r
  }
  tab <- NULL
  if (length(results)) {
    tab <- data.frame(variable = names(results), n = vapply(results, `[[`, 0, "n"),
                      p_value = vapply(results, `[[`, 0, "p_value"), stringsAsFactors = FALSE)
    tab$p_adj <- stats::p.adjust(tab$p_value, "BH")
    tab$reject <- tab$p_adj <= alpha
    rownames(tab) <- NULL
  }
  structure(list(table = tab, results = results, skipped = skipped, alpha = alpha), class = "ivspatial_screen")
}

#' @method format ivspatial_screen
#' @export
format.ivspatial_screen <- function(x, ...) {
  L <- "Screening summary (Benjamini-Hochberg adjustment over the variables tested)"
  t <- x$table
  if (!is.null(t))
    L <- c(L, sprintf("%-24s%7s%10s%10s  %s", "variable", "n", "p", "adj. p", "decision"),
           sprintf("%-24s%7d%10.4f%10.4f  %s", substr(t$variable, 1, 23), as.integer(t$n), t$p_value, t$p_adj,
                   ifelse(t$reject, "reject", "-")))
  for (nm in names(x$skipped)) L <- c(L, sprintf("[%s] skipped: %s", nm, x$skipped[[nm]]))
  L
}

#' @method print ivspatial_screen
#' @export
print.ivspatial_screen <- function(x, ...) {
  cat(format(x), sep = "\n")
  invisible(x)
}

#' @method summary ivspatial_screen
#' @export
summary.ivspatial_screen <- function(object, ...) {
  for (r in object$results) {
    print(r)
    cat("\n")
  }
  print(object)
  invisible(object)
}

#' Test the residuals of a regression
#'
#' Tests the ordinary least squares residuals of \code{y} on the columns of
#' \code{X}. The test is valid without correction when the regressors are
#' exogenous (the residuals are asymptotically exchangeable). It is not valid
#' for spatial-lag models, in which \code{y} is regressed on its own spatial
#' lag; use a model-based bootstrap there.
#'
#' @param y Response at the sites.
#' @param X Matrix or data frame of regressors (one row per site).
#' @param graph Neighbourhood graph.
#' @param intercept If \code{TRUE} (default), an intercept is added.
#' @param ... Options of \code{\link{spatial_test}}.
#' @return An object of class \code{"ivspatial"}.
#' @examples
#' b <- baltimore
#' g <- knn_graph(cbind(b$X, b$Y), k = 4)
#' residual_test(log(b$PRICE), b[, c("NROOM", "NBATH", "SQFT", "AGE")], g, B = 499, seed = 1)
#' @export
residual_test <- function(y, X, graph, intercept = TRUE, ...) {
  y <- as.numeric(y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != length(y)) stop("X must have one row per site", call. = FALSE)
  k <- ncol(X)
  if (intercept) X <- cbind(1, X)
  fit <- stats::lm.fit(X, y)
  e <- fit$residuals
  r2 <- 1 - mean((e - mean(e))^2) / mean((y - mean(y))^2)
  spatial_test(e, graph, residuals = sprintf("OLS residuals on %d regressor(s)%s, R^2 = %.3f", k,
                                              if (intercept) " and an intercept" else "", r2), ...)
}

#' Local feature contributions
#'
#' The contributions \eqn{f_{v,\tau}} of the features of \eqn{T_V(d,2)} at
#' every site; under independence each has mean 0 and variance 1, and their
#' scaled sum is the feature \eqn{W_\tau}. Like local indicators of spatial
#' association they show where a feature receives its contribution, and are
#' meant for maps, not as tests.
#'
#' @param x The variable at the sites.
#' @param graph Neighbourhood graph.
#' @param d Number of quantile classes (default 4).
#' @param seed Optional seed (ties are broken at random).
#' @return A data frame with one row per site and one column per feature.
#' @examples
#' g <- lattice_graph(20, 25)
#' lf <- local_features(mercer_hall$grain, g)
#' head(lf)
#' @export
local_features <- function(x, graph, d = 4, seed = NULL) {
  check_graph(graph)
  x <- as.numeric(x)
  if (length(x) != graph$n || any(!is.finite(x))) stop("x must have one finite value per site", call. = FALSE)
  r <- with_seed(seed, rank0(x))
  eng <- tv_engine(graph, d, 2)
  L <- tv_local(eng, matrix((d * r) %/% graph$n))
  out <- as.data.frame(t(matrix(L, dim(L)[1], dim(L)[2])))
  names(out) <- eng$names
  out
}

#' Null limit of the statistic T_V(d, s) on a graph
#'
#' Under independence \eqn{T_V(d,s)} converges to \eqn{\sum_i \lambda_i
#' \chi^2_{1,i}}, where \eqn{\lambda_i} are the eigenvalues of the exact null
#' covariance matrix of the features. \code{null_limit} gives the weights,
#' mean, variance, critical values at the levels 0.10, 0.05 and 0.01, and the
#' asymptotic P-value of \code{value}; tail probabilities are computed by
#' numerical inversion of the characteristic function (Imhof, 1961).
#'
#' @param graph Neighbourhood graph.
#' @param d Number of quantile classes.
#' @param s Highest order of the features (2 or 3).
#' @param value Optional value of the statistic.
#' @return An object of class \code{"ivspatial_null"}.
#' @references Imhof, J. P. (1961). Computing the distribution of quadratic
#'   forms in normal variables. \emph{Biometrika}, 48, 419--426.
#' @examples
#' null_limit(lattice_graph(20, 25), d = 4, s = 2, value = 35)
#' @export
null_limit <- function(graph, d = 4, s = 2, value = NULL) {
  check_graph(graph)
  eng <- tv_engine(graph, d, s)
  lam <- eigen(tv_sigma(eng), symmetric = TRUE, only.values = TRUE)$values
  lam <- lam[lam > 1e-10]
  out <- list(d = d, s = s, graph = graph_summary(graph), features = length(eng$feats), eigenvalues = lam,
              mean = sum(lam), variance = 2 * sum(lam^2),
              critical = vapply(c(0.10, 0.05, 0.01), function(a) weighted_chi2_quantile(a, lam), numeric(1)),
              value = value, p_value = if (is.null(value)) NULL else weighted_chi2_sf(as.numeric(value), lam))
  names(out$critical) <- c("0.10", "0.05", "0.01")
  structure(out, class = "ivspatial_null")
}

#' @method print ivspatial_null
#' @export
print.ivspatial_null <- function(x, ...) {
  cat(x$graph, "\n", sep = "")
  cat(sprintf("T_V(%d,%d): %d features; limit = weighted sum of %d chi2(1)\n", x$d, x$s, x$features,
              length(x$eigenvalues)))
  cat(sprintf("mean = %.3f, variance = %.3f\n", x$mean, x$variance))
  for (k in names(x$critical)) cat(sprintf("critical value at level %s: %.3f\n", k, x$critical[[k]]))
  if (!is.null(x$value)) cat(sprintf("P(limit > %g) = %.5f\n", x$value, x$p_value))
  invisible(x)
}
