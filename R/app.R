#' Graphical interface
#'
#' Opens the graphical interface of the package in the web browser (or the
#' RStudio viewer); it is built with \pkg{shiny}, which is installed with the
#' package. Everything available in code is also available in the window:
#' \describe{
#'   \item{Test}{open a CSV, TXT or Excel file with one row per site, or a
#'     polygon file (GeoJSON or GeoPackage, or the \code{.shp}, \code{.dbf}
#'     and \code{.shx} files of a shapefile selected together), or use one of
#'     the example data sets; build the neighbourhood graph; choose a
#'     variable (or tick \emph{test all variables}), optionally a log
#'     transform and regressors whose residuals are tested; set the options,
#'     run the test and download the report (text) or the results (CSV).}
#'   \item{Map}{maps of the variable and of its local level and scale
#'     features.}
#'   \item{Null distribution}{the weighted chi-square null limit of
#'     \eqn{T_V(d,s)} on the current graph, its critical values and the
#'     P-value of any value of the statistic.}
#'   \item{User guide}{the user guide of the package.}
#' }
#' Close the browser tab and press Esc (or Stop) in R to finish.
#'
#' @param launch.browser Passed to \code{shiny::runApp}.
#' @return Called for its side effect; returns when the app is closed.
#' @examples
#' if (interactive()) ivspatial_app()
#' @export
ivspatial_app <- function(launch.browser = TRUE) {
  shiny::runApp(ivspatial_app_object(), launch.browser = launch.browser)
}

ivspatial_app_object <- function() {
  num <- function(id, label, value, width = "110px", ...) shiny::numericInput(id, label, value, width = width, ...)
  txt <- function(id, label, value = "", width = "110px") shiny::textInput(id, label, value, width = width)
  sel <- function(id, label, width = "160px", ...) shiny::selectInput(id, label, choices = character(), width = width, ...)
  row <- function(...) shiny::div(style = "display:flex;flex-wrap:wrap;gap:14px;align-items:flex-end;", ...)
  box <- function(title, ...) shiny::div(style = "border:1px solid #ddd;border-radius:6px;padding:8px 12px;margin:8px 0;",
                                         shiny::tags$b(title), ...)
  cond <- function(types, ...) shiny::conditionalPanel(
    sprintf("[%s].indexOf(input.gtype) >= 0", paste0("'", types, "'", collapse = ",")), ...)
  guide <- system.file("guide", "ivspatial-guide.html", package = "ivspatial")
  if (nzchar(guide)) shiny::addResourcePath("ivspatialdoc", dirname(guide))
  has_sf <- requireNamespace("sf", quietly = TRUE)

  ui <- shiny::fluidPage(
    shiny::tags$head(shiny::tags$style(shiny::HTML(
      "pre{font-size:12.5px;white-space:pre;overflow-x:auto} .help{color:#555;max-width:900px}
       .form-group{margin-bottom:6px} .shiny-input-container{margin-bottom:4px}"))),
    shiny::titlePanel(paste0("ivspatial ", utils::packageVersion("ivspatial"),
                             " - nonparametric test of spatial independence")),
    shiny::tabsetPanel(
      id = "tabs",
      shiny::tabPanel(
        "Test",
        box("1. Data",
            shiny::fileInput("file", "Data file: CSV, TXT, Excel (one row per site), or polygons",
                             multiple = TRUE, width = "420px",
                             accept = c(".csv", ".txt", ".dat", ".xlsx", ".xls", ".geojson", ".json", ".gpkg",
                                        ".shp", ".dbf", ".shx", ".prj", ".cpg")),
            row(shiny::actionButton("ex_mh", "Example: wheat yields (grid)"),
                shiny::actionButton("ex_balt", "Example: house prices (points)"),
                if (has_sf) shiny::actionButton("ex_nc", "Example: North Carolina counties (polygons)")),
            shiny::p(class = "help", "For a shapefile select the .shp, .dbf and .shx files together."),
            shiny::textOutput("fileinfo")),
        box("2. Neighbourhood graph",
            row(shiny::selectInput("gtype", "Graph", width = "300px",
                                   choices = c("k nearest neighbours (points)" = "knn",
                                               "Delaunay triangulation (points)" = "delaunay",
                                               "distance band (points)" = "band",
                                               "lattice from row/column indices" = "lattice",
                                               "grid stored row by row" = "grid",
                                               "polygon contiguity" = "polygons",
                                               "edge list file" = "edges")),
                cond(c("knn", "delaunay", "band"), row(sel("gx", "x coordinate"), sel("gy", "y coordinate"))),
                cond("knn", num("gk", "k", 4, min = 1)),
                cond("band", num("gradius", "distance", 1, min = 0)),
                cond("lattice", row(sel("grow", "row index"), sel("gcol", "column index"))),
                cond("grid", row(num("gnrow", "rows", 20, min = 1), num("gncol", "columns", 25, min = 1))),
                cond(c("lattice", "grid", "polygons"),
                     shiny::radioButtons("gcont", "Contiguity", c("rook", "queen"), inline = TRUE)),
                cond("edges", shiny::fileInput("gedges", "Edge list (CSV: two columns of site numbers from 1)",
                                               width = "320px")),
                shiny::actionButton("gbuild", "Build graph", class = "btn-info")),
            shiny::textOutput("graphinfo")),
        box("3. Variable and options",
            row(sel("column", "Variable", "220px"),
                shiny::checkboxInput("all", "test all variables (screening)", FALSE),
                shiny::radioButtons("transform", "Transform", c("none" = "none", "logarithm" = "log"), inline = TRUE),
                shiny::selectInput("regs", "Test residuals on regressors (optional)", choices = character(),
                                   multiple = TRUE, width = "320px")),
            row(num("B", "Permutations B", 9999, min = 1),
                txt("seed", "Seed (blank: random)", "1"),
                num("alpha", "Level", 0.05, min = 0.0001, max = 0.5, step = 0.01),
                num("lags", "Spatial lags", 1, min = 1, max = 4),
                shiny::selectInput("resolution", "Resolution", width = "260px",
                                   choices = c("C4 (proposed)" = "C4", "data-driven: minimum P" = "auto",
                                               "data-driven: Schwarz" = "schwarz")))),
        row(shiny::actionButton("run", "Run test", class = "btn-primary"),
            shiny::downloadButton("save_txt", "Save report (.txt)"),
            shiny::downloadButton("save_csv", "Save results (.csv)"),
            shiny::downloadButton("save_local", "Save local features (.csv)")),
        shiny::br(),
        shiny::verbatimTextOutput("report")
      ),
      shiny::tabPanel(
        "Map",
        shiny::br(),
        shiny::p(class = "help",
                 "Ranks of the last variable tested and its local level and scale features (contributions of the",
                 "features c.lin|n.lin and c.quad|n.quad of T_V(4,2) at every site; mean 0 and variance 1 under",
                 "independence). Like local indicators of spatial association they show where the dependence is;",
                 "they are descriptive, not tests."),
        shiny::plotOutput("map", height = "460px"),
        shiny::downloadButton("save_map", "Save map (.png)")
      ),
      shiny::tabPanel(
        "Null distribution",
        shiny::br(),
        shiny::p(class = "help",
                 "Under independence T_V(d, s) converges to a weighted sum of chi-square variables with one degree",
                 "of freedom, whose weights are the eigenvalues of the exact covariance matrix of the features on",
                 "the current graph. Enter d, s and, optionally, a value of the statistic for its asymptotic P-value."),
        row(num("n_d", "d (classes)", 4, min = 2, max = 8), num("n_s", "s (order)", 2, min = 2, max = 3),
            txt("n_q", "T_V (optional)", ""),
            shiny::actionButton("n_go", "Compute", class = "btn-primary")),
        shiny::br(),
        shiny::verbatimTextOutput("null")
      ),
      shiny::tabPanel(
        "User guide",
        if (nzchar(guide))
          shiny::tags$iframe(src = "ivspatialdoc/ivspatial-guide.html", style = "width:100%;height:80vh;border:0;")
        else shiny::p("The guide is available with vignette(\"ivspatial\") or ivspatial_guide().")
      )
    )
  )

  server <- function(input, output, session) {
    st <- shiny::reactiveValues(data = NULL, sfdata = NULL, name = "", graph = NULL,
                                report = "Open a data file (or an example), build the graph and press 'Run test'.",
                                results = NULL, last_x = NULL, last_label = "", null = "")

    numeric_cols <- function(df) names(df)[vapply(df, is.numeric, logical(1))]
    set_data <- function(df, name, sfdata = NULL, gtype = NULL, x = NULL, y = NULL, column = NULL, rc = NULL) {
      st$data <- df
      st$sfdata <- sfdata
      st$name <- name
      st$graph <- NULL
      cols <- numeric_cols(df)
      pick <- function(cands, default) {
        m <- intersect(cands, cols)
        if (length(m)) m[1] else default
      }
      shiny::updateSelectInput(session, "gx", choices = cols, selected = x %||% pick(c("x", "X", "lon", "long", "easting"), cols[1]))
      shiny::updateSelectInput(session, "gy", choices = cols,
                               selected = y %||% pick(c("y", "Y", "lat", "northing"), cols[min(2, length(cols))]))
      shiny::updateSelectInput(session, "grow", choices = cols, selected = rc[1] %||% pick(c("row", "ROW", "i"), cols[1]))
      shiny::updateSelectInput(session, "gcol", choices = cols,
                               selected = rc[2] %||% pick(c("col", "column", "COL", "j"), cols[min(2, length(cols))]))
      shiny::updateSelectInput(session, "column", choices = cols, selected = column %||% cols[1])
      shiny::updateSelectInput(session, "regs", choices = cols, selected = character())
      if (!is.null(gtype)) shiny::updateSelectInput(session, "gtype", selected = gtype)
      else shiny::updateSelectInput(session, "gtype", selected = if (is.null(sfdata)) "knn" else "polygons")
      if (!is.null(sfdata) || identical(gtype, "polygons")) shiny::updateRadioButtons(session, "gcont", selected = "queen")
      else shiny::updateRadioButtons(session, "gcont", selected = "rook")
    }
    shiny::observeEvent(input$file, {
      f <- input$file
      st$report <- tryCatch({
        ext <- tolower(tools::file_ext(f$name))
        if (any(ext %in% c("shp", "geojson", "json", "gpkg"))) {
          dir <- tempfile("ivspatial")
          dir.create(dir)
          file.copy(f$datapath, file.path(dir, f$name))
          main <- f$name[ext %in% c("shp", "geojson", "json", "gpkg")][1]
          if (tolower(tools::file_ext(main)) == "shp" && !all(c("dbf", "shx") %in% ext))
            stop("select the .shp, .dbf and .shx files of the shapefile together", call. = FALSE)
          sfd <- read_polygons(file.path(dir, main))
          set_data(sf::st_drop_geometry(sfd), main, sfdata = sfd)
        } else {
          set_data(read_sites(f$datapath[1]), f$name[1])
        }
        "Data read. Build the graph (step 2) and press 'Run test'."
      }, error = function(e) paste("Error:", conditionMessage(e)))
    })
    shiny::observeEvent(input$ex_mh, {
      set_data(read_sites(system.file("extdata", "mercer_hall_wheat.csv", package = "ivspatial")),
               "mercer_hall_wheat.csv (example: 500 plots on a 20 x 25 grid)", gtype = "lattice",
               column = "grain", rc = c("row", "col"))
      st$report <- "Example loaded: lattice graph from the columns row and col. Press 'Build graph', then 'Run test'."
    })
    shiny::observeEvent(input$ex_balt, {
      set_data(read_sites(system.file("extdata", "baltimore_houses.csv", package = "ivspatial")),
               "baltimore_houses.csv (example: 211 house sales)", gtype = "knn", x = "X", y = "Y", column = "PRICE")
      shiny::updateRadioButtons(session, "transform", selected = "log")
      st$report <- "Example loaded: 4-nearest-neighbour graph from X and Y. Press 'Build graph', then 'Run test'."
    })
    shiny::observeEvent(input$ex_nc, {
      sfd <- read_polygons(system.file("shape/nc.shp", package = "sf"))
      sfd$SIDR79 <- 1000 * sfd$SID79 / sfd$BIR79
      set_data(sf::st_drop_geometry(sfd), "nc.shp (example: 100 North Carolina counties; SIDR79 = SIDS per 1000 births)",
               sfdata = sfd, gtype = "polygons", column = "SIDR79")
      st$report <- "Example loaded: queen contiguity of the counties. Press 'Build graph', then 'Run test'."
    })
    output$fileinfo <- shiny::renderText({
      if (is.null(st$data)) "" else sprintf("%s: %d sites, %d numeric columns%s", st$name, nrow(st$data),
                                            length(numeric_cols(st$data)), if (is.null(st$sfdata)) "" else ", polygons")
    })

    make_graph <- function() {
      if (is.null(st$data)) stop("open a data file first", call. = FALSE)
      edges <- if (input$gtype == "edges") {
        if (is.null(input$gedges)) stop("choose the edge list file", call. = FALSE)
        input$gedges$datapath
      }
      build_graph(if (input$gtype == "polygons") st$sfdata %||% st$data else st$data, input$gtype,
                  x = input$gx, y = input$gy, k = input$gk, radius = input$gradius, row = input$grow, col = input$gcol,
                  nrow = input$gnrow, ncol = input$gncol, contiguity = input$gcont, edges = edges)
    }
    shiny::observeEvent(input$gbuild, {
      st$graph <- tryCatch(make_graph(), error = function(e) {
        st$report <- paste("Error:", conditionMessage(e))
        NULL
      })
    })
    output$graphinfo <- shiny::renderText(if (is.null(st$graph)) "No graph built yet." else graph_summary(st$graph))

    get_var <- function(nm) {
      v <- as.numeric(st$data[[nm]])
      if (input$transform == "log") {
        if (any(v[is.finite(v)] <= 0)) stop(sprintf("the logarithm needs positive values ('%s')", nm), call. = FALSE)
        v <- log(v)
      }
      v
    }

    shiny::observeEvent(input$run, {
      out <- tryCatch({
        if (is.null(st$data)) stop("open a data file first", call. = FALSE)
        if (is.null(st$graph)) st$graph <- make_graph()
        g <- st$graph
        seed <- if (trimws(input$seed) == "") NULL else as.integer(input$seed)
        args <- list(B = input$B, seed = seed, alpha = input$alpha, resolution = input$resolution, lags = input$lags)
        msgs <- character()
        keep <- function(w) {
          msgs <<- c(msgs, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
        shiny::withProgress(message = "Running random rank arrangements...", value = 0.5, {
          regs <- setdiff(input$regs, input$column)
          skip <- c(if (input$gtype %in% c("knn", "delaunay", "band")) c(input$gx, input$gy),
                    if (input$gtype == "lattice") c(input$grow, input$gcol))
          if (isTRUE(input$all)) {
            cols <- setdiff(numeric_cols(st$data), c(skip, regs))
            if (length(regs)) {
              X <- st$data[, regs, drop = FALSE]
              res <- lapply(cols, function(nm) tryCatch(withCallingHandlers(
                do.call(residual_test, c(list(get_var(nm), X, g, label = nm), args)), warning = keep),
                error = function(e) NULL))
              names(res) <- cols
              res <- Filter(Negate(is.null), res)
              tab <- data.frame(variable = names(res), p = vapply(res, `[[`, 0, "p_value"))
              tab$adj <- stats::p.adjust(tab$p, "BH")
              scr <- c("Screening summary (Benjamini-Hochberg adjustment over the variables tested)",
                       sprintf("%-24s%10s%10s  %s", "variable", "p", "adj. p", "decision"),
                       sprintf("%-24s%10.4f%10.4f  %s", substr(tab$variable, 1, 23), tab$p, tab$adj,
                               ifelse(tab$adj <= input$alpha, "reject", "-")))
            } else {
              dat <- as.data.frame(lapply(stats::setNames(cols, cols), get_var), check.names = FALSE)
              sc <- withCallingHandlers(do.call(screen_variables, c(list(dat, g), args)), warning = keep)
              res <- sc$results
              scr <- format(sc)
            }
            st$results <- res
            st$last_x <- get_var(cols[1])
            st$last_label <- cols[1]
            txt <- c(graph_summary(g), "", unlist(lapply(res, function(r) c(format(r), ""))), scr)
          } else {
            x <- get_var(input$column)
            r <- withCallingHandlers(
              if (length(regs)) do.call(residual_test, c(list(x, st$data[, regs, drop = FALSE], g, label = input$column), args))
              else do.call(spatial_test, c(list(x, g, label = input$column), args)), warning = keep)
            st$results <- list(r)
            st$last_x <- if (length(regs)) {
              X <- cbind(1, as.matrix(st$data[, regs, drop = FALSE]))
              stats::lm.fit(X, x)$residuals
            } else x
            st$last_label <- paste0(input$column, if (length(regs)) " (residuals)" else "")
            txt <- format(r)
          }
          paste(c(txt, if (length(msgs)) paste("Note:", unique(msgs))), collapse = "\n")
        })
      }, error = function(e) paste("Error:", conditionMessage(e)))
      st$report <- out
    })
    output$report <- shiny::renderText(st$report)
    output$save_txt <- shiny::downloadHandler(
      filename = function() "ivspatial_report.txt",
      content = function(file) writeLines(st$report, file))
    output$save_csv <- shiny::downloadHandler(
      filename = function() "ivspatial_results.csv",
      content = function(file) {
        tab <- if (is.null(st$results)) data.frame() else do.call(rbind, lapply(st$results, as.data.frame))
        if (nrow(tab) > 1) tab$p_adj_BH <- stats::p.adjust(tab$p_value, "BH")
        utils::write.csv(tab, file, row.names = FALSE)
      })
    output$save_local <- shiny::downloadHandler(
      filename = function() "ivspatial_local_features.csv",
      content = function(file) {
        if (is.null(st$last_x)) {
          writeLines("Run a test first.", file)
          return()
        }
        lf <- local_features(st$last_x, st$graph)
        utils::write.csv(cbind(site = seq_len(nrow(lf)), lf), file, row.names = FALSE)
      })

    draw_map <- function() {
      if (is.null(st$last_x) || is.null(st$graph)) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, "Run a test first.")
        return(invisible())
      }
      tryCatch(map_local(st$last_x, st$graph, main = st$last_label), error = function(e) {
        graphics::plot.new()
        graphics::text(0.5, 0.5, conditionMessage(e))
      })
    }
    output$map <- shiny::renderPlot(draw_map(), res = 96)
    output$save_map <- shiny::downloadHandler(
      filename = function() "ivspatial_map.png",
      content = function(file) {
        grDevices::png(file, width = 1500, height = 560, res = 120)
        on.exit(grDevices::dev.off())
        draw_map()
      })

    shiny::observeEvent(input$n_go, {
      st$null <- tryCatch({
        if (is.null(st$graph)) st$graph <- make_graph()
        q <- suppressWarnings(as.numeric(input$n_q))
        s <- shiny::withProgress(message = "Computing the exact covariance...", value = 0.5,
                                 null_limit(st$graph, input$n_d, input$n_s, if (is.finite(q)) q else NULL))
        paste(utils::capture.output(print(s)), collapse = "\n")
      }, error = function(e) paste("Error:", conditionMessage(e)))
    })
    output$null <- shiny::renderText(st$null)
  }

  shiny::shinyApp(ui, server)
}

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
