# ivspatial (R) — nonparametric tests of spatial independence

`ivspatial` tests whether the values of a variable observed at spatial units — cells of a
grid, polygons such as counties, irregularly located points, or the nodes of a network — are
independent and identically distributed (i.i.d.) across the units, with the procedure of

> Ashtari Nezhad, E. (2026). *A nonparametric test of spatial independence on neighbourhood graphs.*

The tests are **rank based and exact**: under independence their P-values are valid for every
sample size, every graph and every continuous distribution, with no weight matrix, bandwidth or
tuning to choose. Besides the P-value the package reports **what kind** of dependence is present
(in the level, in the scale, or only among several neighbours), maps **where** it occurs, and can
let the data choose the resolution with exact post-selection calibration.

The package can be used through a **window** (no programming, `ivspatial_app()`) or through
**R code**; both give the same results. A Python version with the same functions is on PyPI:
`pip install ivspatial` ([github.com/emaad7/ivspatial](https://github.com/emaad7/ivspatial)).

## Installation

```r
install.packages("remotes")                       # once
remotes::install_github("emaad7/ivspatial-r")
install.packages(c("deldir", "sf"))               # optional: Delaunay graphs, polygon files
install.packages("readxl")                        # optional: Excel files
```

or, from the source file: `install.packages("ivspatial_1.0.0.tar.gz", repos = NULL, type = "source")`.

This installs both the functions and the window: `library(ivspatial); ivspatial_app()` opens it.

## Quick start

```r
library(ivspatial)

# a grid: Mercer-Hall wheat yields
g <- lattice_from_indices(mercer_hall$row, mercer_hall$col, "rook")
res <- spatial_test(mercer_hall$grain, g, seed = 1)
res                              # report
res$p_value; res$reject          # P-value of the combined test C4 and the decision
res$components                   # T_V(4,3), T_V(4,2), RM, RS with P-values
res$features                     # standardized features of T_V(3,2)

# points: Baltimore house prices, 4 nearest neighbours (or delaunay_graph())
gk <- knn_graph(cbind(baltimore$X, baltimore$Y), k = 4)
spatial_test(log(baltimore$PRICE), gk, seed = 1)
residual_test(log(baltimore$PRICE), baltimore[, c("NROOM", "NBATH", "SQFT", "AGE")], gk, seed = 1)

# polygons (package sf)
nc <- read_polygons(system.file("shape/nc.shp", package = "sf"))
gq <- polygon_graph(nc, "queen")
spatial_test(nc$SID79 / nc$BIR79, gq, seed = 1)

# many variables, data-driven resolution, lags, maps, null limit
screen_variables(mercer_hall, g, columns = c("grain", "straw"), seed = 1)
spatial_test(log(baltimore$PRICE), gk, resolution = "auto", lags = 2, seed = 1)
map_local(mercer_hall$grain, g)
null_limit(g, d = 4, s = 2, value = 35)
```

## The neighbourhood graph

| Data | Function |
|---|---|
| points | `knn_graph(xy, k)`, `delaunay_graph(xy)`, `distance_band_graph(xy, radius)` |
| grid | `lattice_graph(nrow, ncol, "rook")`, `lattice_from_indices(row, col, "queen")` |
| polygons | `polygon_graph(sf_object_or_file, "queen")` |
| network / your own | `graph_from_edges(i, j)`, `graph_from_neighbors(nb)` (e.g. an spdep `nb` object) |
| spatial lag | `lag_graph(g, 2)` |

`build_graph(data, type, ...)` builds any of them from the columns of a data frame.

## Reading the output

* **C4** is the proposed test: the minimum P-value of four components, calibrated exactly from the
  same random rank arrangements. A P-value below the level means that the values are not i.i.d.
* **T_V(4,3)** sees dependence in the joint distribution of a site and two neighbours, including
  dependence that is invisible in pairs; **T_V(4,2)** pairs with resolved tails; **RM** linear
  association of ranks (level); **RS** association of squared ranks (scale).
* **Standardized features** (about N(0,1) under independence): `c.lin|n.lin` similarity in level
  between a site and its neighbours; `c.quad|n.quad` extreme values next to extreme values;
  `n.lin^2`, `n.quad^2` the same between neighbours. |value| > 3 is notable.
* **Local features** (`local_features`, `map_local`, or the *Map* tab) show where a feature
  receives its contribution; they are descriptive, not tests.

The full user guide is `vignette("ivspatial")`, or `ivspatial_guide()` (with screenshots of the window).

## Notes

* The null hypothesis is that the values are i.i.d. A smooth trend or heteroskedasticity is
  detected as dependence; model the trend and test the residuals if needed.
* Residual tests are valid without correction when the regressors are exogenous. They are not
  valid for spatial-lag models (y regressed on W y).
* Ties (discrete data) are broken at random, which keeps the tests exact.

## Citation

`citation("ivspatial")`:

Ashtari Nezhad, E. (2026). A nonparametric test of spatial independence on neighbourhood graphs. Manuscript.

Ashtari Nezhad, E. (2026). ivspatial: Nonparametric Tests of Spatial Independence. R package
version 1.0.0. <https://github.com/emaad7/ivspatial-r>

License: MIT.
