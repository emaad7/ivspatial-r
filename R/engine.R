# Computational engine: quantile classes, the features of T_V(d, s), their exact null
# covariance, rank Moran statistics and weighted chi-square limits.
#
# For a centre v with neighbours N(v) and a feature tau = (c, mu), where c is the label of
# the basis function carried by the centre (0: not involved) and mu[a] is the number of
# neighbours carrying e_a,
#   f_{v,tau} = e_c(delta_v) |O_v(tau)|^{-1/2} sum_{A in O_v(tau)} prod_{(u,a) in A} e_a(delta_u),
#   W_tau = K_tau^{-1/2} sum_v f_{v,tau},   T_V(d, s) = sum_tau W_tau^2.
# The neighbour sums are monomial symmetric functions, computed from the power sums
# P_S = A %*% prod_{l in S} X_l for many rank arrangements at once.

# E[class + 1, a + 1]: orthonormal polynomials on {0, ..., d-1} (uniform weights), e_0 = 1,
# e_a of degree a, normalized so that e_a(highest class) > 0
basis_matrix <- function(d) {
  V <- outer(0:(d - 1), 0:(d - 1), "^")
  Q <- qr.Q(qr(V))
  E <- Q * sqrt(d)
  E[, 1] <- 1
  for (a in seq_len(d - 1) + 1L) if (E[d, a] < 0) E[, a] <- -E[, a]
  E
}

# multisets of `size` labels from 1..(d-1), as count vectors, in the order of
# itertools.combinations_with_replacement
multisets <- function(d, size) {
  if (size == 0) return(list(integer(d - 1)))
  out <- list()
  rec <- function(start, left, acc) {
    if (left == 0) {
      mu <- integer(d - 1)
      for (a in acc) mu[a] <- mu[a] + 1L
      out[[length(out) + 1L]] <<- mu
      return(invisible())
    }
    for (a in start:(d - 1)) rec(a, left - 1, c(acc, a))
  }
  rec(1L, size, integer())
  out
}

feature_list <- function(d, s_max) {
  feats <- list()
  for (c in 0:(d - 1)) for (s_nb in 0:s_max) {
    s <- (c > 0) + s_nb
    if (s >= 2 && s <= s_max) for (mu in multisets(d, s_nb)) feats[[length(feats) + 1L]] <- list(c = c, mu = mu)
  }
  feats
}

LAB <- c("lin", "quad", "cub", "quart", "quint")
lab_name <- function(a) if (a <= length(LAB)) LAB[a] else paste0("e", a)

feature_name <- function(f) {
  parts <- character()
  for (a in seq_along(f$mu)) if (f$mu[a] > 0)
    parts <- c(parts, paste0(lab_name(a), if (f$mu[a] > 1) paste0("^", f$mu[a]) else ""))
  nbp <- if (length(parts)) paste0("n.", paste(parts, collapse = "*")) else ""
  if (f$c > 0) {
    if (nzchar(nbp)) paste0("c.", lab_name(f$c), "|", nbp) else paste0("c.", lab_name(f$c))
  } else nbp
}

orbit_size <- function(k, mu) {
  s <- sum(mu)
  ifelse(k < s, 0, exp(lfactorial(k) - sum(lfactorial(mu)) - lfactorial(pmax(k - s, 0))))
}

# 0-based ranks with ties broken at random (keeps the tests exact for discrete data)
rank0 <- function(x) {
  o <- order(x, stats::runif(length(x)))
  r <- integer(length(x))
  r[o] <- seq_along(x) - 1L
  r
}

tv_engine <- function(graph, d, s_max) {
  if (d < 2 || !(s_max %in% c(2, 3))) stop("d must be at least 2 and s 2 or 3", call. = FALSE)
  feats <- feature_list(d, s_max)
  deg <- graph$degree
  size <- t(vapply(feats, function(f) round(orbit_size(deg, f$mu)), numeric(graph$n)))
  if (graph$n == 1) size <- matrix(size, ncol = 1)
  keep <- rowSums(size > 0) > 0
  feats <- feats[keep]
  size <- size[keep, , drop = FALSE]
  Kt <- rowSums(size > 0)
  scale <- ifelse(size > 0, 1 / sqrt(pmax(size, 1)), 0)
  list(graph = graph, d = d, s_max = s_max, n = graph$n, A = adjacency(graph), E = basis_matrix(d),
       feats = feats, size = size, Kt = Kt, scale = scale,
       names = vapply(feats, feature_name, character(1)))
}

# local contributions: array (features, n, B) for a class matrix delta (n x B, values 0..d-1);
# with sum = TRUE the (features x B) matrix of their sums over the sites
tv_local <- function(eng, delta, sum = FALSE) {
  delta <- as.matrix(delta)
  n <- nrow(delta)
  B <- ncol(delta)
  E <- eng$E
  X <- lapply(seq_len(eng$d - 1), function(a) matrix(E[delta + 1L, a + 1L], n, B))
  cache <- new.env(hash = TRUE)
  P <- function(ls) {
    ls <- sort(ls)
    key <- paste(ls, collapse = ".")
    v <- cache[[key]]
    if (is.null(v)) {
      prod <- X[[ls[1]]]
      for (l in ls[-1]) prod <- prod * X[[l]]
      v <- as.matrix(eng$A %*% prod)
      assign(key, v, envir = cache)
    }
    v
  }
  msym <- function(mu) {
    labels <- rep(seq_along(mu), mu)
    s <- length(labels)
    if (s == 0) return(NULL)
    if (s == 1) return(P(labels))
    if (s == 2) {
      a <- labels[1]
      b <- labels[2]
      if (a == b) return(0.5 * (P(a)^2 - P(c(a, a))))
      return(P(a) * P(b) - P(c(a, b)))
    }
    a <- labels[1]
    b <- labels[2]
    cc <- labels[3]
    if (a == b && b == cc) return((P(a)^3 - 3 * P(a) * P(c(a, a)) + 2 * P(c(a, a, a))) / 6)
    if (a == b || b == cc) {
      if (b == cc) {
        tmp <- a
        a <- cc
        cc <- tmp
      }
      return(0.5 * (P(a)^2 * P(cc) - P(c(a, a)) * P(cc) - 2 * P(c(a, cc)) * P(a) + 2 * P(c(a, a, cc))))
    }
    P(a) * P(b) * P(cc) - P(c(a, b)) * P(cc) - P(c(a, cc)) * P(b) - P(c(b, cc)) * P(a) + 2 * P(c(a, b, cc))
  }
  out <- if (sum) matrix(0, length(eng$feats), B) else array(0, c(length(eng$feats), n, B))
  for (i in seq_along(eng$feats)) {
    f <- eng$feats[[i]]
    m <- msym(f$mu)
    val <- if (is.null(m)) matrix(E[delta + 1L, f$c + 1L], n, B) else if (f$c == 0) m else
      matrix(E[delta + 1L, f$c + 1L], n, B) * m
    if (sum) out[i, ] <- colSums(val * eng$scale[i, ]) else out[i, , ] <- val * eng$scale[i, ]
  }
  out
}

# feature vectors: (features x B) matrix
tv_W <- function(eng, delta) tv_local(eng, delta, sum = TRUE) / sqrt(eng$Kt)

tv_T <- function(eng, delta) colSums(tv_W(eng, delta)^2)

# exact null covariance of the features (hashing of the labelled sets)
tv_sigma <- function(eng) {
  d <- eng$d
  smax <- eng$s_max
  n <- eng$n
  nbs <- eng$graph$neighbors
  p <- length(eng$feats)
  fkey <- vapply(eng$feats, function(f) paste(f$c, paste(f$mu, collapse = ",")), character(1))
  base_s <- n + 1
  big <- (base_s^smax) * (d^smax) < 2^52
  keys <- list()
  fids <- list()
  wts <- list()
  cnt <- 0L
  for (j in 1:smax) {
    # all j-subsets of neighbours of every centre: matrix U (rows) and centre vector V
    sub <- lapply(seq_len(n), function(v) {
      a <- nbs[[v]]
      if (length(a) < j) return(NULL)
      m <- if (length(a) == 1 || j == 1) matrix(a, ncol = 1) else t(utils::combn(a, j))
      if (j == 1) m <- matrix(a, ncol = 1)
      cbind(v, m)
    })
    sub <- do.call(rbind, sub)
    if (is.null(sub) || !nrow(sub)) next
    V <- sub[, 1]
    U <- sub[, -1, drop = FALSE]
    labs <- as.matrix(expand.grid(rep(list(1:(d - 1)), j)))
    labs <- labs[, rev(seq_len(j)), drop = FALSE]
    for (li in seq_len(nrow(labs))) {
      lab <- labs[li, ]
      mu <- tabulate(lab, d - 1)
      for (c in 0:(d - 1)) {
        fi <- match(paste(c, paste(mu, collapse = ",")), fkey)
        if (is.na(fi)) next
        m <- nrow(U)
        if (c == 0) {
          S <- U
          Lm <- matrix(lab, m, j, byrow = TRUE)
        } else {
          pos <- rowSums(U < V)
          S <- matrix(0L, m, j + 1)
          Lm <- matrix(0L, m, j + 1)
          for (k in 1:(j + 1)) {
            before <- k <= pos
            at <- k == pos + 1
            after <- k > pos + 1
            if (any(before)) {
              S[before, k] <- U[before, k]
              Lm[before, k] <- lab[k]
            }
            S[at, k] <- V[at]
            Lm[at, k] <- c
            if (any(after)) {
              S[after, k] <- U[after, k - 1]
              Lm[after, k] <- lab[k - 1]
            }
          }
        }
        t <- ncol(S)
        if (t < smax) {
          S <- cbind(S, matrix(0L, m, smax - t))
          Lm <- cbind(Lm, matrix(0L, m, smax - t))
        }
        if (big) {
          ks <- 0
          kl <- 0
          for (k in 1:smax) {
            ks <- ks * base_s + S[, k]
            kl <- kl * d + Lm[, k]
          }
          key <- ks * d^smax + kl
        } else {
          key <- do.call(paste, c(as.data.frame(cbind(S, Lm)), sep = "."))
        }
        cnt <- cnt + 1L
        keys[[cnt]] <- key
        fids[[cnt]] <- rep(fi, m)
        wts[[cnt]] <- eng$scale[fi, V] / sqrt(eng$Kt[fi])
      }
    }
  }
  key <- unlist(keys)
  g <- match(key, unique(key))
  C <- Matrix::sparseMatrix(i = g, j = unlist(fids), x = unlist(wts), dims = c(max(g), p))
  S <- as.matrix(Matrix::crossprod(C))
  dimnames(S) <- list(eng$names, eng$names)
  S
}

# number of labelled sets enumerated by tv_sigma (to guard against very large designs)
n_sets <- function(graph, d, s) {
  k <- graph$degree
  sum(vapply(0:s, function(j) sum(choose(k, j)) * (d - 1)^j * d, numeric(1)))
}

# ----------------------------------------------------------------------
# Rank Moran statistics RM = |I(ranks) + 1/(n-1)|, RS = |I(squared standardized ranks) + 1/(n-1)|
# ----------------------------------------------------------------------
rank_moran <- function(graph) {
  deg <- pmax(graph$degree, 1)
  W <- Matrix::Diagonal(x = 1 / deg) %*% adjacency(graph)
  list(W = W, S0 = sum(W), n = graph$n)
}

moran_I <- function(rm, Y) {
  Yc <- sweep(Y, 2, colMeans(Y))
  num <- colSums(Yc * as.matrix(rm$W %*% Yc))
  den <- colSums(Yc * Yc)
  rm$n / rm$S0 * num / den
}

rank_moran_stats <- function(rm, R) {
  R <- as.matrix(R) * 1
  n <- rm$n
  sdv <- sqrt(colMeans(sweep(R, 2, colMeans(R))^2))
  Z <- sweep(R - (n - 1) / 2, 2, sdv, "/")
  e <- 1 / (n - 1)
  list(RM = abs(moran_I(rm, R) + e), RS = abs(moran_I(rm, Z * Z) + e))
}

# ----------------------------------------------------------------------
# Weighted chi-square limits: P(sum lam_i chi2_1 > t) (Imhof, 1961)
# ----------------------------------------------------------------------
group_weights <- function(lam, tol = 1e-8) {
  lam <- sort(lam[lam > tol], decreasing = TRUE)
  w <- numeric()
  df <- numeric()
  for (x in lam) {
    k <- length(w)
    if (k && abs(x - w[k]) <= 1e-6 * max(1, x)) df[k] <- df[k] + 1 else {
      w <- c(w, x)
      df <- c(df, 1)
    }
  }
  list(w = w, df = df)
}

chernoff <- function(t, w, df) {
  th <- 0.25 / max(w)
  min(1, exp(-th * t - 0.5 * sum(df * log(1 - 2 * th * w))))
}

weighted_chi2_sf <- function(t, lam) {
  g <- group_weights(lam)
  w <- g$w
  df <- g$df
  if (!length(w)) return(if (t > 0) 0 else 1)
  if (length(w) == 1) return(stats::pchisq(t / w, df, lower.tail = FALSE))
  mean <- sum(w * df)
  sdv <- sqrt(2 * sum(w^2 * df))
  if (t > mean + 12 * sdv) return(chernoff(t, w, df))
  f <- function(u) vapply(u, function(uu) {
    if (uu == 0) return(0.5 * (sum(df * w) - t))
    theta <- 0.5 * sum(df * atan(w * uu)) - 0.5 * t * uu
    rho <- exp(sum(df / 4 * log1p((w * uu)^2)))
    sin(theta) / (uu * rho)
  }, numeric(1))
  r <- tryCatch(stats::integrate(f, 0, Inf, subdivisions = 1000L, rel.tol = 1e-10, abs.tol = 1e-12,
                                 stop.on.error = FALSE),
                error = function(e) list(value = NA, abs.error = Inf))
  p <- 0.5 + r$value / pi
  if (!is.finite(p) || p < 1e-6 || r$abs.error > 1e-4) {
    bound <- chernoff(t, w, df)
    p <- if (is.finite(p) && p > 0) min(bound, p) else min(bound, 1e-6)
  }
  min(max(p, 0), 1)
}

weighted_chi2_quantile <- function(alpha, lam) {
  g <- group_weights(lam)
  mean <- sum(g$w * g$df)
  sdv <- sqrt(2 * sum(g$w^2 * g$df))
  stats::uniroot(function(t) weighted_chi2_sf(t, lam) - alpha, c(1e-9, mean + 30 * sdv), tol = 1e-9)$root
}
