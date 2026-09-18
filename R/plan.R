# Which rule answers a request (d, p): a stored one, or the cheapest tensor product.
# A plan is a list of atoms; an atom is list(d, n, info), info = NULL for the 1-d Gauss rule.

as_count <- function(x, name, min) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x != round(x))
    stop(name, " must be a whole number", call. = FALSE)
  if (x < min) stop(name, " must be at least ", min, ", got ", x, call. = FALSE)
  as.integer(x)
}

# The degree a call asks for: p itself, or 2q - 1. Exactly one of q, p must be given.
degree <- function(q, p) {
  if (is.null(q) == is.null(p)) stop("give exactly one of q and p", call. = FALSE)
  if (!is.null(p)) as_count(p, "p", 0L) else 2L * as_count(q, "q", 1L) - 1L
}

# Nodes of the one-dimensional Gauss rule of degree >= p.
gaussq <- function(p) p %/% 2L + 1L

# The smallest stored rule of dimension d and degree >= p (ties: the lowest degree), or NULL.
# A rule exact to degree p' >= p is a rule of degree p.
best_stored <- function(family, d, p) {
  x <- catalog()
  x <- x[x$family == family & x$d == d & x$p >= p, ]
  if (nrow(x) == 0L) return(NULL)
  x[order(x$n, x$p)[1L], ]
}

# The single-rule answer for (d, p): Gauss in one dimension, a stored rule otherwise.
atom <- function(family, d, p) {
  if (d == 1L) return(list(d = 1L, n = gaussq(p), info = NULL))
  info <- best_stored(family, d, p)
  if (is.null(info)) NULL else list(d = d, n = info$n, info = info)
}

# The cheapest product of atoms covering d dimensions at degree p, by dynamic programming over
# the dimension (a product of products is a product, so splitting in two is enough). A single
# stored rule wins ties. Costs are doubles: q^d overflows integers early.
cheapest <- function(family, d, p) {
  cost <- numeric(d); parts <- vector("list", d)
  for (k in seq_len(d)) {
    a <- atom(family, k, p)
    if (!is.null(a)) { cost[k] <- a$n; parts[[k]] <- list(a) }
    for (j in seq_len(k %/% 2L)) {
      c <- cost[j] * cost[k - j]
      if (is.null(parts[[k]]) || c < cost[k]) { cost[k] <- c; parts[[k]] <- c(parts[[j]], parts[[k - j]]) }
    }
  }
  list(parts = parts[[d]], n = cost[d])
}

norule <- function(family, d, p) {
  x <- catalog(); ps <- x$p[x$family == family & x$d == d]
  have <- if (length(ps) == 0L) paste0("no rules are stored for d = ", d)
          else paste0("stored rules for d = ", d, " reach q = ", (max(ps) + 1L) %/% 2L, ", p = ", max(ps))
  at <- if (p %% 2L == 1L) paste0("q = ", (p + 1L) %/% 2L, " (p = ", p, ")") else paste0("p = ", p)
  cond <- structure(class = c("quadriceps_no_rule", "error", "condition"),
                    list(message = paste0("no stored positive-weight ", if (family == "gh") "GH" else "Le",
                                          " rule for d = ", d, ", ", at, ": ", have,
                                          "; pragmatic = TRUE returns the cheapest tensor product of ",
                                          "lower-dimensional rules instead"), call = NULL))
  stop(cond)
}

plan <- function(family, d, q, p, pragmatic) {
  check_family(family)
  d <- as_count(d, "d", 1L)
  p <- degree(q, p)
  if (isTRUE(pragmatic)) return(cheapest(family, d, p))
  a <- atom(family, d, p)
  if (is.null(a)) norule(family, d, p)
  list(parts = list(a), n = a$n)
}

# The one-dimensional Gauss rule in the normalized frame (N(0,1) for GH, uniform on [0,1] for
# Le), by the Golub-Welsch eigenvalue method: the nodes are the eigenvalues of the Jacobi
# matrix of the orthonormal polynomials, the weights the squared first components of its
# eigenvectors (both weights have total mass 1).
gauss1d <- function(family, q) {
  if (q == 1L) return(list(nodes = matrix(if (family == "gh") 0 else 0.5, 1L, 1L), weights = 1))
  k <- seq_len(q - 1L)
  off <- if (family == "gh") sqrt(k) else k / sqrt(4 * k^2 - 1)
  J <- matrix(0, q, q)
  J[cbind(k, k + 1L)] <- off; J[cbind(k + 1L, k)] <- off
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  x <- e$values[o]; w <- e$vectors[1L, o]^2
  x <- (x - rev(x)) / 2; w <- (w + rev(w)) / 2         # the rule is symmetric; remove rounding asymmetry
  w <- w / sum(w)
  if (family == "le") x <- (x + 1) / 2
  list(nodes = matrix(x, ncol = 1L), weights = w)
}

# Tensor product of rules; the first factor varies slowest.
tensor <- function(rules) {
  n <- prod(vapply(rules, function(r) length(r$weights), 0))
  d <- sum(vapply(rules, function(r) ncol(r$nodes), 0L))
  X <- matrix(0, n, d); w <- rep(1, n)
  rep_ <- n; col <- 0L; idx <- seq_len(n) - 1
  for (r in rules) {
    nk <- length(r$weights); rep_ <- rep_ / nk
    j <- (idx %/% rep_) %% nk + 1
    w <- w * r$weights[j]
    X[, col + seq_len(ncol(r$nodes))] <- r$nodes[j, , drop = FALSE]
    col <- col + ncol(r$nodes)
  }
  list(nodes = X, weights = w)
}

materialize <- function(family, pl) {
  d <- sum(vapply(pl$parts, function(a) a$d, 0L))
  if (pl$n * (d + 1) > 2^31 - 1)
    stop("the cheapest rule for this request has ", format(pl$n, scientific = FALSE, big.mark = ","),
         " nodes, which is more than this function will build", call. = FALSE)
  rules <- lapply(pl$parts, function(a) if (is.null(a$info)) gauss1d(family, a$n) else stored(a$info))
  if (length(rules) == 1L) rules[[1L]] else tensor(rules)
}
