# One-dimensional moments of the normalized weights: E Z^e for N(0,1), E U^e for U[0,1].
moments1d <- function(family, p) {
  if (family == "le") return(1 / (seq_len(p + 1L)))
  m <- numeric(p + 1L); m[1L] <- 1
  for (e in seq(2L, by = 2L, length.out = p %/% 2L)) m[e + 1L] <- m[e - 1L] * (e - 1)     # (e-1)!!
  m
}

#' Largest relative monomial error of a rule
#'
#' The maximum, over all monomials `x^a` of total degree `<= p`, of
#' `|sum_i w_i x_i^a - E x^a| / max(sum_i |w_i| |x_i^a|, 1)`, for the normalized weight of
#' `family` (`"gh"`: `N(0, I_d)`; `"le"`: uniform on `[0,1]^d`). The denominator is the scale
#' of the sum being computed, so a value near machine epsilon means the rule is exact to
#' rounding. Rules returned with `normalize = FALSE` must be checked in the normalized frame.
#'
#' @param nodes An `n x d` matrix of nodes, one per row (or a rule: a list with components
#'   `nodes` and `weights`, in which case `weights` is taken from it).
#' @param weights The `n` weights.
#' @param p Degree up to which exactness is checked.
#' @param family `"gh"` or `"le"`.
#' @return The error, a single number.
#' @examples
#' r <- ghpos(4, 5)                       # degree 2 * 5 - 1 = 9
#' exactness_error(r, p = 9, family = "gh")
#' @export
exactness_error <- function(nodes, weights = NULL, p, family) {
  check_family(family)
  if (is.list(nodes)) { weights <- nodes$weights; nodes <- nodes$nodes }
  X <- as.matrix(nodes); w <- as.numeric(weights)
  n <- nrow(X); d <- ncol(X)
  if (length(w) != n) stop("nodes has ", n, " rows, weights has length ", length(w), call. = FALSE)
  p <- as_count(p, "p", 0L)
  m <- moments1d(family, p)
  P <- lapply(seq_len(d), function(k) outer(X[, k], 0:p, `^`))        # P[[k]][i, e+1] = x_ik^e
  A <- abs(P[[d]])
  if (d == 1L) {
    s <- drop(crossprod(P[[1L]], w)); sa <- drop(crossprod(A, abs(w)))
    return(max(abs(s - m) / pmax(sa, 1)))
  }
  low <- outer(0:p, 0:p, `+`)                                         # exponent sums of the last two coordinates
  descend <- function(k, prev, mom, r) {
    # prev: w * x_1^a_1 ... x_(k-1)^a_(k-1); mom: its exact moment; r: degree left for x_k, ..., x_d
    i <- seq_len(r + 1L)
    if (k == d - 1L) {                    # last two coordinates: all their exponent pairs in one product
      B <- prev * P[[k]][, i, drop = FALSE]
      S <- crossprod(B, P[[d]][, i, drop = FALSE])
      SA <- crossprod(abs(B), A[, i, drop = FALSE])
      E <- abs(S - mom * outer(m[i], m[i])) / pmax(SA, 1)
      return(max(E[low[i, i, drop = FALSE] <= r]))
    }
    max(vapply(0:r, function(a) descend(k + 1L, prev * P[[k]][, a + 1L], mom * m[a + 1L], r - a), 0))
  }
  descend(1L, w, 1, p)
}
