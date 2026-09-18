#' Positive-weight cubature rule for the Gaussian weight
#'
#' The smallest positive-weight rule the package has for the Gaussian weight in `d`
#' dimensions.
#'
#' `ghpos(d, q)` follows one-dimensional Gauss-Hermite functions such as
#' `statmod::gauss.quad(q, kind = "hermite")`: there, `q` is the number of nodes of the
#' one-dimensional Gauss rule, which is exact to degree `2q - 1`. Here, `ghpos(d, q)` returns
#' a rule with the same exactness in `d` dimensions, degree `p = 2q - 1`, that replaces the
#' `q^d`-node product grid; for `d = 1` it is the `q`-node Gauss-Hermite rule itself.
#' Alternatively give the degree: `ghpos(d, q)` is `ghpos(d, p = 2 * q - 1)`.
#'
#' With `r <- ghpos(d, q)`, `sum(r$weights * f(r$nodes))` (for a vectorized `f` of the rows)
#' equals the integral of `f` against the weight for every polynomial `f` of total degree
#' `<= p`, up to rounding; [ruleinfo()] gives the measured error of each stored rule.
#'
#' Rules are stored for `2 <= d <= 5` at odd degrees, up to a ceiling that depends on `d`
#' (see [available_rules()]). A request is served by the smallest stored rule of degree
#' `>= p`, so an even `p` gets the rule for `p + 1`. For `d = 1` the Gauss rule is computed on
#' the fly (Golub-Welsch) and is never an error.
#'
#' @param d Dimension, `d >= 1`.
#' @param q Number of nodes of the one-dimensional Gauss rule whose exactness is wanted,
#'   `q >= 1`. Give `q` or `p`, not both.
#' @param p Degree of exactness, `p >= 0`, as an alternative to `q`.
#' @param normalize With `TRUE` (the default), the weight is the standard normal density
#'   `(2 pi)^(-d/2) exp(-|x|^2 / 2)`: the weights sum to 1 and the rule computes `E f(Z)` for
#'   `Z ~ N(0, I_d)`. With `FALSE`, the weight is `exp(-|x|^2)`, the convention of
#'   `statmod::gauss.quad(kind = "hermite")`, and the weights sum to `pi^(d/2)`.
#' @param pragmatic What to do when no stored rule covers the request. With `FALSE` (the
#'   default), signal an error of class `quadriceps_no_rule`. With `TRUE`, return the cheapest
#'   (fewest nodes) tensor product of lower-dimensional rules instead: stored rules and
#'   one-dimensional Gauss-Hermite rules, combined over the split of `d` that minimizes the
#'   number of nodes. Such a product is a valid positive-weight rule of the requested degree;
#'   it is just not small. With `TRUE` the product is also returned in the rare case that it
#'   has strictly fewer nodes than the stored rule, so the result is always the cheapest the
#'   package can build.
#' @return A list with components `nodes`, an `n x d` matrix with one node per row (also for
#'   `d = 1`), and `weights`, a vector of `n` strictly positive weights.
#' @seealso [lepos()], [nnodes()], [ruleinfo()], [available_rules()]
#' @examples
#' r <- ghpos(3, 4)                                   # degree 7: 27 nodes instead of 4^3 = 64
#' sum(r$weights * r$nodes[, 1]^2 * r$nodes[, 2]^4)   # E[Z1^2 Z2^4] = 3
#' identical(ghpos(3, p = 7), r)
#'
#' try(ghpos(3, 21))                                  # nothing stored at that degree
#' r <- ghpos(3, 21, pragmatic = TRUE)                # tensor product of lower-dimensional rules
#' r <- ghpos(7, 5, pragmatic = TRUE)                 # d = 7 as (d = 2) x (d = 5)
#' @export
ghpos <- function(d, q = NULL, p = NULL, normalize = TRUE, pragmatic = FALSE) {
  pl <- plan("gh", d, q, p, pragmatic)
  r <- materialize("gh", pl)
  if (!isTRUE(normalize)) {               # int f(x) exp(-|x|^2) dx = pi^(d/2) E f(Z / sqrt(2))
    r$nodes <- r$nodes / sqrt(2)
    r$weights <- r$weights * pi^(ncol(r$nodes) / 2)
  }
  r
}

#' Positive-weight cubature rule for the uniform weight on the cube
#'
#' The smallest positive-weight rule the package has for the uniform weight on a
#' `d`-dimensional cube.
#'
#' `lepos(d, q)` follows one-dimensional Gauss-Legendre functions such as
#' `statmod::gauss.quad(q, kind = "legendre")`: `q` is the number of nodes of the
#' one-dimensional Gauss rule, and `lepos(d, q)` returns a rule with the same exactness in `d`
#' dimensions, degree `p = 2q - 1`, that replaces the `q^d`-node product grid; for `d = 1` it
#' is the `q`-node Gauss-Legendre rule itself. Alternatively give the degree `p`.
#'
#' The details given for [ghpos()] apply here too. All nodes of every stored Le rule lie
#' inside the cube.
#'
#' @inheritParams ghpos
#' @param normalize With `TRUE` (the default), the weight is the uniform density on `[0,1]^d`:
#'   the weights sum to 1 and the rule computes `E f(U)` for `U` uniform on the unit cube.
#'   With `FALSE`, the rule is for the plain integral over `[-1,1]^d`, the convention of
#'   `statmod::gauss.quad(kind = "legendre")`, and the weights sum to `2^d`.
#' @param pragmatic As for [ghpos()]: `FALSE` signals an error of class `quadriceps_no_rule`
#'   when no stored rule covers the request; `TRUE` returns the cheapest tensor product of
#'   lower-dimensional rules (stored rules and one-dimensional Gauss-Legendre rules).
#' @return A list with components `nodes` (`n x d` matrix, one node per row) and `weights`
#'   (`n` strictly positive weights).
#' @seealso [ghpos()], [nnodes()], [ruleinfo()], [available_rules()]
#' @examples
#' r <- lepos(2, 5)                                   # degree 9: 17 nodes on [0,1]^2 instead of 25
#' sum(r$weights * r$nodes[, 1]^3 * r$nodes[, 2]^2)   # 1/4 * 1/3
#'
#' r <- lepos(2, 5, normalize = FALSE)                # the same rule on [-1,1]^2, weights sum to 4
#' r <- lepos(6, 6, pragmatic = TRUE)                 # d = 6 from a product of stored rules
#' @export
lepos <- function(d, q = NULL, p = NULL, normalize = TRUE, pragmatic = FALSE) {
  pl <- plan("le", d, q, p, pragmatic)
  r <- materialize("le", pl)
  if (!isTRUE(normalize)) {               # [0,1]^d -> [-1,1]^d
    r$nodes <- 2 * r$nodes - 1
    r$weights <- r$weights * 2^ncol(r$nodes)
  }
  r
}

#' Number of nodes of a rule, without building it
#'
#' The number of nodes of the rule that [ghpos()] (`family = "gh"`) or [lepos()]
#' (`family = "le"`) returns for the same arguments. Signals the same error when there is no
#' rule.
#'
#' @param family `"gh"` or `"le"`.
#' @inheritParams ghpos
#' @return The node count, as a double (tensor products in high dimensions exceed the
#'   integer range).
#' @examples
#' nnodes("gh", 5, 7)
#' nnodes("gh", 7, 5, pragmatic = TRUE)
#' 5^7                                    # the product grid
#' @export
nnodes <- function(family, d, q = NULL, p = NULL, pragmatic = FALSE) plan(family, d, q, p, pragmatic)$n

#' Describe a rule: its factors and their origins
#'
#' Describes the rule that [ghpos()] or [lepos()] returns for the same arguments.
#'
#' Use it to find out whom to cite: `origin` names the published source of every rule that is
#' not the package author's own.
#'
#' @inheritParams nnodes
#' @return A data frame with one row per tensor factor, in the order of the columns of
#'   `nodes`, and the columns of [available_rules()]. A one-dimensional Gauss factor has
#'   `origin = "Gauss"` and `NA` in the columns that do not apply. A request answered by a
#'   single stored rule gives one row.
#' @examples
#' ruleinfo("le", 3, p = 41)$origin
#' ruleinfo("gh", 7, 5, pragmatic = TRUE)[, c("d", "p", "n", "origin")]
#' @export
ruleinfo <- function(family, d, q = NULL, p = NULL, pragmatic = FALSE) {
  pl <- plan(family, d, q, p, pragmatic)
  rows <- lapply(pl$parts, function(a) {
    if (!is.null(a$info)) return(a$info)
    data.frame(family = family, d = 1L, q = a$n, p = 2L * a$n - 1L, n = a$n, moller = NA_integer_,
               relerr = NA_real_, minweight = NA_real_, interior = NA, origin = "Gauss", file = NA_character_,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
