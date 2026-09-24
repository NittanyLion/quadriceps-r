#' quadriceps: positive-weight cubature rules for the Gaussian weight and the cube
#'
#' Positive-weight cubature rules in `d >= 1` dimensions for two weight functions:
#'
#' * **GH** ([ghpos()]): the Gaussian weight, by default the standard normal density
#'   `N(0, I_d)`;
#' * **Le** ([lepos()]): the uniform weight, by default the uniform density on `[0,1]^d`.
#'
#' A rule of degree `p` integrates every polynomial of total degree `<= p` exactly. All
#' weights are strictly positive. The rules shipped with the package are the smallest ones
#' known to its author; [available_rules()] lists them, and `LICENSE.note` and the files `RULES.md` and
#' `NOTICE.md` of the repository say where each one comes from.
#'
#' The second argument `q` of [ghpos()] and [lepos()] is the number of nodes of the
#' one-dimensional Gauss rule whose exactness is wanted, as in
#' `statmod::gauss.quad(q, kind = "hermite")`; the rule returned has degree `p = 2q - 1`.
#' The degree can be given instead, as the argument `p`.
#'
#' This package is the R twin of the Julia package Quadriceps.jl and of the Python package
#' quadriceps; the three share their data and their conventions.
#'
#' @examples
#' r <- ghpos(3, 4)        # as exact as the 4x4x4 Gauss-Hermite grid (degree 7), with 27 nodes
#' dim(r$nodes)
#' r <- ghpos(3, p = 7)    # the same rule, requested by degree
#' r <- lepos(2, 5)        # degree 9 on the unit square: 17 nodes instead of 25
#' @keywords internal
"_PACKAGE"
