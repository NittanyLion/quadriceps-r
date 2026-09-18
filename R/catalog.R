# The catalog of stored rules (extdata/index.tsv) and the loader for the rule files.

.state <- new.env(parent = emptyenv())

families <- c("gh", "le")

check_family <- function(family) {
  if (!is.character(family) || length(family) != 1L || !(family %in% families))
    stop("unknown family; use \"gh\" or \"le\"", call. = FALSE)
  family
}

# Directory of the data: the installed package's extdata, or options(quadriceps.datadir = …)
# when the code is sourced from a checkout (tests/run_tests.R does that).
datadir <- function() {
  dir <- getOption("quadriceps.datadir")
  if (is.null(dir)) dir <- system.file("extdata", package = "quadriceps")
  if (!nzchar(dir) || !dir.exists(dir)) stop("quadriceps: data directory not found", call. = FALSE)
  dir
}

catalog <- function() {
  if (is.null(.state$index)) {
    x <- utils::read.delim(file.path(datadir(), "index.tsv"), comment.char = "#", quote = "",
                           stringsAsFactors = FALSE, encoding = "UTF-8")
    x$interior <- x$interior == "yes"
    x$q <- (x$p + 1L) %/% 2L
    x <- x[order(x$family, x$d, x$p), c("family", "d", "q", "p", "n", "moller", "relerr", "minweight",
                                        "interior", "origin", "file")]
    rownames(x) <- NULL
    .state$index <- x
    .state$cache <- new.env(parent = emptyenv())
  }
  .state$index
}

# The stored rule behind a catalog row, in the normalized frame.
stored <- function(info) {
  cat <- catalog()                                   # makes sure the cache exists
  key <- paste(info$family, info$d, info$p, sep = "_")
  if (is.null(.state$cache[[key]])) {
    v <- scan(file.path(datadir(), info$family, info$file), sep = ",", comment.char = "#", quiet = TRUE)
    a <- matrix(v, ncol = info$d + 1L, byrow = TRUE)
    if (nrow(a) != info$n) stop(info$file, ": ", nrow(a), " nodes, the catalog says ", info$n, call. = FALSE)
    .state$cache[[key]] <- list(nodes = a[, seq_len(info$d), drop = FALSE], weights = a[, info$d + 1L])
  }
  .state$cache[[key]]
}

#' List the stored rules
#'
#' All stored rules of a family, one row per rule, sorted by dimension and degree. Tensor
#' products, which `pragmatic = TRUE` builds on demand, are not listed.
#'
#' @param family `"gh"` (Gaussian weight) or `"le"` (uniform weight on the cube).
#' @return A data frame with columns
#'   * `family`, `d` (dimension), `q`, `p` (degree of exactness, `p = 2q - 1`), `n` (number of
#'     nodes);
#'   * `moller`: Möller's lower bound on `n`, or `-1` where it is not tabulated; a rule with
#'     `n == moller` is proven minimal;
#'   * `relerr`: largest relative monomial error of the stored double-precision rule, measured
#'     when the data were built (see [exactness_error()]);
#'   * `minweight`: smallest weight (always `> 0`);
#'   * `interior`: `TRUE` when every node lies inside the integration domain;
#'   * `origin`: who the rule belongs to: `"own"`, or `"derived: ..."`, `"same-rule: ..."`,
#'     `"transcribed: ..."` followed by the published source;
#'   * `file`: file name under `extdata/<family>/`.
#' @examples
#' head(available_rules("gh")[, c("d", "q", "p", "n")])
#' @export
available_rules <- function(family) {
  check_family(family)
  x <- catalog()
  x <- x[x$family == family, ]
  rownames(x) <- NULL
  x
}
