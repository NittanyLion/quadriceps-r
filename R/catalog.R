# The catalog of stored rules (extdata/index.tsv) and the reader of extdata/rules.bin.
# All rules are in the one binary file rules.bin, format QUADRICEPS1 (see FORMAT.md): an ASCII
# header line, an index of little-endian 64-bit integers, then flat little-endian doubles.

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
                                        "interior", "origin", "source_id")]
    rownames(x) <- NULL
    bin <- read_bin_index(file.path(datadir(), "rules.bin"))
    m <- match(paste(x$family, x$d, x$p), paste(bin$family, bin$d, bin$p))
    if (anyNA(m) || any(bin$n[m] != x$n)) stop("quadriceps: index.tsv and rules.bin disagree", call. = FALSE)
    .state$index <- x
    .state$bin <- bin[m, ]                             # row i of the catalog <-> row i of the bin index
    .state$cache <- new.env(parent = emptyenv())
  }
  .state$index
}

# 64-bit little-endian integers, read as pairs of 32-bit halves (R has no 64-bit integer type);
# every value in the index is far below 2^53, so doubles hold them exactly.
read_int64 <- function(con, n) {
  lo <- readBin(con, "integer", n = 2L * n, size = 4L, endian = "little")
  hi <- lo[c(FALSE, TRUE)]; lo <- lo[c(TRUE, FALSE)]
  hi * 2^32 + ifelse(lo < 0, lo + 2^32, lo)
}

# The index of rules.bin as a data frame: family, d, p, q, n, offset, nbytes, source_id.
read_bin_index <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  header <- readLines(con, n = 1L)
  tok <- strsplit(header, " ", fixed = TRUE)[[1L]]
  if (length(tok) == 0L || tok[1L] != "QUADRICEPS1") stop(path, ": not a QUADRICEPS1 file", call. = FALSE)
  kv <- strsplit(tok[-1L], "=", fixed = TRUE)
  kv <- stats::setNames(vapply(kv, `[`, "", 2L), vapply(kv, `[`, "", 1L))
  if (kv[["fmt"]] != "1" || kv[["endian"]] != "little" || kv[["float"]] != "binary64" || kv[["index_fields"]] != "8")
    stop(path, ": unsupported QUADRICEPS1 variant", call. = FALSE)
  cells <- as.integer(kv[["cells"]])
  seek(con, nchar(header, type = "bytes") + 1L)        # readLines may have buffered past the header
  v <- matrix(read_int64(con, 8L * cells), ncol = 8L, byrow = TRUE)
  bin <- data.frame(family = families[v[, 1L] + 1L], d = v[, 2L], p = v[, 3L], q = v[, 4L], n = v[, 5L],
                    offset = v[, 6L], nbytes = v[, 7L], source_id = v[, 8L], stringsAsFactors = FALSE)
  if (any(bin$nbytes != bin$n * (bin$d + 1) * 8)) stop(path, ": inconsistent index", call. = FALSE)
  bin
}

# The stored rule behind a catalog row, in the normalized frame.
stored <- function(info) {
  x <- catalog()
  i <- which(x$family == info$family & x$d == info$d & x$p == info$p)
  key <- paste(info$family, info$d, info$p, sep = "_")
  if (is.null(.state$cache[[key]])) {
    b <- .state$bin[i, ]
    con <- file(file.path(datadir(), "rules.bin"), "rb"); on.exit(close(con))
    seek(con, b$offset)
    v <- readBin(con, "double", n = b$n * (b$d + 1), size = 8L, endian = "little")
    a <- matrix(v, ncol = b$d + 1L, byrow = TRUE)      # row-major on disk
    .state$cache[[key]] <- list(nodes = a[, seq_len(b$d), drop = FALSE], weights = a[, b$d + 1L])
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
#'   * `source_id`: the origin as the small integer stored in `rules.bin` (0: own; 3: derived
#'     from Diallo and Worku; 10 and up: a published rule; see `FORMAT.md` in the repository).
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
