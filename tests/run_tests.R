# Tests, without a test framework: every check is an expect(); the script stops at the end if
# any failed. Runs against the installed package (R CMD check), or, from a checkout
# (`Rscript tests/run_tests.R` in the package root), against the sources in R/.
if (file.exists("R/api.R") && dir.exists("inst/extdata")) {
  options(quadriceps.datadir = "inst/extdata")
  for (f in list.files("R", full.names = TRUE)) source(f)
} else {
  library(quadriceps)
}

GATE <- 1e-11      # every stored rule passed this gate when the data were built
failures <- character(); checks <- 0L
expect <- function(ok, what) {
  checks <<- checks + 1L
  if (!isTRUE(ok)) failures <<- c(failures, what)
}
expect_error <- function(expr, what, class = "error") {
  r <- tryCatch({ force(expr); "no error" }, condition = function(e) class(e))
  expect(class %in% r, what)
}
near <- function(a, b, tol = 1e-12) abs(a - b) <= tol * max(1, abs(b))
F <- list(gh = ghpos, le = lepos)

# --- catalog
for (fam in names(F)) {
  x <- available_rules(fam)
  expect(nrow(x) > 0 && all(x$p %% 2 == 1) && all(x$d >= 2) && all(x$q == (x$p + 1) / 2), paste("catalog", fam))
}
expect_error(available_rules("la"), "unknown family")

# --- every stored rule
for (fam in names(F)) {
  x <- available_rules(fam)
  for (i in seq_len(nrow(x))) {
    r <- x[i, ]; id <- sprintf("%s d=%d p=%d", fam, r$d, r$p)
    u <- F[[fam]](r$d, p = r$p)
    expect(identical(u, F[[fam]](r$d, r$q)), paste(id, "q form"))
    expect(identical(dim(u$nodes), c(r$n, r$d)) && length(u$weights) == r$n, paste(id, "size"))
    expect(all(u$weights > 0), paste(id, "positive weights"))
    expect(near(sum(u$weights), 1), paste(id, "weights sum to 1"))
    expect(exactness_error(u, p = r$p, family = fam) < GATE, paste(id, "exactness"))
    if (fam == "le") expect(all(u$nodes > 0 & u$nodes < 1), paste(id, "interior"))
    expect(nnodes(fam, r$d, p = r$p) <= r$n, paste(id, "monotone"))
    expect(nnodes(fam, r$d, p = r$p, pragmatic = TRUE) == nnodes(fam, r$d, p = r$p), paste(id, "no product beats it"))
  }
}

# --- one dimension (reference values: the 3-node Gauss rules in closed form)
u <- ghpos(1, 3)
expect(identical(dim(u$nodes), c(3L, 1L)) && all(near(u$nodes[, 1], c(-sqrt(3), 0, sqrt(3)))) &&
       all(near(u$weights, c(1, 4, 1) / 6)), "GH 3-node rule")
u <- ghpos(1, p = 4, normalize = FALSE)                       # even degree: next odd
expect(all(near(u$nodes[, 1], c(-sqrt(1.5), 0, sqrt(1.5)))) && all(near(u$weights, sqrt(pi) * c(1, 4, 1) / 6)),
       "GH 3-node rule, exp(-x^2)")
u <- lepos(1, 3, normalize = FALSE)
expect(all(near(u$nodes[, 1], c(-sqrt(0.6), 0, sqrt(0.6)))) && all(near(u$weights, c(5, 8, 5) / 9)), "Le 3-node rule")
u <- lepos(1, 30)
expect(near(sum(u$weights), 1) && all(u$nodes > 0 & u$nodes < 1) && exactness_error(u, p = 59, family = "le") < 1e-13,
       "Le 30-node rule")
expect(exactness_error(ghpos(1, 25), p = 49, family = "gh") < 1e-12, "GH 25-node rule")

# --- normalize = FALSE
u <- ghpos(2, 4, normalize = FALSE)
expect(near(sum(u$weights), pi), "GH mass pi")
expect(near(sum(u$weights * u$nodes[, 1]^2 * u$nodes[, 2]^4), (sqrt(pi) / 2) * (3 * sqrt(pi) / 4)), "GH x1^2 x2^4")
u <- lepos(3, 3, normalize = FALSE)
expect(near(sum(u$weights), 8) && all(u$nodes > -1 & u$nodes < 1), "Le mass 8")
expect(near(sum(u$weights * u$nodes[, 1]^2 * u$nodes[, 3]^2), 8 / 9), "Le x1^2 x3^2")
expect(abs(sum(u$weights * u$nodes[, 1] * u$nodes[, 2]^2)) < 1e-14, "Le odd moment")

# --- q and p
expect(identical(ghpos(3, p = 6), ghpos(3, p = 7)) && identical(ghpos(3, p = 7), ghpos(3, 4)), "even p")
expect(identical(lepos(2, p = 0), lepos(2, 1)), "p = 0")
expect(nnodes("le", 2, p = 10) == nnodes("le", 2, 6), "nnodes q and p")
expect_error(ghpos(2), "neither q nor p")
expect_error(ghpos(2, 3, p = 5), "both q and p")

# --- no rule: error unless pragmatic
for (fam in names(F)) {
  f <- F[[fam]]; x <- available_rules(fam); q <- max(x$q[x$d == 3]) + 1
  expect_error(f(3, q), paste(fam, "beyond the degrees"), "quadriceps_no_rule")
  expect_error(f(3, p = 2 * q - 2), paste(fam, "beyond the degrees, even p"), "quadriceps_no_rule")
  expect_error(nnodes(fam, 3, q), paste(fam, "nnodes beyond the degrees"), "quadriceps_no_rule")
  expect_error(f(6, 3), paste(fam, "six dimensions"), "quadriceps_no_rule")
  expect(nnodes(fam, 3, q, pragmatic = TRUE) <= q^3, paste(fam, "fallback no dearer than the grid"))
  expect_error(f(0, 3), "d = 0"); expect_error(f(2, 0), "q = 0"); expect_error(f(2, p = -1), "p = -1")
  expect_error(f(2.5, 3), "fractional d")
}
msg <- tryCatch(ghpos(6, 3), error = conditionMessage)
expect(grepl("pragmatic = TRUE", msg, fixed = TRUE), "the error points to pragmatic")

# --- pragmatic fallback
expect(identical(ghpos(4, 5, pragmatic = TRUE), ghpos(4, 5)), "stored GH cell unchanged")
expect(identical(lepos(2, 11, pragmatic = TRUE), lepos(2, 11)), "stored Le cell unchanged")
# d = 4: the split 2 + 2 into stored d = 2 rules beats the grid (at d = 3 the fallback is the bare grid,
# since the stored d = 2 rules end at the same degree as the d = 3 ones)
x <- available_rules("gh"); q <- max(x$q[x$d == 4]) + 1
u <- ghpos(4, q, pragmatic = TRUE)
expect(nrow(u$nodes) == nnodes("gh", 4, q, pragmatic = TRUE) && nrow(u$nodes) < q^4, "fallback beyond the degrees: size")
expect(all(u$weights > 0) && near(sum(u$weights), 1) && exactness_error(u, p = 2 * q - 1, family = "gh") < GATE,
       "fallback beyond the degrees: valid")
for (fam in names(F)) for (dq in list(c(6, 4), c(7, 3), c(8, 2))) {
  d <- dq[1]; q <- dq[2]; id <- sprintf("fallback %s d=%d q=%d", fam, d, q)
  u <- F[[fam]](d, q, pragmatic = TRUE); n <- nnodes(fam, d, q, pragmatic = TRUE)
  expect(identical(dim(u$nodes), as.integer(c(n, d))), paste(id, "size"))
  expect(all(u$weights > 0) && near(sum(u$weights), 1), paste(id, "weights"))
  expect(exactness_error(u, p = 2 * q - 1, family = fam) < GATE, paste(id, "exactness"))
  parts <- ruleinfo(fam, d, q, pragmatic = TRUE)
  expect(sum(parts$d) == d && prod(parts$n) == n, paste(id, "ruleinfo"))
  expect(all(vapply(seq_len(d %/% 2), function(j)
    n <= nnodes(fam, j, q, pragmatic = TRUE) * nnodes(fam, d - j, q, pragmatic = TRUE), NA)), paste(id, "cheapest"))
}
u <- lepos(6, 2, normalize = FALSE, pragmatic = TRUE)
expect(near(sum(u$weights), 2^6) && near(sum(u$weights * u$nodes[, 6]^2), 2^6 / 3), "fallback with normalize = FALSE")
expect(nnodes("gh", 60, 16, pragmatic = TRUE) > 2^53, "huge counts are reported")
expect_error(ghpos(60, 16, pragmatic = TRUE), "huge rules are refused")

# --- the data file
dir <- if (is.null(getOption("quadriceps.datadir"))) system.file("extdata", package = "quadriceps") else getOption("quadriceps.datadir")
expect(identical(sort(list.files(dir, recursive = TRUE)), c("index.tsv", "rules.bin")), "one binary file and its catalog, nothing else")
n_rules <- nrow(available_rules("gh")) + nrow(available_rules("le"))
expect(startsWith(readLines(file.path(dir, "rules.bin"), n = 1L, warn = FALSE),
                  sprintf("QUADRICEPS1 fmt=1 endian=little cells=%d index_fields=8 float=binary64", n_rules)), "header of rules.bin")

cat(sprintf("%d checks, %d failed\n", checks, length(failures)))
if (length(failures) > 0) stop("failed:\n  ", paste(failures, collapse = "\n  "))
