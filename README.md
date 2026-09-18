# quadriceps (R)

> **Paper:** J. Pinkse, *Positive weight Hermite and Legendre quadrature rules* — arXiv: **[ARXIV-LINK-TBA](https://arxiv.org/abs/ARXIV-LINK-TBA)** (link to be filled in on publication)
>
> **Data deposit:** Zenodo — DOI: **[ZENODO-DOI-TBA](https://doi.org/ZENODO-DOI-TBA)** (link to be filled in on publication)

Positive-weight cubature rules in several dimensions, for two weights:

| function | weight (default) | one-dimensional cousin |
|---|---|---|
| `ghpos(d, q)` | standard normal density `N(0, I_d)` on `R^d` | `statmod::gauss.quad(q, "hermite")` |
| `lepos(d, q)` | uniform density on `[0,1]^d` | `statmod::gauss.quad(q, "legendre")` |

A rule of degree `p` is a set of `n` nodes `x_i` in `R^d` and weights `w_i > 0` with
`sum_i w_i f(x_i) = integral of f(x) ω(x) dx` for every polynomial `f` of total degree `<= p`.
The product of `q`-node one-dimensional Gauss rules does this for `p = 2q - 1` with `q^d` nodes.
The rules stored here, the smallest positive-weight rules known to the author, do it with far
fewer; at `d = 5` the saving is more than a factor of ten. They cover `2 <= d <= 5`.
[`RULES.md`](RULES.md) lists every rule with its node count, Möller's lower bound, measured
accuracy and origin.

This is the R twin of [Quadriceps.jl](https://github.com/NittanyLion/Quadriceps.jl) (Julia) and
[quadriceps-py](https://github.com/NittanyLion/quadriceps-py) (Python). The three packages share
their data, their conventions and their function names; the data are refreshed from the Julia
package whenever a smaller rule is found.

## Installation

The repository is private. With access to it, either

```r
remotes::install_git("git@github.com:NittanyLion/quadriceps-r.git")
```

or clone it and run `R CMD INSTALL quadriceps-r`. The package needs base R only.

## Use

```r
library(quadriceps)

r <- ghpos(3, 4)                # d = 3, q = 4 (degree 7): r$nodes is 27 x 3, r$weights has length 27
f <- function(x) x[, 1]^2 * x[, 2]^4
sum(r$weights * f(r$nodes))     # E[Z1^2 Z2^4] = 3

r <- lepos(2, 5)                # q = 5 (degree 9): 17 nodes on the unit square instead of 25
sum(r$weights * r$nodes[, 1]^3 * r$nodes[, 2]^2)      # 1/4 * 1/3

r <- ghpos(3, p = 7)            # the first rule again, requested by its degree
```

Both functions return a list with components `nodes` (an `n x d` matrix, one node per row) and
`weights` (`n` positive weights), like `statmod::gauss.quad`. The second argument `q` has its
one-dimensional meaning: the number of nodes of the one-dimensional Gauss rule, which is exact
to degree `2q - 1`. `ghpos(d, q)` returns a `d`-dimensional rule of that same degree
`p = 2q - 1`: a replacement for the `q^d`-node product grid, and for `d = 1` the `q`-node Gauss
rule itself.

To ask for a degree instead, name the argument `p`: `ghpos(d, p = 7)`, `lepos(d, p = 12)`. Any
`p >= 0` is accepted. Rules are stored at odd degrees and a request is served by the smallest
stored rule of degree `>= p`, so an even `p` gets the rule for `p + 1`. Give `q` or `p`, not
both.

### `normalize`

`statmod::gauss.quad(q, "hermite")` integrates against `exp(-x^2)`. The rules here are made for
the standard normal density, which is what an expectation needs, so `normalize = TRUE` is the
default:

| | `normalize = TRUE` (default) | `normalize = FALSE` (`gauss.quad`'s convention) |
|---|---|---|
| `ghpos` | weight `(2 pi)^(-d/2) exp(-‖x‖^2/2)`; weights sum to 1 | weight `exp(-‖x‖^2)`; weights sum to `pi^(d/2)` |
| `lepos` | uniform density on `[0,1]^d`; weights sum to 1 | plain integral over `[-1,1]^d`; weights sum to `2^d` |

For `Y ~ N(mu, L L')` use the nodes `sweep(r$nodes %*% t(L), 2, mu, "+")` with the same weights;
for a box, rescale the columns of the Le nodes.

### `pragmatic`

Rules are stored for `2 <= d <= 5`, up to a degree that depends on the family and on `d` (see
[`RULES.md`](RULES.md) or `available_rules()`). For any other request:

* `pragmatic = FALSE` (the default) signals an error of class `quadriceps_no_rule`, whose
  message says how far the stored rules go;
* `pragmatic = TRUE` returns the cheapest tensor product of lower-dimensional rules: the split
  of `d` into stored rules and one-dimensional Gauss rules that needs the fewest nodes. The
  result is a valid positive-weight rule of the requested degree. It is not small, but it is
  much smaller than the plain product grid whenever a stored rule can be a factor.

```r
ghpos(7, 5)                             # error: no stored rule in seven dimensions
r <- ghpos(7, 5, pragmatic = TRUE)      # (d = 2) x (d = 5): a few thousand nodes; the grid has 78125

nnodes("gh", 10, 3, pragmatic = TRUE)   # the node count, without building the rule
ruleinfo("gh", 7, 5, pragmatic = TRUE)  # the factors, with their origins
```

With `pragmatic = TRUE` a request that a stored rule covers returns that stored rule, as without
it. (If a product were ever strictly cheaper than the stored rule, the product would be
returned; the data contain no such case, and the tests check that.)

### Other functions

* `available_rules(family)` lists the stored rules as a data frame; `family` is `"gh"` or `"le"`.
* `nnodes(family, d, q, p, pragmatic)` gives a node count without building the rule.
* `ruleinfo(family, d, q, p, pragmatic)` describes the rule and its origin.
* `exactness_error(nodes, weights, p, family)` measures how exact a rule is.

Every function has a help page (`?ghpos`, `?quadriceps`).

All rules are stored in one binary file, `inst/extdata/rules.bin`, with `index.tsv` next to it
as the catalog; [`FORMAT.md`](FORMAT.md) specifies the format, which the Julia and Python
packages share byte for byte.

## Accuracy

Rules are stored in double precision. Every stored rule was checked when the data were built:
all weights positive, and the largest relative monomial error over all monomials of degree
`<= p` below `1e-11`. Most rules sit at `1e-16` to `1e-15`; the largest GH rules at `d = 2, 3`
are the least accurate, at `1e-12` to `1e-11`. The measured value of each rule is in the catalog
(`relerr`) and in [`RULES.md`](RULES.md), and the tests repeat the check for every rule. All
nodes of the Le rules lie strictly inside the cube. One-dimensional Gauss rules are computed by
the Golub–Welsch eigenvalue method in base R.

## Whose rules these are

Most rules were computed by the author. Some are rules from the literature (copied in, or
found again by the author's search and recognized), and some Le rules were obtained by node
elimination started from Diallo and Worku's published rules. `ruleinfo()` and the `origin`
column of [`RULES.md`](RULES.md) say which is which; cite the source named there when you use
such a rule. [`NOTICE.md`](NOTICE.md) has the details and the license notice that travels with
the derived files.

## Tests

```
Rscript tests/run_tests.R       # from the package root; about 20 s
```
