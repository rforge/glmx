.ab_to_scale <- function(a, b)
  (qgl(3/4, c(0, 1, a - b, a + b)) -  qgl(1/4, c(0, 1, a - b, a + b)))/2.197224

.ab_to_lambda <- function(a, b)
  c(0, 1, a - b, a + b)

qpregibon <- function(p, a = 0, b = 0)
{
  stopifnot(require("gld"))
  qgl(p, .ab_to_lambda(a, b))/.ab_to_scale(a, b)
}

ppregibon <- function(q, a = 0, b = 0, tol = 1e-12)
{
  stopifnot(require("gld"))
  pgl(q * .ab_to_scale(a, b), .ab_to_lambda(a, b), inverse.eps = tol)
}

dpregibon <- function(x, a = 0, b = 0, tol = 1e-12)
{
  stopifnot(require("gld"))
  s <- .ab_to_scale(a, b)
  dgl(x * s, .ab_to_lambda(a, b), inverse.eps = tol) * s
}

rpregibon <- function(n, a = 0, b = 0)
{
  stopifnot(require("gld"))
  qpregibon(runif(n), a = a, b = b)
}
