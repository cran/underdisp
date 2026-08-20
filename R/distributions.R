## ---------------------------------------------------------------------------
## Distribution functions (d / p / q) for the CPB and the GEC (Katz) family.
## The CPB simulator rcpb() and the CPB mixtures rhurdle_cpb()/rzicpb() are in
## hurdle.R / zi.R. The zero-truncated / hurdle / zero-inflated Poisson and
## negative-binomial helpers below are INTERNAL (not exported): CRAN already
## serves these standard families (VGAM, an Import, provides dpospois, dzipois,
## and relatives); ours exist so the test suite and the replication scripts can
## simulate and score every DV in one idiom, reached via underdisp:::.
## ---------------------------------------------------------------------------

## --- Continuous Parameter Binomial: density, CDF, quantile ----------------

#' Continuous parameter binomial distribution functions
#'
#' Density, distribution function, and quantile function for the continuous
#' parameter binomial (CPB), and its zero-truncated form. The CPB has a hard
#' ceiling at \eqn{\lfloor \lambda/(1-\alpha) \rfloor}; probability above it is
#' zero. Complements the simulator [rcpb()].
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param lambda Mean parameter (scalar or vector, recycled).
#' @param alpha Shape parameter in (0, 1).
#' @param truncated If `TRUE`, use the zero-truncated CPB.
#' @param log,log.p If `TRUE`, probabilities are given as log.
#' @param lower.tail If `TRUE` (default), probabilities are \eqn{P(X \le x)}.
#' @return `dcpb` a density, `pcpb` a distribution function, `qcpb` a quantile.
#' @examples
#' dcpb(0:5, lambda = 3, alpha = 0.5)
#' pcpb(3, lambda = 3, alpha = 0.5)
#' qcpb(0.9, lambda = 3, alpha = 0.5)
#' @seealso [rcpb()]
#' @name cpb-distribution
#' @export
dcpb <- function(x, lambda, alpha, truncated = FALSE, log = FALSE) {
  n <- max(length(x), length(lambda)); x <- rep_len(x, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    xi <- x[i]
    if (xi < 0 || xi != floor(xi)) return(0)
    pm <- .cpb_pmf1(lambda[i], alpha, kmax = xi, truncated = truncated)
    pm[xi + 1L]
  }, numeric(1))
  if (log) log(out) else out
}

#' @rdname cpb-distribution
#' @export
pcpb <- function(q, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(lambda)); q <- rep_len(q, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    qi <- floor(q[i]); if (qi < 0) return(0)
    sum(.cpb_pmf1(lambda[i], alpha, kmax = qi, truncated = truncated))
  }, numeric(1))
  out <- pmin(pmax(out, 0), 1)
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname cpb-distribution
#' @export
qcpb <- function(p, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(lambda)); p <- rep_len(p, n); lambda <- rep_len(lambda, n)
  vapply(seq_len(n), function(i) {
    K <- floor(lambda[i] / (1 - alpha))
    cdf <- cumsum(.cpb_pmf1(lambda[i], alpha, kmax = K, truncated = truncated))
    as.numeric(which(cdf >= p[i] - 1e-12)[1L] - 1L)
  }, numeric(1))
}

## --- zero-truncated Poisson / negative binomial ---------------------------
## Inverse-CDF sampling: draw U on (F(0), 1) and invert the base quantile, which
## is exact and vectorized (no rejection loop).

dztpois <- function(x, lambda, log = FALSE) {
  d <- ifelse(x < 1 | x != floor(x), 0, stats::dpois(x, lambda) / -expm1(-lambda))  # -expm1 avoids 1-1 cancellation
  if (log) log(d) else d
}
rztpois <- function(n, lambda) {
  lambda <- rep_len(lambda, n); p0 <- exp(-lambda)
  stats::qpois(p0 + stats::runif(n) * (1 - p0), lambda)
}
dztnbinom <- function(x, mu, size, log = FALSE) {
  logp0 <- size * (log(size) - log(size + mu))                       # log P(0), stable for tiny mu
  d <- ifelse(x < 1 | x != floor(x), 0, stats::dnbinom(x, mu = mu, size = size) / -expm1(logp0))
  if (log) log(d) else d
}
rztnbinom <- function(n, mu, size) {
  mu <- rep_len(mu, n); p0 <- stats::dnbinom(0, mu = mu, size = size)
  stats::qnbinom(p0 + stats::runif(n) * (1 - p0), mu = mu, size = size)
}

## --- hurdle Poisson / negative binomial -----------------------------------

## Hurdle count: zero w.p. 1-p, otherwise a zero-truncated Poisson/NB draw.
dhpois <- function(x, mu, p, log = FALSE) {
  d <- ifelse(x == 0, 1 - p, ifelse(x >= 1 & x == floor(x), p * dztpois(x, mu), 0))
  if (log) log(d) else d
}
rhpois <- function(n, mu, p) {
  d <- stats::rbinom(n, 1, p); y <- integer(n); i <- d == 1
  if (any(i)) y[i] <- rztpois(sum(i), rep_len(mu, n)[i]); y
}
dhnbinom <- function(x, mu, p, size, log = FALSE) {
  d <- ifelse(x == 0, 1 - p, ifelse(x >= 1 & x == floor(x), p * dztnbinom(x, mu, size), 0))
  if (log) log(d) else d
}
rhnbinom <- function(n, mu, p, size) {
  d <- stats::rbinom(n, 1, p); y <- integer(n); i <- d == 1
  if (any(i)) y[i] <- rztnbinom(sum(i), rep_len(mu, n)[i], size); y
}

## --- zero-inflated Poisson / negative binomial ----------------------------

## Zero-inflated count: structural zero w.p. pi, otherwise a Poisson/NB draw.
dzipois <- function(x, mu, pi, log = FALSE) {
  d <- ifelse(x == 0, pi + (1 - pi) * stats::dpois(0, mu),
              ifelse(x >= 1 & x == floor(x), (1 - pi) * stats::dpois(x, mu), 0))
  if (log) log(d) else d
}
rzipois <- function(n, mu, pi) ifelse(stats::rbinom(n, 1, pi) == 1, 0L, stats::rpois(n, mu))
dzinbinom <- function(x, mu, pi, size, log = FALSE) {
  d <- ifelse(x == 0, pi + (1 - pi) * stats::dnbinom(0, mu = mu, size = size),
              ifelse(x >= 1 & x == floor(x), (1 - pi) * stats::dnbinom(x, mu = mu, size = size), 0))
  if (log) log(d) else d
}
rzinbinom <- function(n, mu, pi, size)
  ifelse(stats::rbinom(n, 1, pi) == 1, 0L, stats::rnbinom(n, mu = mu, size = size))

## --- generalized event count (Katz) family: density, CDF, quantile, random --
## The GEC mean is the rate lambda and its variance-to-mean ratio is delta
## (delta < 1 underdispersed, = 1 Poisson, > 1 overdispersed). Uses the same
## compiled pmf (gec_pmf_cpp) as the estimator, so d/r are exactly consistent.

.gec_kmax_qr <- function(lambda, delta, max.support) {
  min(as.integer(max.support),
      as.integer(ceiling(max(lambda) + 12 * sqrt(max(lambda) * max(delta, 0.1)) + 30)))
}

#' Generalized event count (Katz) distribution functions
#'
#' Density, distribution, quantile, and random generation for the King (1989)
#' generalized event count model with rate `lambda` (the mean) and dispersion
#' `delta` (the variance-to-mean ratio: `< 1` underdispersed, `= 1` Poisson,
#' `> 1` overdispersed). Consistent with the estimator [gec()].
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param n Number of draws.
#' @param lambda Rate/mean parameter (scalar or vector, recycled).
#' @param delta Dispersion (variance-to-mean ratio).
#' @param max.support Guard on the evaluated support.
#' @param log,log.p Return log probabilities.
#' @param lower.tail If `TRUE` (default), \eqn{P(X \le x)}.
#' @return `dgec` a density, `pgec` a CDF, `qgec` a quantile, `rgec` a numeric
#'   vector of count draws.
#' @examples
#' dgec(0:5, lambda = 3, delta = 0.7)
#' var(rgec(2000, lambda = 3, delta = 0.7)) / mean(rgec(2000, lambda = 3, delta = 0.7))
#' @seealso [gec()]
#' @name gec-distribution
#' @export
dgec <- function(x, lambda, delta, max.support = 500, log = FALSE) {
  n <- max(length(x), length(lambda)); x <- rep_len(x, n); lambda <- rep_len(lambda, n)
  km <- max(0L, as.integer(max(x)))
  P <- gec_pmf_cpp(lambda, delta, km, as.integer(max.support))
  d <- vapply(seq_len(n), function(i) if (x[i] < 0 || x[i] != floor(x[i])) 0 else P[i, x[i] + 1L], numeric(1))
  if (log) log(d) else d
}
#' @rdname gec-distribution
#' @export
pgec <- function(q, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(lambda)); q <- rep_len(q, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    qi <- floor(q[i]); if (qi < 0) return(0)
    sum(gec_pmf_cpp(lambda[i], delta, as.integer(qi), as.integer(max.support)))
  }, numeric(1))
  out <- pmin(pmax(out, 0), 1); if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname gec-distribution
#' @export
qgec <- function(p, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(lambda)); p <- rep_len(p, n); lambda <- rep_len(lambda, n)
  vapply(seq_len(n), function(i) {
    km <- .gec_kmax_qr(lambda[i], delta, max.support)
    cdf <- cumsum(as.numeric(gec_pmf_cpp(lambda[i], delta, km, as.integer(max.support))))
    idx <- which(cdf >= p[i] - 1e-12)[1L]; as.numeric(if (is.na(idx)) km else idx - 1L)
  }, numeric(1))
}
#' @rdname gec-distribution
#' @export
rgec <- function(n, lambda, delta, max.support = 500) {
  lambda <- rep_len(lambda, n)
  vapply(seq_len(n), function(i) {
    km <- .gec_kmax_qr(lambda[i], delta, max.support)
    cdf <- cumsum(as.numeric(gec_pmf_cpp(lambda[i], delta, km, as.integer(max.support))))
    idx <- which(cdf >= stats::runif(1))[1L]; as.numeric(if (is.na(idx)) km else idx - 1L)
  }, numeric(1))
}
