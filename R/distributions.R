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
#' ceiling: King's CPB (the default) puts no probability above
#' \eqn{\lceil \lambda/(1-\alpha) \rceil}, the floor support (`support =
#' "floor"`, the definition before version 0.2.0) none above
#' \eqn{\lfloor \lambda/(1-\alpha) \rfloor}. Complements the simulator [rcpb()].
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param lambda Rate parameter (scalar or vector, recycled). It is the mean
#'   where `lambda / (1 - alpha)` is an integer; otherwise the mean lies above it
#'   on King's support and below it on the floor support.
#' @param alpha Shape parameter in (0, 1).
#' @param truncated If `TRUE`, use the zero-truncated CPB. On the floor support,
#'   where the implied ceiling `lambda / (1 - alpha)` is below 1 the CPB puts all
#'   its mass on 0, and its zero-truncated form is the point mass at 1, the
#'   zero-truncated distribution for every ceiling in [1, 2); King's support
#'   always contains 1.
#' @param support The CPB's support, `"ceiling"` (King's CPB, the default) or
#'   `"floor"`; see [cpb()].
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
dcpb <- function(x, lambda, alpha, truncated = FALSE, log = FALSE, support = c("ceiling", "floor")) {
  support <- match.arg(support)
  if (!length(x) || !length(lambda)) return(numeric(0))
  if (length(alpha) != 1L) stop("'alpha' must be a single value.")
  n <- max(length(x), length(lambda)); x <- rep_len(x, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    xi <- x[i]
    if (xi < 0 || xi != floor(xi)) return(0)
    pm <- .cpb_pmf1(lambda[i], alpha, kmax = xi, truncated = truncated, support = support)
    pm[xi + 1L]
  }, numeric(1))
  out <- .cpb_warn_cap(out)
  if (log) log(out) else out
}

## a CPB pmf refused at the 1e5-term cap (a rate of order 1e5 on King's support) is NA, with a warning
.cpb_warn_cap <- function(out) {
  if (anyNA(out)) warning("the CPB support reaches the 1e5-term cap for ", sum(is.na(out)), " value(s) (a rate of ",
                          "order 1e5): NA returned.", call. = FALSE)
  out
}

#' @rdname cpb-distribution
#' @export
pcpb <- function(q, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE,
                 support = c("ceiling", "floor")) {
  support <- match.arg(support)
  if (!length(q) || !length(lambda)) return(numeric(0))
  if (length(alpha) != 1L) stop("'alpha' must be a single value.")
  n <- max(length(q), length(lambda)); q <- rep_len(q, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    qi <- floor(q[i]); if (qi < 0) return(0)
    sum(.cpb_pmf1(lambda[i], alpha, kmax = qi, truncated = truncated, support = support))
  }, numeric(1))
  out <- pmin(pmax(.cpb_warn_cap(out), 0), 1)
  if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}

#' @rdname cpb-distribution
#' @export
qcpb <- function(p, lambda, alpha, truncated = FALSE, lower.tail = TRUE, log.p = FALSE,
                 support = c("ceiling", "floor")) {
  support <- match.arg(support)
  if (!length(p) || !length(lambda)) return(numeric(0))
  if (length(alpha) != 1L) stop("'alpha' must be a single value.")
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(lambda)); p <- rep_len(p, n); lambda <- rep_len(lambda, n)
  .cpb_warn_cap(vapply(seq_len(n), function(i) {
    K <- if (.cpb_fl(support)) floor(lambda[i] / (1 - alpha) + 1e-9)          # tolerance matches the pmf core
         else { pr <- .cpb_pmf_core(lambda[i], alpha, support)                 # the support, to negligible mass
                if (anyNA(pr)) return(NA_real_); length(pr) - 1L }
    cdf <- cumsum(.cpb_pmf1(lambda[i], alpha, kmax = K, truncated = truncated, support = support))
    as.numeric(which(cdf >= p[i] - 1e-12)[1L] - 1L)
  }, numeric(1)))
}

## --- zero-truncated Poisson / negative binomial densities -------------------

dztpois <- function(x, lambda, log = FALSE) {
  d <- ifelse(x < 1 | x != floor(x), 0, stats::dpois(x, lambda) / -expm1(-lambda))  # -expm1 avoids 1-1 cancellation
  if (log) log(d) else d
}
dztnbinom <- function(x, mu, size, log = FALSE) {
  logp0 <- size * (log(size) - log(size + mu))                       # log P(0), stable for tiny mu
  d <- ifelse(x < 1 | x != floor(x), 0, stats::dnbinom(x, mu = mu, size = size) / -expm1(logp0))
  if (log) log(d) else d
}

## --- generalized event count (Katz) family: density, CDF, quantile, random --
## On an unbounded support (delta >= 1) the GEC mean is the rate lambda and its
## variance-to-mean ratio is delta; for delta < 1 the support is finite and the
## exact moments differ (see dgec's help). Uses the same
## compiled pmf (gec_pmf_cpp) as the estimator, so d/r are exactly consistent.

.gec_kmax_qr <- function(lambda, delta, max.support) {
  min(as.integer(max.support),
      as.integer(ceiling(max(lambda) + 12 * sqrt(max(lambda) * max(delta, 0.1)) + 30)))
}
## the smallest count whose cdf reaches u: evaluated to the rough cut above, then to the guard when an
## overdispersed tail reaches past it; NA when the pmf is unavailable (the guard binds at this rate)
.gec_quantile1 <- function(u, lambda, delta, max.support, tol = 1e-12) {
  for (km in unique(c(.gec_kmax_qr(lambda, delta, max.support), as.integer(max.support)))) {
    pr <- as.numeric(gec_pmf_cpp(lambda, delta, km, as.integer(max.support)))
    if (anyNA(pr)) return(NA_real_)
    idx <- which(cumsum(pr) >= u - tol)[1L]
    if (!is.na(idx)) return(as.numeric(idx - 1L))
  }
  as.numeric(max.support)
}
.gec_warn_guard <- function(out, max.support) {
  if (anyNA(out)) warning("'max.support' (", max.support, ") binds for ", sum(is.na(out)), " rate(s): NA returned; raise it.",
                          call. = FALSE)
  out
}

#' Generalized event count (Katz) distribution functions
#'
#' Density, distribution, quantile, and random generation for the King (1989)
#' generalized event count model with rate `lambda` and Katz dispersion
#' `delta` (`< 1` underdispersed, `= 1` Poisson, `> 1` overdispersed).
#' On an unbounded support (`delta >= 1`) the mean is `lambda` and the
#' variance-to-mean ratio is `delta` exactly; for `delta < 1` the support is
#' finite and the renormalized distribution's mean and variance equal those
#' values only when `lambda/(1 - delta)` is an integer (the exact moments are
#' finite sums of the pmf, as `count_reg()`'s siblings report them).
#' Consistent with the estimator [gec()].
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param n Number of draws.
#' @param lambda Rate parameter (scalar or vector, recycled).
#' @param delta Katz dispersion parameter (scalar).
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
  .ud_check_whole(max.support, "max.support")
  if (!length(x) || !length(lambda)) return(numeric(0))
  if (length(delta) != 1L) stop("'delta' must be a single value.")
  n <- max(length(x), length(lambda)); x <- rep_len(x, n); lambda <- rep_len(lambda, n)
  km <- max(0L, as.integer(max(x)))
  P <- gec_pmf_cpp(lambda, delta, km, as.integer(max.support))
  d <- vapply(seq_len(n), function(i) if (x[i] < 0 || x[i] != floor(x[i])) 0 else P[i, x[i] + 1L], numeric(1))
  if (anyNA(d)) warning("'max.support' (", max.support, ") binds for ", sum(is.na(d)), " rate(s): NA returned; raise it.", call. = FALSE)
  if (log) log(d) else d
}
#' @rdname gec-distribution
#' @export
pgec <- function(q, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE) {
  .ud_check_whole(max.support, "max.support")
  if (!length(q) || !length(lambda)) return(numeric(0))
  if (length(delta) != 1L) stop("'delta' must be a single value.")
  n <- max(length(q), length(lambda)); q <- rep_len(q, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    qi <- floor(q[i]); if (qi < 0) return(0)
    sum(gec_pmf_cpp(lambda[i], delta, as.integer(qi), as.integer(max.support)))
  }, numeric(1))
  if (anyNA(out)) warning("'max.support' (", max.support, ") binds for ", sum(is.na(out)), " rate(s): NA returned; raise it.", call. = FALSE)
  out <- pmin(pmax(out, 0), 1); if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname gec-distribution
#' @export
qgec <- function(p, lambda, delta, max.support = 500, lower.tail = TRUE, log.p = FALSE) {
  .ud_check_whole(max.support, "max.support")
  if (!length(p) || !length(lambda)) return(numeric(0))
  if (length(delta) != 1L) stop("'delta' must be a single value.")
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(lambda)); p <- rep_len(p, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) .gec_quantile1(p[i], lambda[i], delta, max.support), numeric(1))
  .gec_warn_guard(out, max.support)
}
#' @rdname gec-distribution
#' @export
rgec <- function(n, lambda, delta, max.support = 500) {
  .ud_check_whole(max.support, "max.support")
  if (n == 0) return(numeric(0))
  if (length(delta) != 1L) stop("'delta' must be a single value.")
  lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) .gec_quantile1(stats::runif(1), lambda[i], delta, max.support, tol = 0),
                numeric(1))
  .gec_warn_guard(out, max.support)
}
