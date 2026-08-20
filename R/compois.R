## ---------------------------------------------------------------------------
## Native Conway--Maxwell--Poisson (COM-Poisson): the soft-tail counterpart to
## the CPB's hard ceiling. Implemented as a count family so that count_reg(),
## hurdle_count(), and zi_count() give it fixed effects (one- and two-way),
## zero-truncation, hurdle, zero-inflation, offsets, and robust/clustered SEs
## through the same machinery as Poisson/NB; the Mundlak device (mundlak())
## supplies correlated random effects. Distribution functions d/p/q/r below.
##
## Parameterization: log(lambda) = x'beta with dispersion nu (the standard
## COM-Poisson rate parameterization). nu > 1 is underdispersion, nu = 1 is the
## Poisson, nu < 1 is overdispersion. The rate lambda is not the mean; the mean
## E(Y) is recovered numerically (see fam$meanfun), so predict()/first_difference()
## report effects on the mean, while the coefficients are on the log-rate scale.
## ---------------------------------------------------------------------------

## upper support for truncating the normalizing sum Z = sum_j lambda^j/(j!)^nu.
## The unnormalized log-terms j*log(lambda) - nu*lgamma(j+1) are unimodal in j;
## we grow the support past the mode until a term falls ~40 log-units below the
## peak (mass negligible), so the truncation is convergence-based rather than a
## fixed cap. This keeps Z correctly normalized for overdispersion (nu < 1),
## where the support can be large, as well as the underdispersed regime (nu > 1),
## where the mode is small and the loop stops almost immediately. A hard ceiling
## of 100000 bounds a transient extreme lambda during optimization.
.compois_kmax <- function(lambda, nu) {
  logl <- log(min(max(lambda, 1e-12), 1e10))                # bounded log-rate (guards optim excursions)
  mode <- min(max(1, exp(logl / max(nu, 0.05))), 1e5)       # bounded mode of the terms
  term <- function(k) k * logl - nu * lgamma(k + 1)
  peak <- term(mode); if (!is.finite(peak)) peak <- 0
  K <- as.integer(min(ceiling(mode), 100000L))
  step <- max(16L, as.integer(min(ceiling(mode * 0.25) + 1, 10000L)))
  while (K < 100000L) {
    tk <- term(K + step)
    if (!is.finite(tk) || tk < peak - 40) break
    K <- K + step
  }
  as.integer(max(50L, min(K + step, 100000L)))
}

## log normalizing constant, vectorized over lambda (scalar nu), via log-sum-exp
.compois_logZ <- function(lambda, nu, kmax = NULL) {
  K <- if (is.null(kmax)) .compois_kmax(lambda, nu) else kmax
  j <- 0:K; lj <- -nu * lgamma(j + 1)
  vapply(lambda, function(l) {
    lt <- j * log(l) + lj; m <- max(lt); m + log(sum(exp(lt - m)))
  }, numeric(1))
}

## pmf over 0..kmax for a single lambda (normalized by the full Z)
.compois_pmf1 <- function(lambda, nu, kmax) {
  logZ <- .compois_logZ(lambda, nu, kmax = max(kmax, .compois_kmax(lambda, nu)))
  j <- 0:kmax
  exp(j * log(lambda) - nu * lgamma(j + 1) - logZ)
}

## mean E(Y), vectorized over lambda (scalar nu)
.compois_mean <- function(lambda, nu) {
  vapply(lambda, function(l) {
    K <- .compois_kmax(l, nu); j <- 0:K
    p <- .compois_pmf1(l, nu, K); sum(j * p)
  }, numeric(1))
}

#' Conway--Maxwell--Poisson distribution functions
#'
#' Density, distribution, quantile, and random generation for the COM-Poisson
#' with rate `lambda` and dispersion `nu` (`nu > 1` underdispersed, `nu = 1`
#' Poisson, `nu < 1` overdispersed). Complements the estimator
#' `count_reg(..., family = "compois")`.
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param n Number of draws.
#' @param lambda Rate parameter (scalar or vector, recycled).
#' @param nu Dispersion parameter (scalar).
#' @param log,log.p Return log probabilities.
#' @param lower.tail If `TRUE` (default), \eqn{P(X \le x)}.
#' @return `dcompois` a density, `pcompois` a CDF, `qcompois` a quantile,
#'   `rcompois` a numeric vector of count draws.
#' @examples
#' dcompois(0:5, lambda = 3, nu = 1.5)
#' mean(rcompois(1000, lambda = 3, nu = 1.5))
#' @name compois-distribution
#' @export
dcompois <- function(x, lambda, nu, log = FALSE) {
  n <- max(length(x), length(lambda)); x <- rep_len(x, n); lambda <- rep_len(lambda, n)
  logZ <- .compois_logZ(lambda, nu)
  ld <- ifelse(x < 0 | x != floor(x), -Inf, x * log(lambda) - nu * lgamma(x + 1) - logZ)
  if (log) ld else exp(ld)
}
#' @rdname compois-distribution
#' @export
pcompois <- function(q, lambda, nu, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(lambda)); q <- rep_len(q, n); lambda <- rep_len(lambda, n)
  out <- vapply(seq_len(n), function(i) {
    qi <- floor(q[i]); if (qi < 0) return(0)
    sum(.compois_pmf1(lambda[i], nu, qi))
  }, numeric(1))
  out <- pmin(pmax(out, 0), 1); if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
#' @rdname compois-distribution
#' @export
qcompois <- function(p, lambda, nu, lower.tail = TRUE, log.p = FALSE) {
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(lambda)); p <- rep_len(p, n); lambda <- rep_len(lambda, n)
  vapply(seq_len(n), function(i) {
    K <- .compois_kmax(lambda[i], nu); cdf <- cumsum(.compois_pmf1(lambda[i], nu, K))
    as.numeric(which(cdf >= p[i] - 1e-12)[1L] - 1L)
  }, numeric(1))
}
#' @rdname compois-distribution
#' @export
rcompois <- function(n, lambda, nu) {
  lambda <- rep_len(lambda, n)
  vapply(seq_len(n), function(i) {
    K <- .compois_kmax(lambda[i], nu); cdf <- cumsum(.compois_pmf1(lambda[i], nu, K))
    as.numeric(which(cdf >= stats::runif(1))[1L] - 1L)
  }, numeric(1))
}
