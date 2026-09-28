## ---------------------------------------------------------------------------
## Count families for count_reg(), hurdle_count(), and zi_count().
##
## One family object per distribution. Each carries the pmf on the log scale
## (`logpmf`), the pmf vector over 0..kmax (`pvec`), P(Y = 0) (`p0`), the mean
## and variance at the natural parameter (`meanfun`, `varfun`), the support used
## by simulate() (`kmax`), and the shape parameter's parameterization for the
## optimizer (`shape_link`, `shape_inv`, `shape_starts`). Estimation, scoring,
## prediction, simulation, and every two-part variant consume only these
## fields, so a family added here inherits fixed effects (one- and two-way),
## zero-truncation, hurdle and zero-inflated forms, offsets, weights,
## analytic / robust / cluster standard errors, proper scores, DHARMa
## simulation, and the QoI methods.
##
## Natural parameter `mu`: the mean for poisson, negbin, mpcmp (Huang's
## mean-parameterized COM-Poisson), genpois, gammacount, and doublepois (for the
## last two the mean is the large-sample rate; the exact mean is returned by
## `meanfun`), and the rate lambda for the classical COM-Poisson ("compois").
##
## Underdispersed members and their dispersion parameters:
##   compois     nu > 1        classical COM-Poisson, rate-parameterized
##   mpcmp       nu > 1        COM-Poisson, mean-parameterized (Huang 2017)
##   genpois     lambda < 0    Consul-Jain generalized Poisson, constant lambda,
##                             Var/Mean = 1/(1-lambda)^2; finite support when
##                             lambda < 0, renormalized on that support
##   gammacount  alpha > 1     Winkelmann's (1995) gamma-count renewal process
##                             (regular event timing), Var/Mean ~ 1/alpha
##   doublepois  theta > 1     Efron's (1986) double Poisson with the exact
##                             normalizing constant, Var/Mean ~ 1/theta
## The Poisson is the equidispersed member of every family (nu = 1, lambda = 0,
## alpha = 1, theta = 1) and the negative binomial is the overdispersed default.
## ---------------------------------------------------------------------------

## log-sum-exp
.lse <- function(v) { v <- v[is.finite(v)]; if (!length(v)) return(-Inf); m <- max(v); m + log(sum(exp(v - m))) }

## upper support for a unimodal log-term sequence: grow past `K0` until the
## term falls 40 log units below the running maximum (or the cap binds)
.grow_support <- function(lterm, K0, cap = 100000L) {
  if (!is.finite(K0)) return(cap)
  K <- max(50L, as.integer(min(K0, cap))); step <- max(16L, as.integer(min(ceiling(K * 0.25), 10000)))
  lt <- lterm(0:K); mx <- max(lt[is.finite(lt)])
  while (K < cap) {
    tk <- lterm(K + step)
    if (!is.finite(tk) || tk < mx - 40) break
    K <- K + step
  }
  min(K + step, cap)
}

## ---- gamma-count (Winkelmann 1995) ----------------------------------------
## Waiting times iid Gamma(alpha, beta) with beta = alpha * mu, so the long-run
## event rate is mu. Events in (0, 1]: P(Y = y) = G(alpha y, alpha mu) -
## G(alpha (y + 1), alpha mu), G(a, x) the regularized lower incomplete gamma
## function with G(0, x) = 1. alpha > 1 underdispersion, alpha = 1 Poisson.
.gc_G <- function(a, x, lower = TRUE) {
  out <- stats::pgamma(x, shape = pmax(a, 1e-300), lower.tail = lower)
  out[a <= 0] <- if (lower) 1 else 0
  out
}
.gc_pmf <- function(y, mu, alpha) {
  x <- alpha * mu
  lo <- .gc_G(alpha * y, x)
  ## lower tails cancel near 1; switch to upper tails there
  d <- ifelse(lo < 0.5, lo - .gc_G(alpha * (y + 1), x),
              .gc_G(alpha * (y + 1), x, lower = FALSE) - .gc_G(alpha * y, x, lower = FALSE))
  d[y < 0 | y != floor(y)] <- 0
  pmax(d, 0)
}
.gc_kmax <- function(mu, alpha) {
  k <- mu + 12 * sqrt(mu / min(alpha, 1) + 1) + 30
  as.integer(if (!is.finite(k)) 100000 else min(100000, ceiling(k)))
}
.gc_moments <- function(mu, alpha) {
  vapply(mu, function(m) {
    K <- .gc_kmax(m, alpha); k <- 0:K; p <- .gc_pmf(k, m, alpha)
    s <- sum(p); e1 <- sum(k * p) / s; c(e1, sum(k * k * p) / s - e1 * e1)
  }, numeric(2))
}

## ---- double Poisson (Efron 1986) -------------------------------------------
## Unnormalized log-density  0.5 log theta - theta mu + [-y + y log y - log y!]
## + theta y (1 + log mu - log y)  (the bracket and the last term vanish at
## y = 0); the normalizing constant is the exact sum, not Efron's approximation.
.dp_lterm <- function(k, mu, theta) {
  lk <- log(pmax(k, 1))
  t1 <- ifelse(k == 0, 0, -k + k * lk - lgamma(k + 1))
  t2 <- ifelse(k == 0, 0, theta * k * (1 + log(mu) - lk))
  0.5 * log(theta) - theta * mu + t1 + t2
}
.dp_kmax <- function(mu, theta)
  .grow_support(function(k) .dp_lterm(k, mu, theta), ceiling(mu + 12 * sqrt(mu / theta + 1) + 30))
## the normalizer and the exact moments are summed in C++ (src/family_norm.cpp),
## one pass per observation with the same terms and support rule
.dp_logZ <- function(mu, theta) dp_norm_cpp(as.numeric(mu), theta, 100000L)[, 1]
.dp_moments <- function(mu, theta) {
  m <- dp_norm_cpp(as.numeric(mu), theta, 100000L)
  rbind(m[, 2], m[, 3])
}

## ---- generalized Poisson (Consul & Jain 1973), constant lambda ----------------
## P(y) = th (th + lambda y)^(y-1) exp(-th - lambda y) / y!, th_i = mu_i (1 - lambda),
## so E(Y) = mu and Var(Y) = mu / (1 - lambda)^2 on the untruncated support.
## For lambda < 0 the terms are positive only while th + lambda y > 0, so the
## support is finite (y <= K = ceiling(th / -lambda) - 1); the pmf is
## renormalized on that support so the likelihood is proper, and the exact mean
## and variance are computed from the normalized pmf (they equal mu and
## mu/(1-lambda)^2 up to the mass beyond the support, which is negligible away
## from the feasibility boundary).
.gp_lterm <- function(k, mu, lambda) {
  th <- mu * (1 - lambda); a <- th + lambda * k
  out <- log(th) + (k - 1) * log(pmax(a, 1e-300)) - th - lambda * k - lgamma(k + 1)
  out[a <= 0] <- -Inf
  out
}
.gp_kmax <- function(mu, lambda) {
  ## the terms are unimodal in k, so the tail rule applies for every lambda; for
  ## lambda < 0 the finite support end (beyond which the terms are -Inf) caps it
  K <- .grow_support(function(k) .gp_lterm(k, mu, lambda), ceiling(mu + 12 * sqrt(mu / (1 - lambda)^2 + 1) + 30))
  if (lambda < 0) K <- min(K, as.integer(max(0, ceiling(mu * (1 - lambda) / (-lambda)) - 1)))
  as.integer(K)
}
.gp_logZ <- function(mu, lambda) gp_norm_cpp(as.numeric(mu), lambda, 100000L)[, 1]
.gp_moments <- function(mu, lambda) {
  if (lambda >= 0) return(rbind(mu, mu / (1 - lambda)^2))
  m <- gp_norm_cpp(as.numeric(mu), lambda, 100000L)
  rbind(m[, 2], m[, 3])
}

## ---- mean-parameterized COM-Poisson (Huang 2017) ---------------------------
## For every observation the rate solving E[Y | lambda_i, nu] = mu_i is found in
## C++ (cmp_loglambda_cpp); the pmf is then the classical COM-Poisson pmf at
## that rate. The mean is mu by construction; the variance is exact.
.mp_lterm <- function(k, loglam, nu) k * loglam - nu * lgamma(k + 1)
.mp_kmax <- function(mu, nu) {
  v <- cmp_moments_cpp(cmp_loglambda_cpp(mu, nu), nu)[, 2]
  k <- mu + 12 * sqrt(v + 1) + 30
  as.integer(if (!is.finite(k)) 100000 else min(100000, ceiling(k)))
}

## ---- the family table ---------------------------------------------------------
.count_families <- list(
  poisson = list(tag = "poisson", label = "Poisson", nshape = 0L,
    shape_name = NULL, shape_link = identity, shape_inv = identity, shape_starts = numeric(0),
    shape_line = function(theta) "",
    null_shape = NA_real_,
    logpmf  = function(y, mu, theta) stats::dpois(y, mu, log = TRUE),
    pvec    = function(mu, theta, kmax) stats::dpois(0:kmax, mu),
    p0      = function(mu, theta) exp(-mu),
    meanfun = function(mu, theta) mu,
    varfun  = function(mu, theta) mu,
    kmax    = function(mu, theta) stats::qpois(1 - 1e-10, mu)),
  negbin = list(tag = "negbin", label = "NegBinomial", nshape = 1L,
    shape_name = "log(theta)", shape_link = log, shape_inv = exp, shape_starts = c(1, 1e4),
    shape_line = function(theta) sprintf("\ntheta (NB size): %.4f%s\n", theta,
                                         if (is.finite(theta) && theta > 1e6) "  [at the Poisson boundary: no overdispersion]" else ""),
    null_shape = Inf,
    logpmf  = function(y, mu, theta) stats::dnbinom(y, mu = mu, size = theta, log = TRUE),
    pvec    = function(mu, theta, kmax) stats::dnbinom(0:kmax, mu = mu, size = theta),
    p0      = function(mu, theta) stats::dnbinom(0, mu = mu, size = theta),
    meanfun = function(mu, theta) mu,
    varfun  = function(mu, theta) mu + mu^2 / theta,
    kmax    = function(mu, theta) stats::qnbinom(1 - 1e-10, mu = mu, size = theta)),
  compois = list(tag = "compois", label = "COM-Poisson", nshape = 1L,
    shape_name = "log(nu)", shape_link = log, shape_inv = exp, shape_starts = c(1, 2),
    shape_line = function(theta) sprintf("\nnu (COM-Poisson dispersion; > 1 underdispersed): %.4f\n", theta),
    null_shape = 1,
    log_mean = function(eta, theta) eta / theta,     # log lambda^(1/nu), the mean scale (see .count_excursion)
    logpmf  = function(y, mu, theta) y * log(mu) - theta * lgamma(y + 1) - .compois_logZ(mu, theta),
    pvec    = function(mu, theta, kmax) .compois_pmf1(mu, theta, kmax),
    p0      = function(mu, theta) exp(-.compois_logZ(mu, theta)),
    meanfun = function(mu, theta) .compois_mean(mu, theta),
    varfun  = function(mu, theta) cmp_moments_cpp(log(mu), theta)[, 2],
    kmax    = function(mu, theta) .compois_kmax(mu, theta)),
  mpcmp = list(tag = "mpcmp", label = "COM-Poisson (mean-parameterized)", nshape = 1L,
    shape_name = "log(nu)", shape_link = log, shape_inv = exp, shape_starts = c(1, 2),
    shape_line = function(theta) sprintf("\nnu (COM-Poisson dispersion; > 1 underdispersed): %.4f\n", theta),
    null_shape = 1,
    logpmf  = function(y, mu, theta) { ll <- cmp_loglambda_cpp(mu, theta)
      y * ll - theta * lgamma(y + 1) - cmp_logZ_cpp(ll, theta) },
    pvec    = function(mu, theta, kmax) { ll <- cmp_loglambda_cpp(mu, theta)
      exp(.mp_lterm(0:kmax, ll, theta) - cmp_logZ_cpp(ll, theta)) },
    p0      = function(mu, theta) exp(-cmp_logZ_cpp(cmp_loglambda_cpp(mu, theta), theta)),
    meanfun = function(mu, theta) mu,
    varfun  = function(mu, theta) cmp_moments_cpp(cmp_loglambda_cpp(mu, theta), theta)[, 2],
    kmax    = function(mu, theta) .mp_kmax(mu, theta)),
  genpois = list(tag = "genpois", label = "Generalized Poisson", nshape = 1L,
    shape_name = "atanh(lambda)", shape_link = atanh, shape_inv = tanh, shape_starts = c(0, -0.3),
    shape_line = function(theta) sprintf("\nlambda (generalized-Poisson dispersion; < 0 underdispersed): %.4f  [Var/Mean = %.3f]\n",
                                         theta, 1 / (1 - theta)^2),
    null_shape = 0,
    logpmf  = function(y, mu, theta) .gp_lterm(y, mu, theta) - .gp_logZ(mu, theta),
    pvec    = function(mu, theta, kmax) { K <- .gp_kmax(mu, theta); lt <- .gp_lterm(0:kmax, mu, theta)
      out <- exp(lt - .lse(.gp_lterm(0:K, mu, theta))); out[0:kmax > K] <- 0; out },
    p0      = function(mu, theta) exp(.gp_lterm(0, mu, theta) - .gp_logZ(mu, theta)),
    meanfun = function(mu, theta) .gp_moments(mu, theta)[1, ],
    varfun  = function(mu, theta) .gp_moments(mu, theta)[2, ],
    kmax    = function(mu, theta) .gp_kmax(mu, theta)),
  gammacount = list(tag = "gammacount", label = "Gamma-count", nshape = 1L,
    shape_name = "log(alpha)", shape_link = log, shape_inv = exp, shape_starts = c(1, 2),
    shape_line = function(theta) sprintf("\nalpha (gamma-count dispersion; > 1 underdispersed): %.4f  [Var/Mean ~ %.3f]\n",
                                         theta, 1 / theta),
    null_shape = 1,
    logpmf  = function(y, mu, theta) log(.gc_pmf(y, mu, theta)),
    pvec    = function(mu, theta, kmax) .gc_pmf(0:kmax, mu, theta),
    p0      = function(mu, theta) .gc_G(theta, theta * mu, lower = FALSE),
    meanfun = function(mu, theta) .gc_moments(mu, theta)[1, ],
    varfun  = function(mu, theta) .gc_moments(mu, theta)[2, ],
    kmax    = function(mu, theta) .gc_kmax(mu, theta)),
  doublepois = list(tag = "doublepois", label = "Double Poisson", nshape = 1L,
    shape_name = "log(theta)", shape_link = log, shape_inv = exp, shape_starts = c(1, 2),
    shape_line = function(theta) sprintf("\ntheta (double-Poisson dispersion; > 1 underdispersed): %.4f  [Var/Mean ~ %.3f]\n",
                                         theta, 1 / theta),
    null_shape = 1,
    logpmf  = function(y, mu, theta) .dp_lterm(y, mu, theta) - .dp_logZ(mu, theta),
    pvec    = function(mu, theta, kmax) exp(.dp_lterm(0:kmax, mu, theta) - .dp_logZ(mu, theta)),
    p0      = function(mu, theta) exp(.dp_lterm(0, mu, theta) - .dp_logZ(mu, theta)),
    meanfun = function(mu, theta) .dp_moments(mu, theta)[1, ],
    varfun  = function(mu, theta) .dp_moments(mu, theta)[2, ],
    kmax    = function(mu, theta) .dp_kmax(mu, theta))
)
.count_family_names <- c("poisson", "negbin", "compois", "mpcmp", "genpois", "gammacount", "doublepois")
.count_family_aliases <- c(compois_mean = "mpcmp", double_poisson = "doublepois", gamma_count = "gammacount",
                           generalized_poisson = "genpois", nb = "negbin", negative_binomial = "negbin")
.count_fam <- function(family) {
  if (length(family) > 1L) family <- family[1L]
  if (family %in% names(.count_family_aliases)) family <- .count_family_aliases[[family]]
  family <- match.arg(family, .count_family_names)
  .count_families[[family]]
}

## dispersion/shape display line for print/summary
.shape_line <- function(family, theta) .count_fam(family)$shape_line(theta)

## ---------------------------------------------------------------------------
## Distribution functions for the added families
## ---------------------------------------------------------------------------

## shared d/p/q/r engine from a family object: `mu` recycled against x/q/p
.fam_d <- function(fam, x, mu, theta, log) {
  if (!length(x) || !length(mu)) return(numeric(0))
  if (length(theta) != 1L) stop("the dispersion parameter must be a single value.")
  n <- max(length(x), length(mu)); x <- rep_len(x, n); mu <- rep_len(mu, n)
  ok <- x >= 0 & x == floor(x)
  ld <- rep(-Inf, n)
  if (any(ok)) ld[ok] <- fam$logpmf(x[ok], mu[ok], theta)
  if (log) ld else exp(ld)
}
.fam_p <- function(fam, q, mu, theta, lower.tail, log.p) {
  if (!length(q) || !length(mu)) return(numeric(0))
  if (length(theta) != 1L) stop("the dispersion parameter must be a single value.")
  n <- max(length(q), length(mu)); q <- rep_len(q, n); mu <- rep_len(mu, n)
  out <- vapply(seq_len(n), function(i) { qi <- floor(q[i]); if (qi < 0) return(0)
    sum(fam$pvec(mu[i], theta, as.integer(qi))) }, numeric(1))
  out <- pmin(pmax(out, 0), 1); if (!lower.tail) out <- 1 - out
  if (log.p) log(out) else out
}
.fam_q <- function(fam, p, mu, theta, lower.tail, log.p) {
  if (!length(p) || !length(mu)) return(numeric(0))
  if (length(theta) != 1L) stop("the dispersion parameter must be a single value.")
  if (log.p) p <- exp(p); if (!lower.tail) p <- 1 - p
  n <- max(length(p), length(mu)); p <- rep_len(p, n); mu <- rep_len(mu, n)
  vapply(seq_len(n), function(i) {
    K <- fam$kmax(mu[i], theta); cdf <- cumsum(fam$pvec(mu[i], theta, K))
    idx <- which(cdf >= p[i] - 1e-12)[1L]; as.numeric(if (is.na(idx)) K else idx - 1L) }, numeric(1))
}
.fam_r <- function(fam, n, mu, theta) {
  if (n == 0) return(numeric(0))
  if (length(theta) != 1L) stop("the dispersion parameter must be a single value.")
  mu <- rep_len(mu, n)
  vapply(seq_len(n), function(i) {
    K <- fam$kmax(mu[i], theta); pr <- fam$pvec(mu[i], theta, K); pr[!is.finite(pr) | pr < 0] <- 0
    as.numeric(sample.int(K + 1L, 1L, prob = pr) - 1L) }, numeric(1))
}

#' Gamma-count distribution functions
#'
#' Density, distribution, quantile, and random generation for Winkelmann's
#' (1995) gamma-count distribution: the number of events in a unit interval when
#' the waiting times between events are independent Gamma variables with shape
#' `alpha` and rate `alpha * mu`, so that `mu` is the long-run event rate.
#' `alpha > 1` (waiting times more regular than exponential) gives
#' underdispersion, `alpha = 1` is the Poisson, `alpha < 1` overdispersion; the
#' variance-to-mean ratio is approximately `1/alpha`. The exact mean and variance
#' are finite sums of incomplete-gamma terms and are what
#' `count_reg(family = "gammacount")` reports through `predict()` and
#' `fitted()`.
#'
#' @param x,q Vector of quantiles (non-negative integers).
#' @param p Vector of probabilities.
#' @param n Number of draws.
#' @param mu Rate parameter (scalar or vector, recycled).
#' @param alpha Dispersion parameter (scalar, positive).
#' @param log,log.p Return log probabilities.
#' @param lower.tail If `TRUE` (default), \eqn{P(X \le x)}.
#' @return `dgammacount` a density, `pgammacount` a CDF, `qgammacount` a
#'   quantile, `rgammacount` a numeric vector of count draws.
#' @references Winkelmann, R. (1995). Duration dependence and dispersion in
#'   count-data models. \emph{Journal of Business & Economic Statistics}, 13(4),
#'   467-474.
#' @examples
#' dgammacount(0:5, mu = 3, alpha = 2)
#' var(rgammacount(2000, mu = 3, alpha = 2)) / 3
#' @name gammacount-distribution
#' @seealso [count_reg()]
#' @export
dgammacount <- function(x, mu, alpha, log = FALSE) .fam_d(.count_families$gammacount, x, mu, alpha, log)
#' @rdname gammacount-distribution
#' @export
pgammacount <- function(q, mu, alpha, lower.tail = TRUE, log.p = FALSE)
  .fam_p(.count_families$gammacount, q, mu, alpha, lower.tail, log.p)
#' @rdname gammacount-distribution
#' @export
qgammacount <- function(p, mu, alpha, lower.tail = TRUE, log.p = FALSE)
  .fam_q(.count_families$gammacount, p, mu, alpha, lower.tail, log.p)
#' @rdname gammacount-distribution
#' @export
rgammacount <- function(n, mu, alpha) .fam_r(.count_families$gammacount, n, mu, alpha)

#' Double Poisson distribution functions
#'
#' Density, distribution, quantile, and random generation for Efron's (1986)
#' double Poisson distribution with location `mu` and dispersion `theta`,
#' normalized by the exact sum rather than Efron's closed-form approximation.
#' `theta > 1` gives underdispersion, `theta = 1` is the Poisson, `theta < 1`
#' overdispersion; the variance-to-mean ratio is approximately `1/theta`, and the
#' mean is close to (but not exactly) `mu`. The exact mean and variance are what
#' `count_reg(family = "doublepois")` reports through `predict()` and `fitted()`.
#'
#' @inheritParams gammacount-distribution
#' @param mu Location parameter (scalar or vector, recycled).
#' @param theta Dispersion parameter (scalar, positive).
#' @return `ddoublepois` a density, `pdoublepois` a CDF, `qdoublepois` a
#'   quantile, `rdoublepois` a numeric vector of count draws.
#' @references Efron, B. (1986). Double exponential families and their use in
#'   generalized linear regression. \emph{Journal of the American Statistical
#'   Association}, 81(395), 709-721.
#' @examples
#' ddoublepois(0:5, mu = 3, theta = 2)
#' sum(ddoublepois(0:60, mu = 3, theta = 2))
#' @name doublepois-distribution
#' @seealso [count_reg()]
#' @export
ddoublepois <- function(x, mu, theta, log = FALSE) .fam_d(.count_families$doublepois, x, mu, theta, log)
#' @rdname doublepois-distribution
#' @export
pdoublepois <- function(q, mu, theta, lower.tail = TRUE, log.p = FALSE)
  .fam_p(.count_families$doublepois, q, mu, theta, lower.tail, log.p)
#' @rdname doublepois-distribution
#' @export
qdoublepois <- function(p, mu, theta, lower.tail = TRUE, log.p = FALSE)
  .fam_q(.count_families$doublepois, p, mu, theta, lower.tail, log.p)
#' @rdname doublepois-distribution
#' @export
rdoublepois <- function(n, mu, theta) .fam_r(.count_families$doublepois, n, mu, theta)

#' Generalized Poisson distribution functions
#'
#' Density, distribution, quantile, and random generation for the Consul-Jain
#' generalized Poisson distribution in its constant-dispersion form: mean `mu`
#' and dispersion `lambda` in (-1, 1), with
#' \eqn{P(Y = y) = \theta(\theta + \lambda y)^{y-1} e^{-\theta - \lambda y}/y!}
#' and \eqn{\theta = \mu(1 - \lambda)}, so that the variance-to-mean ratio is
#' \eqn{1/(1-\lambda)^2}. `lambda < 0` gives underdispersion, `lambda = 0` the
#' Poisson, `lambda > 0` overdispersion. For `lambda < 0` the support is finite
#' (the terms are positive only while \eqn{\theta + \lambda y > 0}) and the pmf is
#' renormalized on that support, so it is a proper distribution and the exact
#' mean and variance (which then differ from `mu` and `mu/(1-lambda)^2` only by
#' the negligible mass beyond the support) are what
#' `count_reg(family = "genpois")` reports through `predict()` and `fitted()`.
#'
#' @inheritParams gammacount-distribution
#' @param mu Mean parameter (scalar or vector, recycled).
#' @param lambda Dispersion parameter in (-1, 1) (scalar).
#' @return `dgenpois` a density, `pgenpois` a CDF, `qgenpois` a quantile,
#'   `rgenpois` a numeric vector of count draws.
#' @references Consul, P. C. and Jain, G. C. (1973). A generalization of the
#'   Poisson distribution. \emph{Technometrics}, 15(4), 791-799. Consul, P. C.
#'   and Famoye, F. (2006). \emph{Lagrangian Probability Distributions}.
#'   Birkhauser.
#' @examples
#' dgenpois(0:5, mu = 3, lambda = -0.3)
#' var(rgenpois(2000, mu = 3, lambda = -0.3)) / 3   # about 1/(1.3)^2
#' @name genpois-distribution
#' @seealso [count_reg()]
#' @export
dgenpois <- function(x, mu, lambda, log = FALSE) .fam_d(.count_families$genpois, x, mu, lambda, log)
#' @rdname genpois-distribution
#' @export
pgenpois <- function(q, mu, lambda, lower.tail = TRUE, log.p = FALSE)
  .fam_p(.count_families$genpois, q, mu, lambda, lower.tail, log.p)
#' @rdname genpois-distribution
#' @export
qgenpois <- function(p, mu, lambda, lower.tail = TRUE, log.p = FALSE)
  .fam_q(.count_families$genpois, p, mu, lambda, lower.tail, log.p)
#' @rdname genpois-distribution
#' @export
rgenpois <- function(n, mu, lambda) .fam_r(.count_families$genpois, n, mu, lambda)
