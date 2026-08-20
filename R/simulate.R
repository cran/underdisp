## ---------------------------------------------------------------------------
## Simulation methods and DHARMa interoperability.
##
## simulate() draws replicate response vectors from every fitted model in the
## family, and fitted() is completed for the two-part and fixed-effects classes,
## so any fit can be handed to DHARMa::createDHARMa() for simulated-residual
## diagnostics:
##
##   sims <- simulate(fit, nsim = 250)
##   DHARMa::createDHARMa(simulatedResponse = as.matrix(sims),
##                        observedResponse  = fit$y,     # $Y for single-equation fits
##                        fittedPredictedResponse = fitted(fit), integerResponse = TRUE)
##
## One engine serves every class: for each observation the model's pmf is
## computed once on its support and nsim counts are drawn from it, so runtime is
## n pmf evaluations regardless of nsim. Zero-truncation is applied by one rule
## (.sim_norm): zero out the origin, renormalize, and if a degenerate support
## survives (a zero-truncated CPB whose ceiling is below 1) put the mass at the
## smallest admissible count, matching rcpb(). gec_fe inherits simulate()/
## fitted() through its c("gec_fe", "gec") class vector.
## ---------------------------------------------------------------------------

## normalize a raw pmf vector, optionally zero-truncating: one rule for every family
.sim_norm <- function(pr, truncated) {
  pr[!is.finite(pr) | pr < 0] <- 0
  if (truncated && length(pr)) pr[1] <- 0
  s <- sum(pr)
  if (s <= 0) { pr <- numeric(max(2L, length(pr))); pr[if (truncated) 2L else 1L] <- 1 }
  else pr <- pr / s
  pr
}

## draw an n x nsim integer matrix; pmf_fun(i) returns observation i's raw pmf over 0..K_i
.sim_from_pmf <- function(n, nsim, pmf_fun, truncated) {
  out <- matrix(0L, n, nsim)
  for (i in seq_len(n)) {
    pr <- .sim_norm(pmf_fun(i), truncated)
    out[i, ] <- as.integer(sample.int(length(pr), nsim, replace = TRUE, prob = pr) - 1L)
  }
  out
}

## raw (untruncated) pmf rows per family
.sim_pmf_cpb <- function(lambda, alpha, max.support) {
  .cpb_pmf_row(lambda, alpha, truncated = FALSE, max.support = max.support)
}
.sim_pmf_gec <- function(rate, delta, max.support) {
  km <- .gec_kmax_qr(rate, delta, max.support)
  as.numeric(gec_pmf_cpp(rate, delta, km, as.integer(max.support)))
}
.sim_pmf_count <- function(fam, mu, theta) {
  kmax <- switch(fam$tag,
    poisson = stats::qpois(1 - 1e-10, mu),
    negbin  = stats::qnbinom(1 - 1e-10, mu = mu, size = theta),
    compois = .compois_kmax(mu, theta))
  fam$pvec(mu, theta, max(1L, as.integer(kmax)))
}

## stats::simulate() conventions: nsim columns named sim_1.., a "seed" attribute,
## and the caller's RNG state restored when an explicit seed is supplied
.sim_wrap <- function(nsim, seed, draw) {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = globalenv())
  } else {
    prev <- get(".Random.seed", envir = globalenv())
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", prev, envir = globalenv()))
  }
  out <- as.data.frame(draw())
  names(out) <- paste0("sim_", seq_len(nsim))
  attr(out, "seed") <- rng_state
  out
}

## binary stage of the two-part draws: n x nsim indicator matrix with row-wise
## probabilities (participation for hurdles, structural zero for mixtures)
.sim_binary <- function(prob, n, nsim) {
  matrix(stats::rbinom(n * nsim, 1L, rep_len(prob, n * nsim)), n, nsim)
}

#' Simulate responses from a fitted underdisp model
#'
#' Draws `nsim` replicate response vectors from the fitted model, at the
#' estimated parameters and the estimation data, in the format of
#' [stats::simulate()]. Available for every model class in the package:
#' `cpb`, `cpb_fe`, `gec`, `gec_fe` (through inheritance), `count_reg`,
#' `hurdle_count`, `zi_count`, `hurdle_cpb`, `hurdle_gec`, `zi_cpb`, and
#' `zi_gec`. Two-part models first draw the binary stage (participation or
#' structural zero), then the count stage from its own distribution, so the
#' replicates carry the model's full zero structure.
#'
#' The main consumer is simulated-residual diagnostics:
#' `DHARMa::createDHARMa(simulatedResponse = as.matrix(simulate(fit, 250)),
#' observedResponse = y, fittedPredictedResponse = fitted(fit),
#' integerResponse = TRUE)` works for any fit in the family.
#'
#' @param object A fitted model from this package.
#' @param nsim Number of replicate response vectors.
#' @param seed Optional seed, handled as in [stats::simulate()]: the caller's
#'   RNG state is restored on exit when a seed is supplied.
#' @param ... Unused.
#' @return A data frame with `nsim` integer columns, one row per observation,
#'   with a `"seed"` attribute.
#' @examples
#' set.seed(1); x <- rnorm(200)
#' N <- pmax(round(exp(1.4 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
#' fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "none")
#' sims <- simulate(fit, nsim = 5)
#' colMeans(sims)
#' @name simulate.underdisp
NULL

#' @rdname simulate.underdisp
#' @method simulate cpb
#' @export
simulate.cpb <- function(object, nsim = 1, seed = NULL, ...) {
  lam <- object$fitted.values; a <- object$alpha
  ms  <- if (is.null(object$max.support)) 500L else object$max.support
  .sim_wrap(nsim, seed, function()
    .sim_from_pmf(object$n, nsim, function(i) .sim_pmf_cpb(lam[i], a, ms),
                  truncated = isTRUE(object$truncated)))
}

#' @rdname simulate.underdisp
#' @method simulate cpb_fe
#' @export
simulate.cpb_fe <- function(object, nsim = 1, seed = NULL, ...) {
  lam <- object$fitted.values; a <- object$alpha
  ms  <- if (is.null(object$max.support)) 500L else object$max.support
  .sim_wrap(nsim, seed, function()
    .sim_from_pmf(object$n, nsim, function(i) .sim_pmf_cpb(lam[i], a, ms),
                  truncated = isTRUE(object$truncated)))
}

#' @rdname simulate.underdisp
#' @method simulate gec
#' @export
simulate.gec <- function(object, nsim = 1, seed = NULL, ...) {
  rate <- exp(object$linear.predictors); d <- object$delta; ms <- object$max.support
  .sim_wrap(nsim, seed, function()
    .sim_from_pmf(object$n, nsim, function(i) .sim_pmf_gec(rate[i], d, ms),
                  truncated = isTRUE(object$truncated)))
}

#' @rdname simulate.underdisp
#' @method simulate count_reg
#' @export
simulate.count_reg <- function(object, nsim = 1, seed = NULL, ...) {
  fam <- .count_fam(object$family); mu <- object$fitted.values; th <- object$theta
  .sim_wrap(nsim, seed, function()
    .sim_from_pmf(object$n, nsim, function(i) .sim_pmf_count(fam, mu[i], th),
                  truncated = isTRUE(object$truncated)))
}

#' @rdname simulate.underdisp
#' @method simulate hurdle_cpb
#' @export
simulate.hurdle_cpb <- function(object, nsim = 1, seed = NULL, ...) {
  n <- object$n; lam <- object$lambda_full; a <- object$intensity$alpha
  ms <- if (is.null(object$intensity$max.support)) 500L else object$intensity$max.support
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_cpb(lam[i], a, ms), truncated = TRUE)
    cnt * .sim_binary(object$p_full, n, nsim)
  })
}

#' @rdname simulate.underdisp
#' @method simulate hurdle_gec
#' @export
simulate.hurdle_gec <- function(object, nsim = 1, seed = NULL, ...) {
  n <- object$n; lam <- object$lambda_full; d <- object$delta; ms <- object$max.support
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_gec(lam[i], d, ms), truncated = TRUE)
    cnt * .sim_binary(object$p_full, n, nsim)
  })
}

#' @rdname simulate.underdisp
#' @method simulate hurdle_count
#' @export
simulate.hurdle_count <- function(object, nsim = 1, seed = NULL, ...) {
  fam <- .count_fam(object$family)
  n <- object$n; lam <- object$lambda_full; th <- object$theta
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_count(fam, lam[i], th), truncated = TRUE)
    cnt * .sim_binary(object$p_full, n, nsim)
  })
}

#' @rdname simulate.underdisp
#' @method simulate zi_cpb
#' @export
simulate.zi_cpb <- function(object, nsim = 1, seed = NULL, ...) {
  n <- object$n; lam <- object$lambda_full; a <- object$alpha
  ms <- if (is.null(object$max.support)) 500L else object$max.support
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_cpb(lam[i], a, ms), truncated = FALSE)
    cnt * (1L - .sim_binary(object$pi_full, n, nsim))
  })
}

#' @rdname simulate.underdisp
#' @method simulate zi_gec
#' @export
simulate.zi_gec <- function(object, nsim = 1, seed = NULL, ...) {
  n <- object$n; lam <- object$lambda_full; d <- object$delta; ms <- object$max.support
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_gec(lam[i], d, ms), truncated = FALSE)
    cnt * (1L - .sim_binary(object$pi_full, n, nsim))
  })
}

#' @rdname simulate.underdisp
#' @method simulate zi_count
#' @export
simulate.zi_count <- function(object, nsim = 1, seed = NULL, ...) {
  fam <- .count_fam(object$family)
  n <- object$n; lam <- object$lambda_full; th <- object$theta
  .sim_wrap(nsim, seed, function() {
    cnt <- .sim_from_pmf(n, nsim, function(i) .sim_pmf_count(fam, lam[i], th), truncated = FALSE)
    cnt * (1L - .sim_binary(object$pi_full, n, nsim))
  })
}

#' Fitted means for the two-part and fixed-effects classes
#'
#' Completes [stats::fitted()] across the family, returning the same quantity as
#' `predict(object, type = "response")` on the estimation data: the mean for the
#' single-equation fits and the marginal expected count (zeros included) for the
#' two-part fits. `gec_fe` inherits `fitted.gec`.
#'
#' @param object A fitted model from this package.
#' @param ... Unused.
#' @return A numeric vector, one element per estimation observation.
#' @name fitted.underdisp
NULL

#' @rdname fitted.underdisp
#' @method fitted cpb_fe
#' @export
fitted.cpb_fe <- function(object, ...) object$fitted.values

#' @rdname fitted.underdisp
#' @method fitted hurdle_cpb
#' @export
fitted.hurdle_cpb <- function(object, ...) object$p_full * object$lambda_full

#' @rdname fitted.underdisp
#' @method fitted hurdle_gec
#' @export
fitted.hurdle_gec <- function(object, ...) {
  p0 <- gec_pmf_cpp(object$lambda_full, object$delta, 0L, object$max.support)[, 1]
  object$p_full * object$lambda_full / pmax(1 - p0, 1e-8)
}

#' @rdname fitted.underdisp
#' @method fitted hurdle_count
#' @export
fitted.hurdle_count <- function(object, ...) {
  fam <- .count_fam(object$family)
  object$p_full * fam$meanfun(object$lambda_full, object$theta) /
    pmax(1 - fam$p0(object$lambda_full, object$theta), 1e-8)
}

#' @rdname fitted.underdisp
#' @method fitted zi_cpb
#' @export
fitted.zi_cpb <- function(object, ...) (1 - object$pi_full) * object$lambda_full

#' @rdname fitted.underdisp
#' @method fitted zi_gec
#' @export
fitted.zi_gec <- function(object, ...)
  (1 - object$pi_full) * gec_mean_cpp(object$lambda_full, object$delta, object$max.support)

#' @rdname fitted.underdisp
#' @method fitted zi_count
#' @export
fitted.zi_count <- function(object, ...) {
  fam <- .count_fam(object$family)
  (1 - object$pi_full) * fam$meanfun(object$lambda_full, object$theta)
}
