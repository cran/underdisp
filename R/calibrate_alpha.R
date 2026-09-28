## ---------------------------------------------------------------------------
## Opt-in calibration of the interval for the CPB dispersion parameter.
##
## The support of the CPB, 0..floor(lambda_i/(1-alpha)), moves with the
## parameters, and it does so only through the log-ceiling coefficients (the
## intercept less log(1-alpha), and the slopes). Those are endpoint-type
## parameters: the likelihood gains each time an unobserved top support point is
## shed and stops at the first count on its ceiling. Given the ceilings, alpha
## is a regular exponential-family parameter, but its estimate inherits a bias
## toward zero from them, so the signed root of the profile likelihood ratio,
##   r(alpha) = sign(alpha-hat - alpha) sqrt(2 (l-hat - l_p(alpha))),
## is not centred at zero and the chi-square cut gives an interval that covers
## about 0.90 to 0.94, its misses on the upper side. The calibration replaces
## the standard normal quantiles of r by its quantiles under the fitted model
## (a parametric bootstrap): { alpha : q_lo <= r(alpha) <= q_hi }.
## ---------------------------------------------------------------------------

#' Calibrate the interval for the CPB dispersion parameter
#'
#' An opt-in step that stores, in the fit, the parametric-bootstrap
#' distribution of the signed root of the profile likelihood ratio for `alpha`.
#' Once it is stored, [alpha_confint()], [confint.cpb()], [implied_ceiling()]
#' and `summary()` report the calibrated interval instead of the first-order
#' one.
#'
#' @details
#' The first-order interval cuts the profile log-likelihood `qchisq(level, 1) / 2`
#' below its maximum on both sides. For the CPB that cut is not calibrated: the
#' support `0, ..., floor(lambda / (1 - alpha))` moves with the parameters, the
#' log-ceiling coefficients behave like endpoint parameters, and `alpha`-hat
#' inherits a bias toward zero from them. In simulations from the CPB the
#' first-order 95% interval covers between 0.90 and 0.94 across a grid of nine
#' designs, and nearly all of its misses have the true `alpha` above the upper
#' limit.
#'
#' The calibration simulates `B` responses from the fitted model (the fitted
#' rates, `alpha`-hat, the same truncation and offset), refits each one, and
#' records the signed root `r* = sign(alpha* - alpha-hat) sqrt(2 (l* - l*_p(alpha-hat)))`.
#' The interval is then `{alpha : q_lo <= r(alpha) <= q_hi}`, with `q_lo` and
#' `q_hi` the order statistics of the `r*` at ranks `(B + 1)(1 - level)/2` from
#' each end (the 5th smallest and 5th largest of 199 at the 95% level), so the
#' upper limit sits `q_lo^2 / 2` and the lower limit `q_hi^2 / 2` below the
#' maximum of the profile. In the same simulations the calibrated 95% interval
#' covers 0.92 to 0.96, with its misses balanced between the two sides. The
#' first-order interval falls short only when the design has many distinct
#' covariate patterns (the endpoint effect needs many distinct ceilings): with an
#' intercept only, or a covariate taking ten or fewer values, it covers at or
#' above its level, and the calibration, which stays close to the nominal level
#' there too (0.95 with a binary covariate), is not needed.
#'
#' The responses are all drawn first, on the calling process, and the refits
#' use no random numbers, so the result depends on `seed` but not on `cores`.
#' With `B = 199` the two cuts carry simulation error (a standard deviation of
#' about 0.2 in `q_lo` and `q_hi`), which moves a limit by a fraction of its
#' distance from the estimate; `B = 999` reduces it and is needed for a 99%
#' interval. Each replicate costs about a third of the original fit.
#'
#' The interval is model-based. It assumes the CPB and independent
#' observations, and it is not cluster-robust: when the fit was given
#' `cluster`, a message says so. Under misspecification `alpha` has no fixed
#' target (on data with a soft upper tail its estimate rises with the sample
#' size), and no interval computed under the fitted CPB repairs that.
#'
#' Calibration is refused, with the reason, where it has nothing to work with: a
#' fit whose ceiling reaches `max.support` (the guard, not the data, bounds
#' `alpha`), weights that are not whole numbers, and a fit saved by version
#' 0.1.0. Whole-number weights are treated as frequency weights: a row of weight
#' `w` contributes `w` independent simulated responses.
#'
#' @param object A pooled [cpb()] fit, or a [hurdle_cpb()] fit whose intensity
#'   is a pooled CPB.
#' @param B Number of bootstrap replicates (default 199).
#' @param cores Worker processes for the refits.
#' @param seed Optional seed for the simulated responses; the caller's random
#'   number state is restored afterwards, as in [stats::simulate()].
#' @param ... Unused.
#' @return `object`, with the component `alpha.cal`: a list holding the signed
#'   roots `r` (`NA` for a replicate that failed to fit), `B`, `B_ok` and `seed`.
#' @examples
#' \donttest{
#' set.seed(7); x <- rnorm(100)
#' d <- data.frame(y = rcpb(100, exp(1.5 + 0.4 * x), 0.5, truncated = TRUE), x = x)
#' fit <- cpb(y ~ x, d, se = "none")
#' alpha_confint(fit)                       # first-order
#' ## 19 replicates keep the example short and can serve a 90% interval; use 199 or more
#' fit <- calibrate_alpha(fit, B = 19, seed = 1)
#' alpha_confint(fit, level = 0.90)         # calibrated
#' }
#' @seealso [alpha_confint()], [dispersion_test()]
#' @export
calibrate_alpha <- function(object, B = 199, cores = 1L, seed = NULL, ...) UseMethod("calibrate_alpha")

#' @rdname calibrate_alpha
#' @method calibrate_alpha cpb
#' @export
calibrate_alpha.cpb <- function(object, B = 199, cores = 1L, seed = NULL, ...) {
  .ud_check_whole(B, "B", lower = 19)
  if (inherits(object, "cpb_fe"))
    stop("calibrate_alpha() applies to the pooled cpb() fit; a fixed-effects fit has no profile interval for alpha.", call. = FALSE)
  if (is.null(object$profile) || is.null(object$rate) || is.null(object$X) || is.null(object$Y))
    stop("This fit was saved by an earlier version of the package and does not carry what the calibration needs; ",
         "refit with the current version.", call. = FALSE)
  if (isTRUE(object$support_binding))
    stop("The fitted ceiling reaches max.support, so alpha is bounded by the guard and not by the data; ",
         "there is no interval to calibrate.", call. = FALSE)
  ex <- .disp_expand(object)
  if (is.null(ex))
    stop("The weights are not whole numbers, so the fit has no replicate data to simulate; the first-order interval stays in use.",
         call. = FALSE)
  ob <- ex$object; X <- ob$X; off <- .disp_offset(ob); n <- nrow(X); p <- ncol(X)
  if (isTRUE(object$clustered))
    message("The calibrated interval for alpha is model-based: it assumes independent observations and does not use 'cluster'.")
  ms <- object$max.support; tr <- isTRUE(object$truncated)
  ahat <- object$alpha; bhat <- unname(object$coefficients); rate <- object$rate[ex$index]

  ## all B responses first, on this process: the refits below use no random numbers
  ystar <- .sim_draws(seed, function() .sim_from_pmf(n, B, function(i) .sim_pmf_cpb(rate[i], ahat, ms), truncated = tr))

  s <- .ud_colscale(X, rep(1, n)); Xs <- sweep(X, 2, s, "/")
  icol <- match("(Intercept)", colnames(X))
  astarts <- unique(c(0.5, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95))
  one <- function(b) {
    ys <- ystar[, b]
    tryCatch(suppressWarnings({
      g <- .cpb_fit_ms(X, ys, off, ms, tr, astarts, 20000, 1e-8, w = NULL)
      if (is.null(g) || !is.finite(g$value) || g$value >= 1e9) return(NA_real_)
      g <- .cpb_alpha_scan(X, ys, off, ms, tr, NULL, g, 20000, 1e-8)
      bs <- tryCatch({ v <- stats::glm.fit(Xs, ys, offset = off, family = stats::poisson())$coefficients; v[!is.finite(v)] <- 0; v },
                     error = function(e) rep(0, p))
      lp <- .cpb_prof_at(ahat, bhat * s, Xs, ys, off, ms, tr, NULL, s, icol, ys > 0, b_cold = bs)$ll
      top <- max(-g$value, lp)                                  # a restricted maximum above the refit's: the refit was short
      sign(stats::plogis(g$par[p + 1L]) - ahat) * sqrt(2 * (top - lp))
    }), error = function(e) NA_real_)
  }
  r <- unlist(.ud_lapply(seq_len(B), one, cores = cores))
  B_ok <- sum(is.finite(r))
  if (B_ok < 0.8 * B)
    warning("calibrate_alpha: ", B - B_ok, " of ", B, " bootstrap replicates failed to fit; they are counted against both tails.",
            call. = FALSE)
  object$alpha.cal <- list(r = r, B = as.integer(B), B_ok = B_ok, seed = seed)
  object
}

#' @rdname calibrate_alpha
#' @method calibrate_alpha hurdle_cpb
#' @export
calibrate_alpha.hurdle_cpb <- function(object, B = 199, cores = 1L, seed = NULL, ...) {
  if (!inherits(object$intensity, "cpb") || inherits(object$intensity, "cpb_fe"))
    stop("calibrate_alpha() needs a pooled CPB intensity; a fixed-effects intensity has no profile interval for alpha.", call. = FALSE)
  object$intensity <- calibrate_alpha(object$intensity, B = B, cores = cores, seed = seed)
  object
}

#' @rdname calibrate_alpha
#' @method calibrate_alpha default
#' @export
calibrate_alpha.default <- function(object, B = 199, cores = 1L, seed = NULL, ...)
  stop("calibrate_alpha() applies to cpb() fits and to hurdle_cpb() fits with a pooled intensity.", call. = FALSE)

## draw with an optional seed, restoring the caller's RNG state when one is given
.sim_draws <- function(seed, draw) {
  if (is.null(seed)) return(draw())
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
  prev <- get(".Random.seed", envir = globalenv())
  on.exit(assign(".Random.seed", prev, envir = globalenv()))
  set.seed(seed)
  draw()
}

## The two cuts of a calibrated interval, as depths below the maximum of the
## profile: order statistics of the stored signed roots at rank
## k = floor((B + 1)(1 - level)/2) from each end, with every failed replicate
## counted beyond BOTH tails (a fit that fails is more likely a tail fit).
.alpha_cal_drops <- function(cal, level) {
  a <- (1 - level) / 2; r <- sort(cal$r[is.finite(cal$r)]); fails <- cal$B - length(r)
  k <- floor((cal$B + 1) * a + 1e-9) - fails
  if (k < 1L)
    stop(sprintf(paste0("%d replicates (%d failed) cannot serve a %s%% interval: it needs at least %d successful ones; ",
                        "rerun calibrate_alpha() with a larger B."), cal$B, fails, format(100 * level), ceiling(1 / a) - 1L + fails),
         call. = FALSE)
  q_lo <- r[k]; q_hi <- r[length(r) - k + 1L]
  c(lower = max(q_hi, 0)^2 / 2, upper = max(-q_lo, 0)^2 / 2)
}

## The interval every accessor reports: calibrated when the fit carries the
## bootstrap signed roots, first-order otherwise; attr "method" says which.
.cpb_alpha_ci <- function(object, level = 0.95) {
  cal <- object$alpha.cal
  if (is.null(cal)) {
    ci <- .cpb_alpha_profile_ci(object, level = level)
    attr(ci, "method") <- "first-order profile likelihood"
    return(ci)
  }
  ci <- .cpb_alpha_profile_ci(object, level = level, drops = .alpha_cal_drops(cal, level))
  attr(ci, "method") <- sprintf("profile likelihood, signed root calibrated by parametric bootstrap (%d of %d replicates)",
                                as.integer(cal$B_ok), as.integer(cal$B))
  ci
}
