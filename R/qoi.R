## ---------------------------------------------------------------------------
## Quantities of interest for CPB fits
## ---------------------------------------------------------------------------

## CPB probability mass for one (lambda, alpha): P(Y = k), k = 0..K. If truncated,
## P(Y = 0) = 0 and the remaining mass is renormalized to sum to one. On the
## floor support the support stops at max.support (a larger ceiling is
## infeasible there); on King's support max.support bounds only the width of the
## likelihood's pmf body, so the shared core's own end applies (NA past its cap).
.cpb_pmf_row <- function(lambda, alpha, truncated, max.support = 500, support) {
  pr <- if (.cpb_fl(support)) .cpb_pmf_core(lambda, alpha, support, kcap = as.integer(max.support))
        else .cpb_pmf_core(lambda, alpha, support)
  if (anyNA(pr)) return(NA_real_)
  if (truncated) { p0 <- pr[1]; pr <- pr / (1 - p0); pr[1] <- 0 }
  pr
}

## P(Y = y) at given (lambda, alpha), vectorized over equal-length lambda and y.
## A pmf refused at the 1e5-term support cap (a rate of order 1e5 on King's support) gives NA, with dcpb()'s warning.
.cpb_prob_at <- function(lambda, alpha, y, truncated, max.support = 500, support) {
  out <- mapply(function(l, yy) {
    if (yy < 0 || yy != floor(yy)) return(0)
    pr <- .cpb_pmf_row(l, alpha, truncated, max.support, support = support)
    if (length(pr) == 1L && is.na(pr)) NA_real_ else if (yy + 1 > length(pr)) 0 else pr[yy + 1]
  }, lambda, y)
  cap <- is.na(out) & is.finite(rep_len(lambda, length(out)))
  if (any(cap)) warning("the CPB support reaches the 1e5-term cap for ", sum(cap), " value(s) (a rate of ",
                        "order 1e5): NA returned.", call. = FALSE)
  out
}

#' Predictions from a CPB fit
#'
#' @param object A `"cpb"` object.
#' @param newdata Optional data frame of new covariate profiles; if omitted, the
#'   fitted data are used.
#' @param type One of `"response"` (the exact mean of the fitted distribution;
#'   for a zero-truncated fit the conditional mean E(Y | Y >= 1)), `"rate"`
#'   (the CPB rate parameter lambda = exp(x'b), which equals the mean only when
#'   lambda/(1-alpha) is an integer), `"link"` (the linear predictor, log
#'   lambda), `"ceiling"` (the implied ceiling lambda/(1-alpha), the CPB's
#'   continuous number of trials; the largest count with positive probability is
#'   its ceiling on King's support and its floor on the floor support), or
#'   `"prob"` (the probability that `Y` equals `at`).
#' @param at For `type = "prob"`, the count value(s) `y` whose probability is
#'   returned (length 1, or one per row of the prediction data).
#' @param offset Optional offset (log scale) for the `newdata` rows: a numeric
#'   vector or a column name in `newdata`. By default a fit made with an offset
#'   column reads that column from `newdata`, and one made with an offset vector
#'   requires this argument. For `newdata = NULL` the fit's own offset is used.
#' @param se.fit,interval,level The uncertainty of the predicted quantity at
#'   the rows of `newdata`: `se.fit = TRUE` adds standard errors, and
#'   `interval = "confidence"` a `level` confidence interval. They follow the
#'   rule of [first_difference()]: the fit's stored bootstrap draws (a
#'   percentile interval, and the draws' standard deviation), otherwise the
#'   delta method over the fit's covariance, otherwise none. They need
#'   `newdata`. A row at a fixed-effects unit's own effect, or at its own
#'   zero-equation dummy, has none, since the draws do not carry that
#'   parameter; a row of a unit the fit has not seen, or `newdata` without the
#'   unit column, enters at the reference unit, which they carry. [ame()] gives
#'   average marginal effects.
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @examples
#' set.seed(1); x <- rnorm(300)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
#' fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "none")
#' predict(fit, newdata = data.frame(x = 0), type = "ceiling")
#' predict(fit, newdata = data.frame(x = c(-1, 1)), type = "prob", at = 5)
#' @method predict cpb
#' @export
predict.cpb <- function(object, newdata = NULL,
                        type = c("response", "rate", "link", "ceiling", "prob"), at = NULL,
                        offset = NULL, se.fit = FALSE, interval = c("none", "confidence"),
                        level = 0.95, ...) {
  type <- match.arg(type)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset, at = at))
  if (is.null(newdata)) {
    X <- object$X
    off <- if (!is.null(object$offset)) object$offset else rep_len(0, nrow(X))
  } else {
    Terms <- delete.response(object$terms)
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss))
      stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "),
           ". Supply every predictor in the model.")
    X  <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    off <- .ud_newdata_offset(object, newdata, offset, nrow(X))
  }
  bad <- if (is.null(newdata)) FALSE else .ud_na_rows(X, off)     # a missing covariate or offset predicts NA
  if (any(bad)) { X[bad, ] <- 0; off[bad] <- 0 }
  eta <- as.numeric(off + X %*% object$coefficients); lam <- exp(eta)
  out <- switch(type,
    link     = eta,
    rate     = lam,
    response = .cpb_mean(lam, object$alpha, isTRUE(object$truncated), .cpb_support(object)),
    ceiling  = lam / (1 - object$alpha),
    prob = {
      if (is.null(at)) stop("For type = \"prob\", supply 'at' (the count value y).")
      yv <- if (length(at) == 1) rep(at, length(lam)) else at
      if (length(yv) != length(lam)) stop("'at' must be length 1 or nrow(newdata).")
      .cpb_prob_at(lam, object$alpha, yv, object$truncated, object$max.support, support = .cpb_support(object))
    })
  .ud_mask(out, bad)
}

#' Implied ceiling, and the ceiling at the limits of alpha's interval
#'
#' Returns the observation- (or profile-) specific ceiling lambda/(1-alpha), and
#' the ceiling at the fitted rate with `alpha` at the limits of its interval
#' (see [alpha_confint()]; calibrated when the fit went through
#' [calibrate_alpha()]). Those two columns hold the rate at its estimate, so
#' they show how far the dispersion's uncertainty alone moves the ceiling and
#' are not a confidence interval for it; [first_difference()] with
#' `quantity = "ceiling"` propagates the coefficients' uncertainty too. The ceiling is the CPB's continuous number of trials. No count above
#' its integer ceiling has positive probability on King's support (the default),
#' and none above its integer floor on the floor support (`support = "floor"`);
#' that largest feasible count is returned as `max_count`. Whether held-out
#' counts respect the bound is tested by [ceiling_exceedance()].
#'
#' The ceiling is the fitted distribution's support end. With the rate and
#' `alpha` it is fixed by the CPB's variance-mean link, so it is a bound on the
#' counts only when the CPB is the data's process. On counts with a common
#' hard bound from another process (binomial counts out of six trials, say),
#' the fitted ceiling follows the mean and the dispersion instead: in one
#' simulation it ran from 3 where the mean was low to 12 where it was high,
#' below the true bound in the first case and well above it in the second.
#'
#' @param object A `"cpb"` object.
#' @param newdata Optional covariate profiles.
#' @param level Confidence level (default 0.95).
#' @return A data frame with `lambda`, `ceiling`, `max_count` (the largest count
#'   with positive probability), and `ceiling_at_alpha_lower` /
#'   `ceiling_at_alpha_upper`, the ceiling at the fitted rate with `alpha` at the
#'   lower and upper limits of its interval (for `cpb` fits, and hurdles with a
#'   `cpb` intensity).
#' @param ... Further arguments passed to methods.
#' @examples
#' set.seed(6); x <- rnorm(300)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
#' implied_ceiling(fit, newdata = data.frame(x = c(-1, 0, 1)))
#' @export
implied_ceiling <- function(object, ...) UseMethod("implied_ceiling")

#' @rdname implied_ceiling
#' @method implied_ceiling cpb
#' @export
implied_ceiling.cpb <- function(object, newdata = NULL, level = 0.95, ...) {
  lam <- predict(object, newdata = newdata, type = "rate")
  aci <- .cpb_alpha_ci(object, level = level)
  data.frame(lambda    = lam,
             ceiling   = lam / (1 - object$alpha),
             max_count = .cpb_top(lam / (1 - object$alpha), .cpb_support(object)),
             ceiling_at_alpha_lower = lam / (1 - aci[["lower"]]),
             ceiling_at_alpha_upper = lam / (1 - aci[["upper"]]),
             row.names = NULL)
}

#' @rdname implied_ceiling
#' @method implied_ceiling cpb_fe
#' @export
implied_ceiling.cpb_fe <- function(object, newdata = NULL, level = 0.95, ...) {
  ## the concentrated fit has no profile interval for alpha (the profile would
  ## re-concentrate every unit effect at each alpha); the ceiling is reported
  ## at the point estimate, at the average unit for newdata
  lam <- predict(object, newdata = newdata, type = "rate")
  data.frame(lambda = lam, ceiling = lam / (1 - object$alpha),
             max_count = .cpb_top(lam / (1 - object$alpha), .cpb_support(object)), row.names = NULL)
}

#' Held-out ceiling exceedance of a CPB fit
#'
#' Whether counts respect the CPB's ceiling can be checked, but not in sample:
#' its likelihood is minus infinity for any count above the bound, so every
#' feasible fit respects the bound on the data it was fit to. The check is
#' one-sided: it detects a ceiling that is too low, where held-out counts reach
#' it, and has no power against one that is too high (a fitted ceiling above a
#' true bound is never exceeded), nor where counts stay well below the ceiling
#' (see [implied_ceiling()]). `ceiling_exceedance()` takes held-out
#' observations and reports how many of their counts exceed the largest count
#' the fitted model allows at their covariates (`max_count` of
#' [implied_ceiling()]); for a two-part model the bound is the count
#' component's, and a fixed-effects unit seen in training enters at its own
#' effect (an unseen one at the mean effect), as in [score()]. A breached ceiling also zeroes the held-out likelihood, which
#' [score()] records with its 1e-12 floor.
#'
#' @param fit A `"cpb"`, `"cpb_fe"`, `"hurdle_cpb"`, or `"zi_cpb"` object.
#' @param newdata Held-out observations, including the response.
#' @param weights Optional frequency weights for the rows of `newdata`: a column
#'   name of `newdata` or a numeric vector with one value per row.
#' @param offset Optional offset (log scale) of the count for the `newdata`
#'   rows; by default a fit made with an offset reads it from the same column of
#'   `newdata` (see [predict.cpb()]).
#' @return A one-row data frame with `n` (held-out rows, weighted), `exceed`
#'   (rows above the bound, weighted), and `share`.
#' @seealso [implied_ceiling()], [score()]
#' @examples
#' set.seed(7); x <- rnorm(400); y <- rcpb(400, exp(1.4 + 0.3 * x), 0.5)
#' d <- data.frame(y = y, x = x)
#' fit <- cpb(y ~ x, d[1:300, ], truncated = FALSE, se = "none")
#' ceiling_exceedance(fit, d[301:400, ])
#' @export
ceiling_exceedance <- function(fit, newdata, weights = NULL, offset = NULL) {
  if (!inherits(fit, c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb")))
    stop("ceiling_exceedance() applies to the CPB family (cpb, cpb_fe, hurdle_cpb, zi_cpb).")
  if (missing(newdata) || !is.data.frame(newdata)) stop("'newdata' must be a data frame of held-out observations.")
  y <- stats::model.response(stats::model.frame(stats::update(fit$formula, . ~ 1), newdata, na.action = NULL))
  ## held-out rates as score() computes them: a unit seen in training at its own fixed effect, else the mean
  if (inherits(fit, "cpb")) {
    lam <- as.numeric(exp(.ud_newdata_offset(fit, newdata, offset) +
                          .ud_newdata_matrix(stats::delete.response(fit$terms), newdata, fit$levels, fit$contrasts,
                                             names(fit$coefficients)) %*% fit$coefficients)); a <- fit$alpha
  } else if (inherits(fit, "cpb_fe")) {
    X <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    fe <- fit$fe[as.character(newdata[[fit$fe_var]])]; fe[is.na(fe)] <- .fe_ref(fit)
    lam <- as.numeric(exp(fe + .ud_newdata_offset(fit, newdata, offset) + X %*% fit$coefficients)); a <- fit$alpha
  } else if (inherits(fit, "hurdle_cpb")) {
    Xb <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$intensity$levels, fit$intensity$contrasts,
                             names(fit$int_beta)) %*% fit$int_beta + .ud_newdata_offset(fit$intensity, newdata, offset)
    lam <- if (is.null(fit$fe)) as.numeric(exp(Xb))
           else { fe <- fit$intensity$fe[as.character(newdata[[fit$fe]])]; fe[is.na(fe)] <- .fe_ref(fit$intensity)
                  as.numeric(exp(fe + Xb)) }
    a <- fit$intensity$alpha
  } else {
    Xb <- .ud_newdata_matrix(fit$int_terms, newdata, fit$int_xlev, fit$int_contrasts, names(fit$coefficients)) %*%
      fit$coefficients
    lam <- if (is.null(fit$fe)) as.numeric(exp(.ud_newdata_offset(fit, newdata, offset) + Xb))
           else { fe <- fit$fe_hat[as.character(newdata[[fit$fe]])]; fe[is.na(fe)] <- fit$int_fe_ref
                  as.numeric(exp(fe + Xb)) }
    a <- fit$alpha
  }
  top <- .cpb_top(lam / (1 - a), .cpb_support(fit))
  if (length(top) != length(y)) stop("'newdata' must give the response and the covariates on the same rows.")
  w <- .ud_weights(weights, newdata, seq_len(nrow(newdata))); if (is.null(w)) w <- rep(1, length(y))
  ok <- is.finite(y) & is.finite(top)
  ex <- y[ok] > top[ok]
  data.frame(n = sum(w[ok]), exceed = sum(w[ok][ex]), share = sum(w[ok][ex]) / sum(w[ok]), row.names = NULL)
}

#' Profile-likelihood interval for the dispersion parameter alpha
#'
#' The interval inverts the likelihood-ratio test of `alpha`, re-maximizing the
#' coefficients at each value. By default the cut is the chi-square one (a
#' first-order interval). On King's support it covers between 0.940 and 0.956 in
#' simulations from the CPB. On the floor support it covers about 0.90 to 0.94,
#' with most misses on the upper side, because the estimate of `alpha` is biased
#' toward zero there ([calibrate_alpha()] explains the mechanism).
#' A fit that went through [calibrate_alpha()] gets the interval calibrated by
#' parametric bootstrap instead. Either interval is model-based: it assumes the
#' CPB and independent observations.
#'
#' @param object A `"cpb"` object.
#' @param level Confidence level (default 0.95).
#' @return A length-2 numeric vector (`lower`, `upper`) with attributes `alpha` (the
#'   point estimate), `method` (first-order or calibrated) and `boundary` (`TRUE`
#'   when the profile does not fall to the cut on the lower side: on King's
#'   support, whose trace runs down to `alpha = 1e-6`, the lower limit is then
#'   0; on the floor support, whose trace stops at `alpha = 0.005`, it is that
#'   bound).
#' @examples
#' set.seed(7); x <- rnorm(200)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
#' alpha_confint(fit)
#' @export
alpha_confint <- function(object, level = 0.95) {
  ## the profile re-maximizes the coefficients of a pooled fit; a fixed-effects fit would also need its unit
  ## effects re-concentrated at each alpha, which this profile does not do
  if (!inherits(object, "cpb"))
    stop("alpha_confint() takes a \"cpb\" fit. For a cpb_fe fit, the same likelihood is a cpb() fit with the ",
         "unit factor as dummies (cpb(y ~ x + factor(unit), truncated = ...)), whose interval this function gives.",
         call. = FALSE)
  aci <- .cpb_alpha_ci(object, level = level)
  structure(aci[c("lower", "upper")], alpha = object$alpha, method = attr(aci, "method"),
            boundary = aci[["boundary"]] == 1)
}

#' Incidence rate ratios for a CPB fit
#'
#' The ratios are \eqn{\exp(\beta)}, ratios of the CPB's rate parameter
#' \eqn{\lambda}. The CPB's mean equals \eqn{\lambda} only where
#' \eqn{\lambda/(1-\alpha)} is an integer, and a zero-truncated fit (the `cpb()`
#' default) models \eqn{E(Y \mid Y \ge 1)}, so a ratio of expected counts can
#' differ from them; [first_difference()] reports the exact mean contrast.
#'
#' @param object A `"cpb"` object. With `se = "bootstrap"` at fit time the
#'   ratios carry bootstrap percentile intervals; otherwise point estimates are
#'   returned with `method = "none"`, as for every other `irr()` method.
#' @param level Confidence level (default 0.95).
#' @return A `"ud_irr"` data frame -- the package-wide rate-ratio contract
#'   (columns `term`, `equation`, `ratio`, `estimate`, `lower`, `upper`,
#'   `method`) shared by every `irr()` method; here with bootstrap percentile
#'   intervals from the stored draws.
#' @param ... Further arguments passed to methods.
#' @examples
#' \donttest{
#' set.seed(8); x <- rnorm(200)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
#' irr(fit)
#' }
#' @export
irr <- function(object, ...) UseMethod("irr")

#' @rdname irr
#' @method irr cpb
#' @export
irr.cpb <- function(object, level = 0.95, ...) .ud_irr_draws(object, object$coefficients, level)

#' First difference for a CPB fit
#'
#' The effect on a quantity of interest of moving one covariate `from` one value `to`
#' another, holding the other covariates at their means (or at the values `at`
#' gives them).
#'
#' Every `first_difference()` method in the package takes its interval by one
#' rule. When the fit carries bootstrap draws (`se = "bootstrap"` at fit
#' time), the quantity is computed at every draw, jointly over the
#' coefficients and the dispersion parameter (and, for a fixed-effects fit, the
#' reference unit's effect), and the percentile interval is taken; the draws of
#' a fixed-effects fit are first recentered on the reported estimate, because
#' its replicates share the incidental-parameters bias. Otherwise the interval
#' is the delta method's over the fit's full covariance, the dispersion
#' parameter included (the matched count families' analytic, robust or
#' clustered covariance; a hurdle's two margins, which are independent). A fit
#' with neither gets the point estimate with `method = "none"`.
#'
#' @param object A `"cpb"` object. With `se = "bootstrap"` at fit time the
#'   difference carries a bootstrap percentile interval; otherwise the point
#'   estimate is returned with `method = "none"`.
#' @param variable Name of a model-matrix column to vary.
#' @param from,to The two values of `variable` to contrast.
#' @param quantity `"mean"` (E(Y); for a zero-truncated fit the conditional
#'   mean E(Y | Y >= 1)), `"ceiling"` (lambda/(1-alpha)), or `"prob"`
#'   (P(Y = `y`)).
#' @param y The count value for `quantity = "prob"`.
#' @param level Confidence level (default 0.95).
#' @param at Optional named numeric vector fixing other model-matrix columns at
#'   given values instead of their means, e.g. `at = c(z = 1, regionB = 0)` (a
#'   factor's dummies are columns, so a reference level sets them to 0).
#' @return A `"ud_fd"` data frame -- the package-wide first-difference contract
#'   (columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`) --
#'   with one row for the requested quantity and a bootstrap percentile
#'   interval on the difference (`method = "bootstrap (stored)"`). Every
#'   `first_difference()` method in the package returns this same shape.
#' @param ... Further arguments passed to methods; unknown arguments error.
#' @examples
#' \donttest{
#' set.seed(1); x <- rnorm(250)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(250, N, 0.5)
#' fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
#' first_difference(fit, "x", from = -1, to = 1, quantity = "mean")
#' }
#' @export
first_difference <- function(object, ...) UseMethod("first_difference")

#' @rdname first_difference
#' @method first_difference cpb
#' @export
first_difference.cpb <- function(object, variable, from, to,
                             quantity = c("mean", "ceiling", "prob"), y = NULL,
                             level = 0.95, at = NULL, ...) {
  .fd_dots(...); quantity <- match.arg(quantity)
  .fd_build(.fd_spec(object, variable, from, to, at = at, quantity = quantity, y = y), level)
}
