## ---------------------------------------------------------------------------
## Quantities of interest for CPB fits
## ---------------------------------------------------------------------------

## CPB probability mass for one (lambda, alpha): P(Y = k), k = 0..K. If truncated,
## P(Y = 0) = 0 and the remaining mass is renormalized to sum to one.
## Support capped at max.support; numerics from the shared .cpb_pmf_core.
.cpb_pmf_row <- function(lambda, alpha, truncated, max.support = 500) {
  pr <- .cpb_pmf_core(lambda, alpha, kcap = as.integer(max.support))
  if (truncated) { p0 <- pr[1]; pr <- pr / (1 - p0); pr[1] <- 0 }
  pr
}

## P(Y = y) at given (lambda, alpha), vectorized over equal-length lambda and y.
.cpb_prob_at <- function(lambda, alpha, y, truncated, max.support = 500) {
  mapply(function(l, yy) {
    if (yy < 0 || yy != floor(yy)) return(0)
    pr <- .cpb_pmf_row(l, alpha, truncated, max.support)
    if (yy + 1 > length(pr)) 0 else pr[yy + 1]
  }, lambda, y)
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
#'   lambda), `"ceiling"` (the implied ceiling lambda/(1-alpha)), or `"prob"`
#'   (the probability that `Y` equals `at`).
#' @param at For `type = "prob"`, the count value(s) `y` whose probability is
#'   returned (length 1, or one per row of the prediction data).
#' @param offset Optional offset (log scale) for the `newdata` branch: a numeric
#'   vector or a column name in `newdata`. For `newdata = NULL` the fit's own
#'   offset is used.
#' @param ... Unused.
#' @return A numeric vector.
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
                        offset = NULL, ...) {
  type <- match.arg(type)
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
    off <- if (is.null(offset)) rep_len(0, nrow(X))
           else as.numeric(if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset)
  }
  bad <- if (is.null(newdata)) FALSE else .ud_na_rows(X, off)     # a missing covariate or offset predicts NA
  if (any(bad)) { X[bad, ] <- 0; off[bad] <- 0 }
  eta <- as.numeric(off + X %*% object$coefficients); lam <- exp(eta)
  out <- switch(type,
    link     = eta,
    rate     = lam,
    response = .cpb_mean(lam, object$alpha, isTRUE(object$truncated)),
    ceiling  = lam / (1 - object$alpha),
    prob = {
      if (is.null(at)) stop("For type = \"prob\", supply 'at' (the count value y).")
      yv <- if (length(at) == 1) rep(at, length(lam)) else at
      if (length(yv) != length(lam)) stop("'at' must be length 1 or nrow(newdata).")
      .cpb_prob_at(lam, object$alpha, yv, object$truncated, object$max.support)
    })
  .ud_mask(out, bad)
}

#' Implied ceiling with a profile-likelihood interval
#'
#' Returns the observation- (or profile-) specific ceiling lambda/(1-alpha), with
#' bounds that carry the limits of the interval for `alpha` (see [alpha_confint()];
#' calibrated when the fit went through [calibrate_alpha()]) to the ceiling at the
#' fitted rate. The rate is held at its estimate: coefficient uncertainty in lambda
#' is not propagated, so the bounds are not a confidence interval for the ceiling;
#' use [first_difference()] with `quantity = "ceiling"` for a fully bootstrapped
#' contrast.
#'
#' @param object A `"cpb"` object.
#' @param newdata Optional covariate profiles.
#' @param level Confidence level (default 0.95).
#' @return A data frame with `lambda`, `ceiling`, and `lower`/`upper` bounds.
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
  data.frame(lambda  = lam,
             ceiling = lam / (1 - object$alpha),
             lower   = lam / (1 - aci["lower"]),
             upper   = lam / (1 - aci["upper"]),
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
  data.frame(lambda = lam, ceiling = lam / (1 - object$alpha), row.names = NULL)
}

#' Profile-likelihood interval for the dispersion parameter alpha
#'
#' The interval inverts the likelihood-ratio test of `alpha`, re-maximizing the
#' coefficients at each value. By default the cut is the chi-square one (a
#' first-order interval); it covers about 0.90 to 0.94 in simulations from the
#' CPB, with nearly all misses on the upper side, because the estimate of
#' `alpha` is biased toward zero ([calibrate_alpha()] explains the mechanism).
#' A fit that went through [calibrate_alpha()] gets the interval calibrated by
#' parametric bootstrap instead. Either interval is model-based: it assumes the
#' CPB and independent observations.
#'
#' @param object A `"cpb"` object.
#' @param level Confidence level (default 0.95).
#' @return A length-2 numeric vector (`lower`, `upper`) with attributes `alpha` (the
#'   point estimate), `method` (first-order or calibrated) and `boundary` (`TRUE`
#'   when the profile has not fallen to the cut by `alpha = 0.005`, so that the
#'   lower limit is the parameter bound).
#' @examples
#' set.seed(7); x <- rnorm(200)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(200, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
#' alpha_confint(fit)
#' @export
alpha_confint <- function(object, level = 0.95) {
  aci <- .cpb_alpha_ci(object, level = level)
  structure(aci[c("lower", "upper")], alpha = object$alpha, method = attr(aci, "method"),
            boundary = aci[["boundary"]] == 1)
}

#' Incidence rate ratios for a CPB fit
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
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 40)
#' irr(fit)
#' }
#' @export
irr <- function(object, ...) UseMethod("irr")

#' @rdname irr
#' @method irr cpb
#' @export
irr.cpb <- function(object, level = 0.95, ...) {
  if (is.null(object$boot))
    return(.ud_irr_wald(object$coefficients, NULL, level, "count", "IRR"))
  ok <- object$boot[complete.cases(object$boot), , drop = FALSE]
  a  <- (1 - level) / 2
  ci <- t(apply(ok[, 1:object$p, drop = FALSE], 2, function(col) quantile(exp(col), c(a, 1 - a))))
  .ud_irr(names(object$coefficients), "count", "IRR", exp(unname(object$coefficients)),
          lower = unname(ci[, 1]), upper = unname(ci[, 2]),
          method = "bootstrap (stored)")
}

#' First difference for a CPB fit
#'
#' The effect on a quantity of interest of moving one covariate `from` one value `to`
#' another, holding the other covariates at their means, with a percentile
#' interval from the model's bootstrap draws when the fit carries them.
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
#' fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 60)
#' first_difference(fit, "x", from = -1, to = 1, quantity = "mean")
#' }
#' @export
first_difference <- function(object, ...) UseMethod("first_difference")

#' @rdname first_difference
#' @method first_difference cpb
#' @export
first_difference.cpb <- function(object, variable, from, to,
                             quantity = c("mean", "ceiling", "prob"), y = NULL,
                             level = 0.95, ...) {
  .fd_dots(...); quantity <- match.arg(quantity)
  if (!variable %in% colnames(object$X))
    stop("'variable' must name a model-matrix column: ",
         paste(colnames(object$X), collapse = ", "))
  .fd_check_terms(variable, colnames(object$X))
  off <- if (is.null(object$offset)) 0 else mean(object$offset)   # offset held at its mean
  x0 <- colMeans(object$X); xf <- x0; xt <- x0
  xf[variable] <- from; xt[variable] <- to
  qfun <- function(beta, alpha) {
    lf <- exp(off + sum(xf * beta)); lt <- exp(off + sum(xt * beta))
    if (quantity == "mean")         .cpb_mean(c(lf, lt), alpha, isTRUE(object$truncated))
    else if (quantity == "ceiling") c(lf / (1 - alpha), lt / (1 - alpha))
    else {
      if (is.null(y)) stop("quantity = \"prob\" requires 'y'.")
      c(.cpb_prob_at(lf, alpha, y, object$truncated, object$max.support),
        .cpb_prob_at(lt, alpha, y, object$truncated, object$max.support))
    }
  }
  pt  <- qfun(object$coefficients, object$alpha)
  if (is.null(object$boot))
    return(.ud_fd(component = quantity, from = pt[1], to = pt[2]))
  ok  <- object$boot[complete.cases(object$boot), , drop = FALSE]
  fdb <- apply(ok, 1, function(r) { q <- qfun(r[1:object$p], r[object$p + 1]); q[2] - q[1] })
  a   <- (1 - level) / 2
  .ud_fd(component = quantity, from = pt[1], to = pt[2],
         lower = unname(quantile(fdb, a)), upper = unname(quantile(fdb, 1 - a)),
         method = "bootstrap (stored)")
}
