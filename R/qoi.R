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
#' @param type One of `"response"` (the mean lambda), `"link"` (the linear
#'   predictor), `"ceiling"` (the implied ceiling lambda/(1-alpha)), or `"prob"`
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
                        type = c("response", "link", "ceiling", "prob"), at = NULL,
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
    mf <- model.frame(Terms, newdata, xlev = object$levels)
    X  <- model.matrix(Terms, mf, contrasts.arg = object$contrasts)
    off <- if (is.null(offset)) rep_len(0, nrow(X))
           else as.numeric(if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset)
  }
  eta <- as.numeric(off + X %*% object$coefficients); lam <- exp(eta)
  switch(type,
    link     = eta,
    response = lam,
    ceiling  = lam / (1 - object$alpha),
    prob = {
      if (is.null(at)) stop("For type = \"prob\", supply 'at' (the count value y).")
      yv <- if (length(at) == 1) rep(at, length(lam)) else at
      if (length(yv) != length(lam)) stop("'at' must be length 1 or nrow(newdata).")
      .cpb_prob_at(lam, object$alpha, yv, object$truncated, object$max.support)
    })
}

#' Implied ceiling with a profile-likelihood interval
#'
#' Returns the observation- (or profile-) specific ceiling lambda/(1-alpha), with an
#' interval propagating the profile-likelihood uncertainty in `alpha` at the fitted
#' mean. (Coefficient uncertainty in lambda is not propagated here; use
#' [first_difference()] with `quantity = "ceiling"` for a fully bootstrapped contrast.)
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
  lam <- predict(object, newdata = newdata, type = "response")
  aci <- .cpb_alpha_profile_ci(object, level = level)
  data.frame(lambda  = lam,
             ceiling = lam / (1 - object$alpha),
             lower   = lam / (1 - aci["lower"]),
             upper   = lam / (1 - aci["upper"]),
             row.names = NULL)
}

#' Profile-likelihood interval for the dispersion parameter alpha
#'
#' @param object A `"cpb"` object.
#' @param level Confidence level (default 0.95).
#' @return A length-2 numeric vector (`lower`, `upper`) with attributes `alpha` (the
#'   point estimate) and `boundary` (`TRUE` if the lower bound is at the feasibility
#'   boundary, i.e. strong underdispersion, where the interval is one-sided).
#' @examples
#' set.seed(7); x <- rnorm(300)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "none")
#' alpha_confint(fit)
#' @export
alpha_confint <- function(object, level = 0.95) {
  aci <- .cpb_alpha_profile_ci(object, level = level)
  structure(aci[c("lower", "upper")], alpha = object$alpha,
            boundary = aci["boundary"] == 1)
}

#' Incidence rate ratios for a CPB fit
#'
#' @param object A `"cpb"` object fit with `se = "bootstrap"`.
#' @param level Confidence level (default 0.95).
#' @return A `"ud_irr"` data frame -- the package-wide rate-ratio contract
#'   (columns `term`, `equation`, `ratio`, `estimate`, `lower`, `upper`,
#'   `method`) shared by every `irr()` method; here with bootstrap percentile
#'   intervals from the stored draws.
#' @param ... Further arguments passed to methods.
#' @examples
#' \donttest{
#' set.seed(8); x <- rnorm(300)
#' N <- pmax(round(exp(1.5 + 0.4 * x) / 0.5), 1); y <- rbinom(300, N, 0.5)
#' fit <- cpb(y ~ x, data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 100)
#' irr(fit)
#' }
#' @export
irr <- function(object, ...) UseMethod("irr")

#' @rdname irr
#' @method irr cpb
#' @export
irr.cpb <- function(object, level = 0.95, ...) {
  if (is.null(object$boot)) stop("Bootstrap required; refit with se = \"bootstrap\".")
  ok <- object$boot[complete.cases(object$boot), , drop = FALSE]
  a  <- (1 - level) / 2
  ci <- t(apply(ok[, 1:object$p, drop = FALSE], 2, function(col) quantile(exp(col), c(a, 1 - a))))
  .ud_irr(names(object$coefficients), "count", "IRR", exp(unname(object$coefficients)),
          lower = unname(ci[, 1]), upper = unname(ci[, 2]),
          method = "bootstrap (stored)")
}

#' King-style first difference for a CPB fit
#'
#' The effect on a quantity of interest of moving one covariate `from` one value `to`
#' another, holding the other covariates at their means, with a bootstrap percentile
#' interval (Tomz, Wittenberg & King style, using the model's bootstrap draws).
#'
#' @param object A `"cpb"` object fit with `se = "bootstrap"`.
#' @param variable Name of a model-matrix column to vary.
#' @param from,to The two values of `variable` to contrast.
#' @param quantity `"mean"` (E(Y)), `"ceiling"` (lambda/(1-alpha)), or `"prob"`
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
#' set.seed(1); x <- rnorm(400)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(400, N, 0.5)
#' fit <- cpb(y ~ x, data = data.frame(y = y, x = x)[y > 0, ], se = "bootstrap", B = 200)
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
  if (is.null(object$boot)) stop("Bootstrap draws required; refit with se = \"bootstrap\".")
  if (!variable %in% colnames(object$X))
    stop("'variable' must name a model-matrix column: ",
         paste(colnames(object$X), collapse = ", "))
  x0 <- colMeans(object$X); xf <- x0; xt <- x0
  xf[variable] <- from; xt[variable] <- to
  qfun <- function(beta, alpha) {
    lf <- exp(sum(xf * beta)); lt <- exp(sum(xt * beta))
    if (quantity == "mean")         c(lf, lt)
    else if (quantity == "ceiling") c(lf / (1 - alpha), lt / (1 - alpha))
    else {
      if (is.null(y)) stop("quantity = \"prob\" requires 'y'.")
      c(.cpb_prob_at(lf, alpha, y, object$truncated, object$max.support),
        .cpb_prob_at(lt, alpha, y, object$truncated, object$max.support))
    }
  }
  pt  <- qfun(object$coefficients, object$alpha)
  ok  <- object$boot[complete.cases(object$boot), , drop = FALSE]
  fdb <- apply(ok, 1, function(r) { q <- qfun(r[1:object$p], r[object$p + 1]); q[2] - q[1] })
  a   <- (1 - level) / 2
  .ud_fd(component = quantity, from = pt[1], to = pt[2],
         lower = unname(quantile(fdb, a)), upper = unname(quantile(fdb, 1 - a)),
         method = "bootstrap (stored)")
}
