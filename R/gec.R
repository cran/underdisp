## ---------------------------------------------------------------------------
## King's generalized event count (Katz-family) regression: one dispersion
## parameter delta (Var/Mean) spanning under-, equi-, and over-dispersion.
## ---------------------------------------------------------------------------

#' Generalized event count (Katz-family) regression
#'
#' Fits King's generalized event count model, a Katz-family count regression
#' whose single dispersion parameter `delta` (the variance-to-mean ratio) is
#' estimated freely and spans underdispersion (`delta < 1`, a finite-support
#' member; the continuous parameter binomial is this cell), equidispersion
#' (`delta = 1`, Poisson), and overdispersion (`delta > 1`, the negative
#' binomial). Unlike [cpb()], which fixes the direction of dispersion to
#' under, `gec()` lets the data choose. The Katz recursion delivers exact first
#' and second moments---the Winkelmann--Signorino--King correction realized
#' directly---and the likelihood is evaluated in C++.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param truncated Logical; if `TRUE`, fit the zero-truncated GEC (all `Y >= 1`),
#'   the intensity model of [hurdle_gec()].
#' @param se Coefficient inference: `"none"` (default; fast) or `"bootstrap"`.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector) for a cluster/block bootstrap; see [cpb()].
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param max.support Guard on the maximum evaluated support.
#' @param maxit,reltol Optimizer controls.
#' @return An object of class `"gec"` with `coefficients`, `delta` (the estimated
#'   dispersion), `loglik`, bootstrap standard errors, and bookkeeping.
#' @examples
#' set.seed(1); x <- rnorm(400)
#' y <- rpois(400, exp(1 + 0.5 * x))
#' gec(y ~ x, data = data.frame(y = y, x = x), se = "none")
#' @seealso [cpb()], [compare_dispersion()]
#' @export
gec <- function(formula, data, truncated = FALSE, se = c("none", "bootstrap"), B = 500, cluster = NULL,
                offset = NULL, max.support = 500, maxit = 20000, reltol = 1e-8) {
  se <- match.arg(se); cl <- match.call()
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y  <- as.integer(stats::model.response(mf)); X <- stats::model.matrix(formula, mf)
  n  <- length(Y); p <- ncol(X)
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")
  off <- rep_len(0, n)                                          # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    off <- as.numeric(ov[match(rownames(mf), rownames(data))])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  clid <- NULL
  if (!is.null(cluster)) {
    cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    if (is.null(cv)) stop("'cluster' must be a column name in 'data' or a vector.")
    clid <- cv[match(rownames(mf), rownames(data))]
    if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
  }
  ms <- as.integer(max.support)
  objf <- if (!truncated) function(par, Xm, Ym, om) gec_nll_cpp(par, Xm, Ym, om, ms)
          else function(par, Xm, Ym, om) {           # zero-truncated Katz likelihood (P(Y>=1))
            m <- gec_lp0_cpp(par, Xm, Ym, om, ms)
            if (any(m[, 2] >= 1 - 1e-12)) return(1e10)
            v <- -sum(pmax(m[, 1], -700) - log1p(-m[, 2])); if (is.finite(v)) v else 1e10
          }
  fitone <- function(Xm, Ym, om) {
    b0 <- tryCatch({ v <- stats::glm.fit(Xm, Ym, offset = om, family = stats::poisson())$coefficients
                     v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
    best <- NULL
    for (ld in c(-0.5, 0, 0.5)) {                     # under / Poisson / over starts
      op <- tryCatch(stats::optim(c(b0, ld), function(par) objf(par, Xm, Ym, om),
                                  method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)),
                     error = function(e) NULL)
      if (!is.null(op) && op$value < 1e9 && (is.null(best) || op$value < best$value)) best <- op
    }
    best
  }
  fit <- fitone(X, Y, off)
  if (is.null(fit)) stop("All optimization starts failed.")
  beta   <- setNames(fit$par[1:p], colnames(X)); delta <- exp(fit$par[p + 1]); loglik <- -fit$value
  mu     <- as.numeric(exp(off + X %*% beta))
  fitted <- if (!truncated) gec_mean_cpp(mu, delta, ms)
            else { p0 <- gec_lp0_cpp(c(beta, log(delta)), X, Y, off, ms)[, 2]; mu / pmax(1 - p0, 1e-8) }

  se.beta <- setNames(rep(NA_real_, p), colnames(X)); se.delta <- NA_real_; boot <- NULL; ci.beta <- NULL
  if (se == "bootstrap") {
    grp <- if (!is.null(clid)) split(seq_len(n), clid) else NULL
    boot <- matrix(NA_real_, B, p + 1)
    for (b in seq_len(B)) {
      idx <- if (is.null(grp)) sample.int(n, n, replace = TRUE)
             else unlist(grp[sample.int(length(grp), length(grp), replace = TRUE)], use.names = FALSE)
      fb <- fitone(X[idx, , drop = FALSE], Y[idx], off[idx])
      if (!is.null(fb)) boot[b, ] <- c(fb$par[1:p], exp(fb$par[p + 1]))
    }
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta  <- setNames(apply(ok[, 1:p, drop = FALSE], 2, stats::sd), colnames(X))
      se.delta <- stats::sd(ok[, p + 1])
      ci.beta  <- t(apply(ok[, 1:p, drop = FALSE], 2, stats::quantile, c(.025, .975)))
      rownames(ci.beta) <- colnames(X)
    } else warning("Too few bootstrap resamples converged for stable inference.")
  }

  structure(list(coefficients = beta, delta = delta, dispersion = delta, loglik = loglik,
                 se.beta = se.beta, se.delta = se.delta, boot = boot, ci.beta = ci.beta,
                 fitted.values = fitted, linear.predictors = as.numeric(off + X %*% beta), offset = off,
                 n = n, df = p + 1, p = p, se.type = se, clustered = !is.null(clid),
                 truncated = truncated, max.support = ms, formula = formula, terms = attr(mf, "terms"),
                 X = X, Y = Y, call = cl), class = "gec")
}

#' @method print gec
#' @export
print.gec <- function(x, ...) {
  disp <- if (x$delta < 0.97) "underdispersed" else if (x$delta > 1.03) "overdispersed" else "~ equidispersed"
  cat("Generalized event count (Katz family) regression\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta)))
    print(round(cbind(Estimate = x$coefficients, `Std. Error` = x$se.beta), 4))
  else print(round(x$coefficients, 4))
  cat(sprintf("\ndispersion delta (Var/Mean) = %.3f  [%s]", x$delta, disp))
  if (is.finite(x$se.delta)) cat(sprintf("  (SE %.3f)", x$se.delta))
  cat(sprintf("\nlogLik = %.2f,  n = %d\n", x$loglik, x$n))
  invisible(x)
}

#' @method coef gec
#' @export
coef.gec <- function(object, ...) object$coefficients

#' @method logLik gec
#' @export
logLik.gec <- function(object, ...) {
  val <- object$loglik; attr(val, "df") <- object$df; attr(val, "nobs") <- object$n
  class(val) <- "logLik"; val
}

#' @method nobs gec
#' @export
nobs.gec <- function(object, ...) object$n

#' @method fitted gec
#' @export
fitted.gec <- function(object, ...) object$fitted.values

#' @method vcov gec
#' @export
vcov.gec <- function(object, ...) {
  if (is.null(object$boot)) return(NULL)                        # no bootstrap requested
  p  <- length(object$coefficients)                            # covariate block (works for gec AND gec_fe)
  ok <- object$boot[stats::complete.cases(object$boot), 1:p, drop = FALSE]
  if (nrow(ok) < 2L) return(NULL)
  v <- stats::cov(ok); dimnames(v) <- list(names(object$coefficients), names(object$coefficients)); v
}

#' Predictions from a generalized event count fit
#'
#' @param object A `"gec"` object.
#' @param newdata Optional covariate profiles.
#' @param type `"response"` (the mean), `"link"` (the log mean), or `"prob"`
#'   (the probability of the count `at`).
#' @param at Count value(s) for `type = "prob"`.
#' @param offset Optional offset (log scale) for the count component when
#'   predicting on `newdata`; a numeric vector or a column name in `newdata`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict gec
#' @export
predict.gec <- function(object, newdata = NULL, type = c("response", "link", "prob"), at = NULL,
                        offset = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    X <- object$X; off <- if (!is.null(object$offset)) object$offset else rep_len(0, nrow(X))
  } else {
    Terms <- stats::delete.response(object$terms)
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X <- stats::model.matrix(Terms, newdata)
    off <- if (is.null(offset)) rep_len(0, nrow(X))
           else as.numeric(if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset)
  }
  mu <- as.numeric(exp(off + X %*% object$coefficients))
  switch(type,
    link = log(mu),
    response = gec_mean_cpp(mu, object$delta, object$max.support),
    prob = {
      if (is.null(at)) stop("For type = \"prob\", supply 'at' (the count value).")
      yv <- if (length(at) == 1) rep(as.integer(at), length(mu)) else as.integer(at)
      if (length(yv) != length(mu)) stop("'at' must be length 1 or nrow(newdata).")
      m <- gec_lp0_cpp(c(object$coefficients, log(object$delta)), X, yv, off, object$max.support)
      pk <- exp(pmax(m[, 1], -700))
      if (isTRUE(object$truncated)) pk <- ifelse(yv == 0L, 0, pk / (1 - m[, 2]))  # condition on Y > 0
      pk
    })
}

## ---------------------------------------------------------------------------
## GEC with high-dimensional unit fixed effects (concentrated likelihood)
## ---------------------------------------------------------------------------

#' GEC (Katz-family) regression with high-dimensional unit fixed effects
#'
#' Fits [gec()] with a full set of unit fixed effects, concentrating (profiling)
#' out the unit intercepts by a one-dimensional inner maximization per unit, so
#' the outer optimizer handles only the covariate coefficients and the free
#' dispersion parameter `delta`. This is the [cpb_fe()] concentration generalized
#' to the whole Katz family: the panel can be under-, equi-, or overdispersed and
#' the direction is estimated, not presumed.
#'
#' Limitation: unlike [cpb_fe()], `gec_fe()` has no `truncated` argument -- a
#' concentrated zero-truncated GEC is not currently implemented. For a
#' zero-truncated GEC with fixed effects, enter the unit factor as dummies in
#' [gec()]'s formula (feasible for moderate unit counts).
#'
#' Short panels: `delta` and the fixed effects carry the incidental-parameters
#' bias of nonlinear fixed-effects estimation, of order 1/T; treat `delta`
#' cautiously below about T = 30. The covariate coefficients are not materially
#' affected. `bias_correct = "jackknife"` removes the leading 1/T term by the
#' split-panel jackknife exactly as in [cpb_fe()]: refit on each unit's temporal
#' halves and report `2 * full - mean(halves)`, with the unit effects and fitted
#' values re-concentrated at the corrected parameters and `logLik`/`AIC` kept at
#' the maximum-likelihood fit (uncorrected estimates in `$uncorrected`). The
#' same time-homogeneity validity gate applies: on panels where the two halves
#' do not estimate a common parameter the correction is REFUSED with a warning
#' and the maximum-likelihood fit is returned (see [cpb_fe()], Details).
#'
#' @param formula A model formula for the covariates only (no unit factor, no
#'   intercept; the fixed effects absorb it).
#' @param data A data frame.
#' @param fe Name of the column holding the unit identifier.
#' @param se Inference for the covariate coefficients: `"none"` (default) or
#'   `"bootstrap"`, a pairs/cluster bootstrap over units.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Cluster for the bootstrap; `NULL` (default) resamples the
#'   fixed-effects units. See [cpb_fe()].
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param max.support Guard on the maximum evaluated support.
#' @param inner_it Golden-section iterations for each unit's inner maximization.
#' @param maxit,reltol Outer optimizer controls.
#' @param bias_correct `"none"` (default) or `"jackknife"`, the split-panel
#'   jackknife correction for the 1/T incidental-parameters bias (see Details).
#' @return An object of class `c("gec_fe", "gec")`.
#' @examples
#' \donttest{
#' set.seed(5)
#' d <- do.call(rbind, lapply(1:30, function(i) {
#'   x <- rnorm(10); N <- pmax(round(exp(1 + rnorm(1, 0, 0.4) + 0.3 * x) / 0.5), 1)
#'   data.frame(unit = i, x = x, y = rbinom(10, N, 0.5))
#' }))
#' gec_fe(y ~ x, data = d, fe = "unit")   # delta ~ 0.5: underdispersed
#' }
#' @seealso [gec()], [cpb_fe()]
#' @export
gec_fe <- function(formula, data, fe, se = c("none", "bootstrap"), B = 500, cluster = NULL,
                   offset = NULL, max.support = 500L, inner_it = 30L, maxit = 3000L, reltol = 1e-7,
                   bias_correct = c("none", "jackknife")) {
  se <- match.arg(se); bias_correct <- match.arg(bias_correct)
  if (!fe %in% names(data)) stop("'fe' must name a column of 'data'.")
  if (!is.null(offset) && !(is.character(offset) && length(offset) == 1L)) {
    data[[".gec_off"]] <- as.numeric(offset); offset <- ".gec_off"   # vector -> column, sorts with data
  }
  data <- data[order(data[[fe]]), , drop = FALSE]
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y  <- as.integer(stats::model.response(mf))
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  off <- if (is.null(offset)) rep_len(0, length(Y))
         else as.numeric(data[[offset]][match(rownames(mf), rownames(data))])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  Xf <- stats::model.matrix(formula, mf)
  keep <- setdiff(colnames(Xf), "(Intercept)")
  if (!length(keep))
    stop("Provide at least one covariate in 'formula'; the fixed effects are the intercepts.")
  X  <- Xf[, keep, drop = FALSE]
  uf <- factor(data[[fe]]); nu <- nlevels(uf)
  ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
  p  <- ncol(X); ms <- as.integer(max.support); it <- as.integer(inner_it)
  bs <- tryCatch({ v <- stats::glm.fit(cbind(1, X), Y, offset = off, family = stats::poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  nll <- function(par) gec_fe_nll_cpp(par, X, Y, off, ustart, nu, ms, it)
  cand <- lapply(c(-0.5, 0, 0.5), function(ld)
    tryCatch(stats::optim(c(bs, ld), nll, method = "Nelder-Mead",
                          control = list(maxit = maxit, reltol = reltol)), error = function(e) NULL))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) stop("All optimization starts failed.")
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  beta <- setNames(best$par[1:p], keep); delta <- exp(best$par[p + 1])
  fe_hat <- gec_fe_intercepts_cpp(best$par, X, Y, off, ustart, nu, ms, it)
  rate <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
  fitted <- gec_mean_cpp(rate, delta, ms)

  se.beta <- setNames(rep(NA_real_, p), keep); ci.beta <- NULL; boot <- NULL; clab <- NULL
  if (se == "bootstrap") {
    fev  <- as.character(data[[fe]])
    cval <- if (is.null(cluster)) fev
            else if (is.character(cluster) && length(cluster) == 1L) as.character(data[[cluster]])
            else as.character(cluster)
    clab <- if (is.null(cluster)) fe else if (is.character(cluster) && length(cluster) == 1L) cluster else "custom"
    grp  <- split(seq_along(Y), cval); boot <- matrix(NA_real_, B, p)
    for (b in seq_len(B)) {
      gs   <- sample.int(length(grp), length(grp), replace = TRUE)
      rows <- unlist(grp[gs], use.names = FALSE)
      bunit <- unlist(lapply(seq_along(gs), function(k) paste0(k, "_", fev[grp[[gs[k]]]])), use.names = FALSE)
      bd <- data[rows, , drop = FALSE]; bd[[".bootunit"]] <- bunit
      fb <- tryCatch(gec_fe(formula, data = bd, fe = ".bootunit", se = "none", offset = offset,
                            max.support = max.support, inner_it = inner_it, maxit = maxit, reltol = reltol),
                     error = function(e) NULL)
      if (!is.null(fb)) boot[b, ] <- fb$coefficients[keep]
    }
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta <- setNames(apply(ok, 2, stats::sd), keep)
      z <- stats::qnorm(0.975)
      ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta); rownames(ci.beta) <- keep
    } else warning("Too few bootstrap resamples converged for stable inference.")
    attr(boot, "nboot_ok") <- nrow(ok)
  }

  ## split-panel jackknife: refit on each unit's temporal halves and remove the
  ## leading 1/T incidental-parameters bias term (see .fe_jackknife)
  uncorrected <- NULL
  if (bias_correct == "jackknife") {
    dest <- data[match(rownames(mf), rownames(data)), , drop = FALSE]
    refit <- function(dd) tryCatch({
      f <- gec_fe(formula, data = dd, fe = fe, se = "none", offset = offset,
                  max.support = max.support, inner_it = inner_it, maxit = maxit, reltol = reltol)
      c(f$coefficients[keep], f$delta)
    }, error = function(e) NULL)
    jk <- .fe_jackknife(dest, uf, refit, c(beta, delta),
                        resid  = (Y - fitted) / sqrt(pmax(delta * fitted, 1e-12)),
                        torder = stats::ave(seq_along(Y), as.integer(uf), FUN = seq_along),
                        disp_lower = 0.01, disp_upper = Inf)
    if (!is.null(jk$refused)) {
      warning("bias_correct = \"jackknife\" REFUSED: ", jk$refused,
              ". Returning the uncorrected maximum-likelihood estimates.")
      bias_correct <- "none"
    } else {
      for (w in jk$warn) warning("bias_correct = \"jackknife\": ", w, ".")
      uncorrected <- list(coefficients = beta, delta = delta)
      beta  <- setNames(jk$par[seq_len(p)], keep)
      delta <- unname(jk$par[p + 1L])
      fe_hat <- gec_fe_intercepts_cpp(c(beta, log(delta)), X, Y, off, ustart, nu, ms, it)
      rate   <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
      fitted <- gec_mean_cpp(rate, delta, ms)
      if (!is.null(ci.beta)) {                        # recenter the normal-approx interval
        z <- stats::qnorm(0.975)
        ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
        rownames(ci.beta) <- keep
      }
    }
  }

  structure(list(coefficients = beta, delta = delta, dispersion = delta, loglik = -best$value,
                 se.beta = se.beta, ci.beta = ci.beta, boot = boot, se.type = se, cluster = clab,
                 fe = setNames(fe_hat, levels(uf)), fitted.values = fitted, rate = rate, offset = off,
                 xref = colMeans(X), fe_ref = mean(fe_hat),
                 linear.predictors = log(rate), n = length(Y), n_units = nu, Y = Y,
                 df = p + nu + 1L, max.support = ms, formula = formula, fe_var = fe,
                 bias_correct = bias_correct, uncorrected = uncorrected,
                 call = match.call()), class = c("gec_fe", "gec"))
}

#' @method print gec_fe
#' @export
print.gec_fe <- function(x, ...) {
  disp <- if (x$delta < 0.97) "underdispersed" else if (x$delta > 1.03) "overdispersed" else "~ equidispersed"
  cat("GEC (Katz family) regression with", x$n_units, "unit fixed effects (concentrated likelihood)\n")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta))) {
    cat("Coefficients (", x$cluster, "-clustered bootstrap SEs):\n", sep = "")
    print(round(cbind(Estimate = x$coefficients, `Std. Error` = x$se.beta[names(x$coefficients)]), 4))
  } else { cat("Coefficients:\n"); print(round(x$coefficients, 4)) }
  cat(sprintf("\ndispersion delta (Var/Mean) = %.3f  [%s],  n = %d\n", x$delta, disp, x$n))
  if (identical(x$bias_correct, "jackknife"))
    cat("Estimates are split-panel jackknife bias-corrected; logLik/AIC refer to the\nuncorrected maximum-likelihood fit. See ?gec_fe.\n")
  else
    cat("Note: delta is subject to incidental-parameters bias for short panels; see ?gec_fe.\n")
  invisible(x)
}

#' @method predict gec_fe
#' @export
predict.gec_fe <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) { rate <- object$rate } else {
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X <- stats::model.matrix(Terms, newdata)[, names(object$coefficients), drop = FALSE]
    rate <- as.numeric(exp(mean(object$fe) + X %*% object$coefficients))
  }
  switch(type, link = log(rate), response = gec_mean_cpp(rate, object$delta, object$max.support))
}
