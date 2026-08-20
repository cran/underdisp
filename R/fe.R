## ---------------------------------------------------------------------------
## CPB regression with high-dimensional unit fixed effects (concentrated likelihood)
## ---------------------------------------------------------------------------

#' CPB regression with high-dimensional unit fixed effects
#'
#' Fits the CPB with a full set of unit fixed effects by concentrating (profiling)
#' out the unit intercepts. Each unit's intercept is solved by a one-dimensional
#' inner maximization, so the outer optimizer handles only the covariate coefficients
#' and the dispersion parameter. This avoids an explicit unit-dummy design matrix and
#' scales to thousands of units, where dummy-based fitting exhausts memory. The
#' concentrated log-likelihood equals the full-dummy log-likelihood exactly.
#'
#' Short panels: the dispersion parameter `alpha` and the fixed effects are subject to
#' the incidental-parameters bias of nonlinear fixed-effects estimation. The bias in
#' `alpha` is downward and of order 1/T, where T is the per-unit number of
#' observations. It is small for T at least about 30 and should be treated cautiously
#' for short panels. The covariate coefficients are not materially affected.
#' `bias_correct = "jackknife"` removes the leading 1/T term by the split-panel
#' (half-panel) jackknife: the model is refit on the first and second temporal
#' halves of every unit's series (rows are split in the order supplied, which
#' should be temporal order) and the corrected estimate is
#' `2 * full - mean(halves)`. In the package's fixed-effects bias Monte Carlo the
#' alpha bias at T = 6/10/20/40 falls from -0.147/-0.096/-0.056/-0.033 to
#' -0.038/-0.022/-0.016/-0.009. With correction on, `coefficients` and `alpha`
#' are the corrected estimates (the maximum-likelihood values are kept in
#' `$uncorrected`), the unit effects and fitted values are re-concentrated at the
#' corrected parameters, and `logLik`/`AIC` continue to refer to the
#' maximum-likelihood fit.
#'
#' Validity gate: the split-panel identity requires the two half-panels to
#' estimate the same pseudo-true parameter (time-homogeneity; Dhaene & Jochmans
#' 2015, Review of Economic Studies). On trending or time-heterogeneous panels
#' the correction is invalid, so the function REFUSES it -- returning the
#' uncorrected maximum-likelihood fit with a warning naming the failed check --
#' whenever the corrected dispersion leaves its feasible space or the two
#' half-panel dispersion estimates disagree beyond what the 1/T bias can
#' explain. A refusal is diagnostic information about the panel, not an error.
#' The gate is deliberately powered over sized: in the package's calibration
#' Monte Carlo it refuses about 9 percent of genuinely time-homogeneous panels
#' (a conservative nuisance; the returned fit is exactly the ordinary
#' maximum-likelihood estimate) while catching 98 percent of dispersion regime
#' changes and 100 percent of smooth unmodeled trends -- the cases where an
#' uncaught correction would be silently wrong.
#'
#' One set of fixed effects is concentrated out. For two-way (e.g. unit and time)
#' fixed effects, add `factor(time)` to `formula` (entered as dummies), or use
#' [count_reg()], which supports two-way fixed effects natively.
#'
#' @param formula A model formula for the covariates only. Do not include the unit
#'   factor; supply it via `fe`. The intercept is absorbed by the fixed effects.
#' @param data A data frame.
#' @param fe Name of the column holding the unit identifier.
#' @param truncated Logical; fit the zero-truncated CPB if `TRUE`. **Note the
#'   sibling default differs:** `cpb_fe()` defaults to `FALSE` (whole-panel
#'   fits including zeros), while [cpb()] defaults to `TRUE` (positives-only
#'   fits). State `truncated` explicitly when moving a specification between
#'   the two.
#' @param se Inference for the covariate coefficients: `"none"` (default) or
#'   `"bootstrap"`, a pairs/cluster bootstrap.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Cluster for the bootstrap. `NULL` (default) resamples the
#'   fixed-effects units themselves (unit-clustered inference, the natural panel
#'   default); a column name or vector resamples those clusters while preserving
#'   the unit fixed effects. Duplicated blocks are relabeled so the concentrated
#'   likelihood treats them as distinct units.
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param max.support Guard on the maximum evaluated support.
#' @param inner_it Golden-section iterations for each unit's inner maximization.
#' @param maxit,reltol Outer optimizer controls.
#' @param bias_correct `"none"` (default) or `"jackknife"`, the split-panel
#'   jackknife correction for the 1/T incidental-parameters bias (see Details).
#' @return An object of class `"cpb_fe"` with `coefficients`, `alpha`, `loglik`, `fe`
#'   (the estimated unit intercepts), `fitted.values`, and the implied per-observation
#'   `ceiling`.
#' @seealso [cpb()], [ud_screen()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- do.call(rbind, lapply(1:60, function(i) {
#'   x <- rnorm(20); lam <- exp(rnorm(1, 0, 0.5) + 0.5 * x)
#'   N <- pmax(round(lam / 0.5), 1)
#'   data.frame(unit = i, x = x, y = rbinom(20, N, 0.5))
#' }))
#' fit <- cpb_fe(y ~ x, data = d, fe = "unit")
#' fit
#' }
#' @export
cpb_fe <- function(formula, data, fe, truncated = FALSE, se = c("none", "bootstrap"),
                   B = 500, cluster = NULL, offset = NULL, max.support = 500L,
                   inner_it = 30L, maxit = 3000L, reltol = 1e-7,
                   bias_correct = c("none", "jackknife")) {
  se <- match.arg(se); bias_correct <- match.arg(bias_correct)
  if (!fe %in% names(data)) stop("'fe' must name a column of 'data'.")
  if (!is.null(offset) && !(is.character(offset) && length(offset) == 1L)) {
    data[[".cpb_off"]] <- as.numeric(offset); offset <- ".cpb_off"   # vector -> column, so it sorts with data
  }
  data <- data[order(data[[fe]]), , drop = FALSE]
  mf <- model.frame(formula, data, na.action = na.omit)
  Y  <- as.integer(model.response(mf))
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")
  off <- if (is.null(offset)) rep_len(0, length(Y))
         else as.numeric(data[[offset]][match(rownames(mf), rownames(data))])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  Xf <- model.matrix(formula, mf)
  keep <- setdiff(colnames(Xf), "(Intercept)")
  if (!length(keep))
    stop("Provide at least one covariate in 'formula'; the fixed effects are the intercepts.")
  X  <- Xf[, keep, drop = FALSE]
  uf <- factor(data[[fe]]); nu <- nlevels(uf)
  ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
  p  <- ncol(X)
  bs <- tryCatch({ v <- glm.fit(cbind(1, X), Y, offset = off, family = poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  nll <- function(par) cpb_fe_nll_cpp(par, X, Y, off, ustart, nu, as.integer(max.support),
                                      truncated, as.integer(inner_it))
  cand <- lapply(c(0.4, 0.6), function(a0)
    tryCatch(optim(c(bs, qlogis(a0)), nll, method = "Nelder-Mead",
                   control = list(maxit = maxit, reltol = reltol)), error = function(e) NULL))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) stop("All optimization starts failed.")
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  beta  <- setNames(best$par[1:p], keep); alpha <- plogis(best$par[p + 1])
  fe_hat <- cpb_fe_intercepts_cpp(best$par, X, Y, off, ustart, nu, as.integer(max.support),
                                  truncated, as.integer(inner_it))
  lam <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))

  ## pairs/cluster bootstrap for the covariate coefficients: resample whole
  ## clusters (default = the fixed-effects unit) and relabel each sampled block
  ## so the concentration treats duplicated draws as distinct units.
  se.beta <- setNames(rep(NA_real_, p), keep); ci.beta <- NULL; boot <- NULL
  clab <- NULL
  if (se == "bootstrap") {
    fev  <- as.character(data[[fe]])
    cval <- if (is.null(cluster)) fev
            else if (is.character(cluster) && length(cluster) == 1L) as.character(data[[cluster]])
            else as.character(cluster)
    clab <- if (is.null(cluster)) fe else if (is.character(cluster) && length(cluster) == 1L) cluster else "custom"
    grp  <- split(seq_along(Y), cval)
    boot <- matrix(NA_real_, B, p)
    for (b in seq_len(B)) {
      gs   <- sample.int(length(grp), length(grp), replace = TRUE)
      rows <- unlist(grp[gs], use.names = FALSE)
      bunit <- unlist(lapply(seq_along(gs), function(k) paste0(k, "_", fev[grp[[gs[k]]]])),
                      use.names = FALSE)
      bd <- data[rows, , drop = FALSE]; bd[[".bootunit"]] <- bunit
      fb <- tryCatch(cpb_fe(formula, data = bd, fe = ".bootunit", truncated = truncated,
                            se = "none", offset = offset, max.support = max.support, inner_it = inner_it,
                            maxit = maxit, reltol = reltol), error = function(e) NULL)
      if (!is.null(fb)) boot[b, ] <- fb$coefficients[keep]
    }
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta <- setNames(apply(ok, 2, stats::sd), keep)
      ## normal-approximation interval centered on the point estimate: FE bootstrap
      ## replicates share the incidental-parameters bias, so a percentile interval
      ## is off-center; the bootstrap SD is the robust quantity to keep.
      z <- stats::qnorm(0.975)
      ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
      rownames(ci.beta) <- keep
    } else warning("Too few bootstrap resamples converged for stable inference.")
    attr(boot, "nboot_ok") <- nrow(ok)
  }

  ## split-panel jackknife: refit on each unit's temporal halves and remove the
  ## leading 1/T incidental-parameters bias term (see .fe_jackknife)
  uncorrected <- NULL
  if (bias_correct == "jackknife") {
    dest <- data[match(rownames(mf), rownames(data)), , drop = FALSE]
    refit <- function(dd) tryCatch({
      f <- cpb_fe(formula, data = dd, fe = fe, truncated = truncated, se = "none",
                  offset = offset, max.support = max.support, inner_it = inner_it,
                  maxit = maxit, reltol = reltol)
      c(f$coefficients[keep], f$alpha)
    }, error = function(e) NULL)
    jk <- .fe_jackknife(dest, uf, refit, c(beta, alpha),
                        resid  = (Y - lam) / sqrt(pmax(alpha * lam, 1e-12)),
                        torder = stats::ave(seq_along(Y), as.integer(uf), FUN = seq_along),
                        disp_lower = 1e-3, disp_upper = 1 - 1e-3)
    if (!is.null(jk$refused)) {
      warning("bias_correct = \"jackknife\" REFUSED: ", jk$refused,
              ". Returning the uncorrected maximum-likelihood estimates.")
      bias_correct <- "none"
    } else {
      for (w in jk$warn) warning("bias_correct = \"jackknife\": ", w, ".")
      uncorrected <- list(coefficients = beta, alpha = alpha)
      beta  <- setNames(jk$par[seq_len(p)], keep)
      alpha <- unname(jk$par[p + 1L])
      fe_hat <- cpb_fe_intercepts_cpp(c(beta, qlogis(alpha)), X, Y, off, ustart, nu,
                                      as.integer(max.support), truncated, as.integer(inner_it))
      lam <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
      if (!is.null(ci.beta)) {                        # recenter the normal-approx interval
        z <- stats::qnorm(0.975)
        ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
        rownames(ci.beta) <- keep
      }
    }
  }

  structure(list(coefficients = beta, alpha = alpha, loglik = -best$value,
                 se.beta = se.beta, ci.beta = ci.beta, boot = boot, se.type = se, cluster = clab,
                 fe = setNames(fe_hat, levels(uf)), fitted.values = lam, offset = off,
                 xref = colMeans(X), fe_ref = mean(fe_hat),
                 ceiling = lam / (1 - alpha), n = length(Y), n_units = nu, Y = Y,
                 df = p + nu + 1L, truncated = truncated, formula = formula,
                 fe_var = fe, bias_correct = bias_correct, uncorrected = uncorrected,
                 call = match.call()),
            class = "cpb_fe")
}

#' @export
print.cpb_fe <- function(x, ...) {
  cat("CPB regression with", x$n_units, "unit fixed effects (concentrated likelihood)\n")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta))) {
    cat("Coefficients (", x$cluster, "-clustered bootstrap SEs):\n", sep = "")
    print(round(cbind(Estimate = x$coefficients,
                      `Std. Error` = x$se.beta[names(x$coefficients)]), 4))
  } else { cat("Coefficients:\n"); print(round(x$coefficients, 4)) }
  cat("\nalpha (shape parameter):", round(x$alpha, 4),
      "  median implied bound:", round(stats::median(x$ceiling), 2), "\n")
  if (identical(x$bias_correct, "jackknife"))
    cat("Estimates are split-panel jackknife bias-corrected; logLik/AIC refer to the\nuncorrected maximum-likelihood fit. See ?cpb_fe.\n")
  else
    cat("Note: alpha is subject to incidental-parameters bias for short panels; see ?cpb_fe.\n")
  invisible(x)
}

#' Predictions from a fixed-effects CPB fit
#'
#' @param object A `"cpb_fe"` object.
#' @param newdata Optional covariate profiles. With `newdata`, predictions use the
#'   average unit (the mean fixed effect), since the panel's own unit effects do
#'   not apply to new rows.
#' @param type `"response"` (the mean), `"link"` (the log mean), or `"ceiling"`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict cpb_fe
#' @export
predict.cpb_fe <- function(object, newdata = NULL, type = c("response", "link", "ceiling"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    lam <- object$fitted.values
  } else {
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss))
      stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X  <- stats::model.matrix(Terms, newdata)[, names(object$coefficients), drop = FALSE]
    lam <- as.numeric(exp(mean(object$fe) + X %*% object$coefficients))
  }
  switch(type, link = log(lam), response = lam, ceiling = lam / (1 - object$alpha))
}
