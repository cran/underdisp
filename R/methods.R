## ---------------------------------------------------------------------------
## S3 methods for class "cpb"
## ---------------------------------------------------------------------------

#' @export
print.cpb <- function(x, ...) {
  cat("Continuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)" else "(untruncated)", "\n")
  cat("Call:  ", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n"); print(round(x$coefficients, 4))
  cat("\nalpha (Var/Mean at an integer ceiling):", round(x$alpha, 4),
      "  median largest feasible count:", stats::median(.cpb_top(x$ceiling, .cpb_support(x))), "\n")
  .cpb_lower_note(x$alpha)
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  if (isTRUE(x$support_binding)) cat("Note: the fitted ceiling reaches max.support; alpha is bounded by the guard, not the data.\n")
  invisible(x)
}

#' @method coef cpb
#' @export
coef.cpb <- function(object, ...) object$coefficients

#' @method vcov cpb
#' @export
vcov.cpb <- function(object, ...) {
  if (is.null(object$vcov)) return(NULL)                 # no bootstrap requested
  object$vcov
}

#' @method logLik cpb
#' @export
logLik.cpb <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")

#' @method nobs cpb
#' @export
nobs.cpb <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted

#' @method fitted cpb
#' @export
fitted.cpb <- function(object, ...) object$fitted.values

#' @method residuals cpb
#' @export
residuals.cpb <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$Y - object$fitted.values
  if (type == "pearson") {                          # exact variance of the fitted pmf
    v <- .cpb_moments(.cpb_rate(object), object$alpha, isTRUE(object$truncated), .cpb_support(object))$var
    r <- r / sqrt(pmax(v, 1e-12))
  }
  r
}

#' Summarize a CPB fit
#'
#' @param object A `"cpb"` object.
#' @param ... Unused.
#' @return An object of class `"summary.cpb"` with the coefficient table, the
#'   dispersion parameter and its profile-likelihood interval (first-order, or
#'   calibrated when the fit went through [calibrate_alpha()]), the implied ceiling,
#'   fit statistics, and the likelihood-ratio test against a (zero-truncated) Poisson.
#' @method summary cpb
#' @export
summary.cpb <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  aci <- .cpb_alpha_ci(object)
  ## the likelihood-ratio statistic against the (zero-truncated) Poisson limit and
  ## its asymptotic boundary-mixture p-value: a value the support guard pushed
  ## below zero is the boundary value 0, and the p-value is flagged as asymptotic
  ## because dispersion_test() calibrates boundary nulls by parametric bootstrap
  lr  <- .disp_boundary_lr(if (is.na(object$loglik.null)) NA_real_ else -2 * (object$loglik.null - object$loglik),
                           isTRUE(object$support_binding), .cpb_at_limit(object))
  if (is.finite(lr$p))
    lr$note <- c(lr$note, paste0("the p-value is asymptotic, which over-rejects in finite samples at this boundary; ",
                                 "dispersion_test() gives the parametric-bootstrap p-value"))
  out <- list(call = object$call, truncated = object$truncated, coefficients = ctab,
              alpha = object$alpha, alpha.ci = aci, ceiling = object$ceiling,
              loglik = object$loglik, aic = -2 * object$loglik + 2 * object$df,
              LR = lr$LR, LR.p = lr$p, LR.note = lr$note,
              support_binding = isTRUE(object$support_binding), converged = object$converged,
              n = object$n, n_dropped = object$n_dropped, nobs = nobs(object), se.type = object$se.type,
              support = .cpb_support(object),
              nboot_ok = if (!is.null(object$boot)) attr(object$boot, "nboot_ok") else NA)
  class(out) <- "summary.cpb"
  out
}

#' @export
print.summary.cpb <- function(x, ...) {
  cat("\nContinuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)\n" else "(untruncated)\n")
  cat("N =", x$n, if (!is.null(x$nobs) && x$nobs != x$n) paste0(" (weight total ", format(x$nobs), ")"),
      if (isTRUE(x$n_dropped > 0)) paste0(" (", x$n_dropped, " rows with a missing model variable dropped)"),
      "   inference:", x$se.type, "\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  calibrated <- isTRUE(grepl("calibrated", attr(x$alpha.ci, "method"), fixed = TRUE))
  cat("\nalpha =", round(x$alpha, 4),
      sprintf("  (%s profile %d%% CI: %.3g to %.3g%s)\n", if (calibrated) "calibrated" else "first-order", 95L,
              x$alpha.ci["lower"], x$alpha.ci["upper"],
              if (x$alpha.ci["boundary"] == 1)
                if (.cpb_fl(.cpb_support(x))) ", lower limit at the parameter bound"
                else ", the interval reaches alpha's limit 0 (the two-point distribution)" else ""))
  cat("Implied ceiling lambda/(1-alpha): median", round(stats::median(x$ceiling), 2),
      "  range", paste(round(range(x$ceiling), 2), collapse = " to "), "\n")
  cat("logLik =", round(x$loglik, 2), "   AIC =", round(x$aic, 2), "\n")
  if (!is.na(x$LR))
    cat(sprintf("LR vs %sPoisson (H0: alpha = 1): %.2f, p %s\n",
                if (x$truncated) "ZT-" else "", x$LR, format.pval(x$LR.p)))
  for (nt in x$LR.note) cat("Note: ", nt, ".\n", sep = "")
  .cpb_lower_note(x$alpha)
  if (isTRUE(x$support_binding))
    cat("Note: the fitted ceiling reaches max.support; alpha is bounded by the guard, not the data.\n")
  if (!is.na(x$nboot_ok))
    cat("(", x$nboot_ok, " bootstrap resamples converged)\n", sep = "")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' Confidence intervals for a CPB fit
#'
#' Coefficient intervals use the cold-multistart bootstrap percentile method. In
#' simulations from King's CPB (n = 300) the 95% interval for a slope covers 0.91
#' at alpha 0.3, 0.95 at 0.5 and 0.94 at 0.7, and on the floor support 0.95. The
#' interval for `alpha` is a
#' profile-likelihood interval. It is first-order by default: on King's support
#' it covers close to its level in simulations from the CPB, and on the floor
#' support about 0.90 to 0.94, with most misses on the upper side (see
#' [calibrate_alpha()]). It is calibrated by parametric bootstrap when the fit
#' went through [calibrate_alpha()].
#'
#' @param object A `"cpb"` object. Coefficient intervals need `se = "bootstrap"` at fit time; the
#'   interval for `"alpha"` does not.
#' @param parm Optional subset of parameters (coefficient names and/or `"alpha"`).
#' @param level Confidence level (default 0.95).
#' @param ... Unused.
#' @return A matrix of lower/upper bounds.
#' @method confint cpb
#' @export
confint.cpb <- function(object, parm, level = 0.95, ...) {
  a   <- (1 - level) / 2
  ok  <- .ud_draws(object, names(object$coefficients))
  coef_ci <- missing(parm) || !all(parm %in% "alpha")
  if (is.null(ok) && is.null(object$boot) && coef_ci)
    stop("Coefficient intervals need the bootstrap: refit with se = \"bootstrap\". The profile interval for alpha ",
         "needs none: confint(fit, \"alpha\").", call. = FALSE)
  if (is.null(ok) && coef_ci) .ud_boot_too_few()
  cib <- if (!is.null(ok)) t(apply(ok, 2, stats::quantile, c(a, 1 - a)))
         else if (coef_ci) matrix(NA_real_, length(object$coefficients), 2L,
                                  dimnames = list(names(object$coefficients), NULL))
  aci <- .cpb_alpha_ci(object, level = level)
  ci  <- rbind(cib, alpha = aci[c("lower", "upper")])
  colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
