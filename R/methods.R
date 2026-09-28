## ---------------------------------------------------------------------------
## S3 methods for class "cpb"
## ---------------------------------------------------------------------------

#' @export
print.cpb <- function(x, ...) {
  cat("Continuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)" else "(untruncated)", "\n")
  cat("Call:  ", deparse(x$call), "\n\n", sep = "")
  cat("Coefficients:\n"); print(round(x$coefficients, 4))
  cat("\nalpha (shape parameter):", round(x$alpha, 4),
      "  median implied bound:", round(stats::median(x$ceiling), 2), "\n")
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
    v <- .cpb_moments(.cpb_rate(object), object$alpha, isTRUE(object$truncated))$var
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
                           isTRUE(object$support_binding))
  if (is.finite(lr$p))
    lr$note <- c(lr$note, paste0("the p-value is asymptotic, which over-rejects in finite samples at this boundary; ",
                                 "dispersion_test() gives the parametric-bootstrap p-value"))
  out <- list(call = object$call, truncated = object$truncated, coefficients = ctab,
              alpha = object$alpha, alpha.ci = aci, ceiling = object$ceiling,
              loglik = object$loglik, aic = -2 * object$loglik + 2 * object$df,
              LR = lr$LR, LR.p = lr$p, LR.note = lr$note,
              support_binding = isTRUE(object$support_binding), converged = object$converged,
              n = object$n, nobs = nobs(object), se.type = object$se.type,
              nboot_ok = if (!is.null(object$boot)) attr(object$boot, "nboot_ok") else NA)
  class(out) <- "summary.cpb"
  out
}

#' @export
print.summary.cpb <- function(x, ...) {
  cat("\nContinuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)\n" else "(untruncated)\n")
  cat("N =", x$n, if (!is.null(x$nobs) && x$nobs != x$n) paste0(" (weight total ", format(x$nobs), ")"),
      "   inference:", x$se.type, "\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  calibrated <- isTRUE(grepl("calibrated", attr(x$alpha.ci, "method"), fixed = TRUE))
  cat("\nalpha =", round(x$alpha, 4),
      sprintf("  (%s profile %d%% CI: %.3f to %.3f%s)\n", if (calibrated) "calibrated" else "first-order", 95L,
              x$alpha.ci["lower"], x$alpha.ci["upper"],
              if (x$alpha.ci["boundary"] == 1) ", lower limit at the parameter bound" else ""))
  cat("Implied ceiling lambda/(1-alpha): median", round(stats::median(x$ceiling), 2),
      "  range", paste(round(range(x$ceiling), 2), collapse = " to "), "\n")
  cat("logLik =", round(x$loglik, 2), "   AIC =", round(x$aic, 2), "\n")
  if (!is.na(x$LR))
    cat(sprintf("LR vs %sPoisson (H0: alpha = 1): %.2f, p %s\n",
                if (x$truncated) "ZT-" else "", x$LR, format.pval(x$LR.p)))
  for (nt in x$LR.note) cat("Note: ", nt, ".\n", sep = "")
  if (isTRUE(x$support_binding))
    cat("Note: the fitted ceiling reaches max.support; alpha is bounded by the guard, not the data.\n")
  if (!is.na(x$nboot_ok))
    cat("(", x$nboot_ok, " bootstrap resamples converged)\n", sep = "")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' Confidence intervals for a CPB fit
#'
#' Coefficient intervals use the cold-multistart bootstrap percentile method
#' (validated to nominal coverage). The interval for `alpha` is a
#' profile-likelihood interval: first-order by default (it covers about 0.90 to
#' 0.94 in simulations from the CPB, with its misses on the upper side; see
#' [calibrate_alpha()] for the reason), and calibrated by parametric bootstrap
#' when the fit went through [calibrate_alpha()].
#'
#' @param object A `"cpb"` object fit with `se = "bootstrap"`.
#' @param parm Optional subset of parameters (coefficient names and/or `"alpha"`).
#' @param level Confidence level (default 0.95).
#' @param ... Unused.
#' @return A matrix of lower/upper bounds.
#' @method confint cpb
#' @export
confint.cpb <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$boot))
    stop("Bootstrap required; refit with se = \"bootstrap\".")
  a  <- (1 - level) / 2
  ok <- object$boot[complete.cases(object$boot), , drop = FALSE]
  cib <- t(apply(ok[, 1:object$p, drop = FALSE], 2, quantile, c(a, 1 - a)))
  rownames(cib) <- names(object$coefficients)
  aci <- .cpb_alpha_ci(object, level = level)
  ci  <- rbind(cib, alpha = aci[c("lower", "upper")])
  colnames(ci) <- c(paste0(format(100 * a), "%"), paste0(format(100 * (1 - a)), "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
