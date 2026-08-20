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
  invisible(x)
}

#' @method coef cpb
#' @export
coef.cpb <- function(object, ...) object$coefficients

#' @method vcov cpb
#' @export
vcov.cpb <- function(object, ...) {
  if (is.null(object$vcov))
    stop("No covariance available; refit with se = \"bootstrap\".")
  object$vcov
}

#' @method logLik cpb
#' @export
logLik.cpb <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")

#' @method nobs cpb
#' @export
nobs.cpb <- function(object, ...) object$n

#' @method fitted cpb
#' @export
fitted.cpb <- function(object, ...) object$fitted.values

#' @method residuals cpb
#' @export
residuals.cpb <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  r <- object$Y - object$fitted.values
  if (type == "pearson") r <- r / sqrt(object$alpha * object$fitted.values)
  r
}

#' Summarize a CPB fit
#'
#' @param object A `"cpb"` object.
#' @param ... Unused.
#' @return An object of class `"summary.cpb"` with the coefficient table, the
#'   dispersion parameter and its profile-likelihood interval, the implied ceiling,
#'   fit statistics, and the likelihood-ratio test against a (zero-truncated) Poisson.
#' @method summary cpb
#' @export
summary.cpb <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * pnorm(-abs(z)))
  aci <- .cpb_alpha_profile_ci(object)
  LR  <- if (is.na(object$loglik.null)) NA_real_ else -2 * (object$loglik.null - object$loglik)
  out <- list(call = object$call, truncated = object$truncated, coefficients = ctab,
              alpha = object$alpha, alpha.ci = aci, ceiling = object$ceiling,
              loglik = object$loglik, aic = -2 * object$loglik + 2 * object$df,
              # alpha = 1 is on the parameter boundary: null is 0.5*chi2_0 + 0.5*chi2_1 (Self & Liang 1987)
              LR = LR, LR.p = if (is.na(LR)) NA_real_ else 0.5 * pchisq(LR, 1, lower.tail = FALSE),
              n = object$n, se.type = object$se.type,
              nboot_ok = if (!is.null(object$boot)) attr(object$boot, "nboot_ok") else NA)
  class(out) <- "summary.cpb"
  out
}

#' @export
print.summary.cpb <- function(x, ...) {
  cat("\nContinuous Parameter Binomial regression",
      if (x$truncated) "(zero-truncated)\n" else "(untruncated)\n")
  cat("N =", x$n, "   inference:", x$se.type, "\n\n")
  printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  cat("\nalpha =", round(x$alpha, 4),
      sprintf("  (profile %d%% CI: %.3f to %.3f%s)\n", 95L, x$alpha.ci["lower"],
              x$alpha.ci["upper"], if (x$alpha.ci["boundary"] == 1) ", at feasibility boundary" else ""))
  cat("Implied ceiling lambda/(1-alpha): median", round(stats::median(x$ceiling), 2),
      "  range", paste(round(range(x$ceiling), 2), collapse = " to "), "\n")
  cat("logLik =", round(x$loglik, 2), "   AIC =", round(x$aic, 2), "\n")
  if (!is.na(x$LR))
    cat(sprintf("LR vs %sPoisson (H0: alpha = 1): %.2f, p %s\n",
                if (x$truncated) "ZT-" else "", x$LR, format.pval(x$LR.p)))
  if (!is.na(x$nboot_ok))
    cat("(", x$nboot_ok, " bootstrap resamples converged)\n", sep = "")
  invisible(x)
}

#' Confidence intervals for a CPB fit
#'
#' Coefficient intervals use the cold-multistart bootstrap percentile method
#' (validated to nominal coverage); the interval for `alpha` uses the
#' profile-likelihood method, which is reliable except under strong underdispersion,
#' where `alpha` sits at the feasibility boundary and the interval is one-sided.
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
  aci <- .cpb_alpha_profile_ci(object, level = level)
  ci  <- rbind(cib, alpha = aci[c("lower", "upper")])
  colnames(ci) <- c(paste0(format(100 * a), "%"), paste0(format(100 * (1 - a)), "%"))
  if (!missing(parm)) ci <- ci[parm, , drop = FALSE]
  ci
}
