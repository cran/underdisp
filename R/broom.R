## ---------------------------------------------------------------------------
## broom (tidyverse) integration for "cpb" objects. Methods are registered for
## broom only when it is installed (it is in Suggests); no hard dependency.
## Marginal quantities of interest are provided natively by predict(),
## first_difference(), implied_ceiling(), and irr(); marginaleffects gates on a
## fixed class whitelist and does not yet recognize "cpb".
## ---------------------------------------------------------------------------

#' Tidy a CPB fit (broom method)
#'
#' @param x A `"cpb"` object.
#' @param conf.int If `TRUE`, add `conf.low`/`conf.high` (requires a fit with
#'   bootstrap standard errors).
#' @param conf.level Confidence level for the interval.
#' @param ... Unused.
#' @return A data frame with one row per coefficient.
#' @exportS3Method broom::tidy
tidy.cpb <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  est <- x$coefficients
  se  <- if (is.null(x$se.beta)) rep(NA_real_, length(est)) else x$se.beta
  z   <- est / se
  out <- data.frame(term = names(est), estimate = as.numeric(est),
                    std.error = as.numeric(se), statistic = as.numeric(z),
                    p.value = as.numeric(2 * stats::pnorm(-abs(z))),
                    row.names = NULL, stringsAsFactors = FALSE)
  if (conf.int && !is.null(x$vcov)) {
    ci <- stats::confint(x, level = conf.level)
    out$conf.low  <- ci[names(est), 1]
    out$conf.high <- ci[names(est), 2]
  }
  out
}

#' Glance at a CPB fit (broom method)
#'
#' @param x A `"cpb"` object.
#' @param ... Unused.
#' @return A one-row data frame of fit statistics.
#' @exportS3Method broom::glance
glance.cpb <- function(x, ...) {
  data.frame(alpha = x$alpha, logLik = x$loglik,
             AIC = -2 * x$loglik + 2 * x$df, df = x$df, nobs = x$n,
             row.names = NULL, stringsAsFactors = FALSE)
}

#' Augment data with CPB fitted values (broom method)
#'
#' @param x A `"cpb"` object.
#' @param ... Unused.
#' @return A data frame with fitted values, response residuals, and the implied
#'   per-observation ceiling.
#' @exportS3Method broom::augment
augment.cpb <- function(x, ...) {
  data.frame(.fitted = x$fitted.values,
             .resid = x$Y - x$fitted.values,
             .ceiling = x$ceiling,
             row.names = NULL)
}

.zdf <- function(est, se) {
  z <- est / se
  data.frame(estimate = as.numeric(est), std.error = as.numeric(se),
             statistic = as.numeric(z), p.value = as.numeric(2 * stats::pnorm(-abs(z))),
             row.names = NULL, stringsAsFactors = FALSE)
}

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.cpb_fe <- function(x, ...) {
  se <- if (is.null(x$se.beta)) rep(NA_real_, length(x$coefficients)) else x$se.beta
  cbind(term = names(x$coefficients), .zdf(x$coefficients, se))
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.cpb_fe <- function(x, ...)
  data.frame(alpha = x$alpha, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, n.units = x$n_units, row.names = NULL)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_cpb <- function(x, ...) {
  ps  <- summary(x$participation)$coefficients
  part <- data.frame(component = "participation", term = rownames(ps),
                     estimate = ps[, 1], std.error = ps[, 2], statistic = ps[, 3],
                     p.value = ps[, 4], row.names = NULL, stringsAsFactors = FALSE)
  ise <- if (is.null(x$intensity$se.beta)) rep(NA_real_, length(x$intensity$coefficients))
         else x$intensity$se.beta
  int <- cbind(component = "intensity", term = names(x$intensity$coefficients),
               .zdf(x$intensity$coefficients, ise))
  rbind(part, int)
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_cpb <- function(x, ...) {
  idf <- if (!is.null(x$intensity$df)) x$intensity$df else length(x$intensity$coefficients) + 1L
  df  <- length(stats::coef(x$participation)) + idf
  ll  <- as.numeric(stats::logLik(x$participation)) + x$intensity$loglik
  data.frame(alpha = x$intensity$alpha, logLik = ll, AIC = -2 * ll + 2 * df,
             df = df, nobs = x$n, row.names = NULL)
}

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_cpb <- function(x, ...) {
  cse <- if (is.null(x$se.beta)) rep(NA_real_, length(x$coefficients)) else x$se.beta
  zse <- if (is.null(x$se.zero)) rep(NA_real_, length(x$zero_coef))    else x$se.zero
  rbind(cbind(component = "count", term = names(x$coefficients), .zdf(x$coefficients, cse)),
        cbind(component = "zero",  term = names(x$zero_coef),    .zdf(x$zero_coef, zse)))
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_cpb <- function(x, ...)
  data.frame(alpha = x$alpha, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.gec <- function(x, ...) {
  se <- if (is.null(x$se.beta)) rep(NA_real_, length(x$coefficients)) else x$se.beta
  cbind(term = names(x$coefficients), .zdf(x$coefficients, se))
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.gec <- function(x, ...)
  data.frame(dispersion = x$delta, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_gec <- function(x, ...) {
  ps  <- summary(x$participation)$coefficients
  part <- data.frame(component = "participation", term = rownames(ps),
                     estimate = ps[, 1], std.error = ps[, 2], statistic = ps[, 3],
                     p.value = ps[, 4], row.names = NULL, stringsAsFactors = FALSE)
  ise <- if (is.null(x$intensity$se.beta)) rep(NA_real_, length(x$intensity$coefficients))
         else x$intensity$se.beta
  int <- cbind(component = "intensity", term = names(x$intensity$coefficients),
               .zdf(x$intensity$coefficients, ise))
  rbind(part, int)
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_gec <- function(x, ...) {
  df <- length(stats::coef(x$participation)) + x$intensity$df
  ll <- as.numeric(stats::logLik(x$participation)) + x$intensity$loglik
  data.frame(dispersion = x$delta, logLik = ll, AIC = -2 * ll + 2 * df, df = df, nobs = x$n,
             row.names = NULL)
}

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_gec <- function(x, ...) {
  cse <- if (is.null(x$se.beta)) rep(NA_real_, length(x$coefficients)) else x$se.beta
  zse <- if (is.null(x$se.zero)) rep(NA_real_, length(x$zero_coef))    else x$se.zero
  rbind(cbind(component = "count", term = names(x$coefficients), .zdf(x$coefficients, cse)),
        cbind(component = "zero",  term = names(x$zero_coef),    .zdf(x$zero_coef, zse)))
}

#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_gec <- function(x, ...)
  data.frame(dispersion = x$delta, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL)

## broom methods for the matched Poisson/NB families (count_reg / zi_count / hurdle_count)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.count_reg <- function(x, ...) {
  est <- x$coefficients; se <- if (is.null(x$se.beta)) rep(NA_real_, length(est)) else x$se.beta
  cbind(term = names(est), .zdf(est, se))
}
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.count_reg <- function(x, ...)
  data.frame(family = x$family, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL, stringsAsFactors = FALSE)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_count <- function(x, ...) {
  cse <- if (is.null(x$se.beta)) rep(NA_real_, length(x$coefficients)) else x$se.beta
  rbind(cbind(component = "count", term = names(x$coefficients), .zdf(x$coefficients, cse)),
        cbind(component = "zero",  term = names(x$zero.coefficients), .zdf(x$zero.coefficients, NA_real_)))
}
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_count <- function(x, ...)
  data.frame(family = x$family, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL, stringsAsFactors = FALSE)

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_count <- function(x, ...) {
  ps <- summary(x$participation)$coefficients
  part <- data.frame(component = "participation", term = rownames(ps),
                     estimate = ps[, 1], std.error = ps[, 2], statistic = ps[, 3],
                     p.value = ps[, 4], row.names = NULL, stringsAsFactors = FALSE)
  ise <- if (is.null(x$intensity$se.beta)) rep(NA_real_, length(x$intensity$coefficients)) else x$intensity$se.beta
  int <- cbind(component = "intensity", term = names(x$intensity$coefficients),
               .zdf(x$intensity$coefficients, ise))
  rbind(part, int)
}
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_count <- function(x, ...)
  data.frame(family = x$family, logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df,
             df = x$df, nobs = x$n, row.names = NULL, stringsAsFactors = FALSE)
