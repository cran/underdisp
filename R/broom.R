## ---------------------------------------------------------------------------
## broom (tidyverse) integration. Methods are registered for broom only when it
## is installed (it is in Suggests); no hard dependency. Every class has the same
## tidy() and glance(): tidy() reports each coefficient's estimate, standard
## error, z and p (the separated terms of a zero-inflation or participation
## equation without a standard error, as summary() shows them) and, with
## conf.int = TRUE, the interval confint() gives; glance() reports the fit
## statistics and the dispersion parameter under one name, with its direction.
## Marginal quantities of interest are provided natively by predict(),
## first_difference(), implied_ceiling(), and irr(); marginaleffects gates on a
## fixed class whitelist and does not yet recognize these classes.
## ---------------------------------------------------------------------------

.zdf <- function(est, se) {
  z <- est / se
  data.frame(estimate = as.numeric(est), std.error = as.numeric(se),
             statistic = as.numeric(z), p.value = as.numeric(2 * stats::pnorm(-abs(z))),
             row.names = NULL, stringsAsFactors = FALSE)
}

## One equation's rows: `prefix` names the component of a two-part fit (its terms are prefixed by it, so that
## modelsummary lays the table out without a `shape` argument, and a `component` column stays for filtering);
## `sep` are separated terms, without a standard error
.tidy_rows <- function(est, se, prefix = NULL, sep = character(0)) {
  se <- if (is.null(se)) rep(NA_real_, length(est)) else unname(as.numeric(se))
  se[names(est) %in% sep] <- NA_real_
  out <- cbind(term = names(est), .zdf(est, se), stringsAsFactors = FALSE)
  if (!is.null(prefix)) out <- cbind(component = prefix, out[, c("term", "estimate", "std.error", "statistic", "p.value")],
                                     stringsAsFactors = FALSE)
  out
}
.tidy_finish <- function(x, out, conf.int, conf.level) {
  cs <- .coef_scale(x)$scale
  sc <- if (is.null(out$component)) rep(cs, nrow(out))
        else ifelse(out$component %in% c("participation", "zero"), if (is.null(x$link)) "logit" else x$link, cs)
  if (!is.null(out$component)) out$term <- paste(out$component, out$term, sep = ":")
  if (conf.int) {
    ci <- tryCatch(stats::confint(x, level = conf.level), error = function(e) NULL)
    k <- if (is.null(ci)) rep(NA_integer_, nrow(out)) else match(out$term, rownames(ci))
    out$conf.low  <- if (is.null(ci)) NA_real_ else unname(ci[k, 1])
    out$conf.high <- if (is.null(ci)) NA_real_ else unname(ci[k, 2])
  }
  out$scale <- sc
  rownames(out) <- NULL
  out
}

#' Tidy a fitted model (broom method)
#'
#' One row per coefficient with its estimate, standard error, z statistic and
#' p-value, for every model class of the package. The two-part classes
#' (`hurdle_*`, `zi_*`) have a `component` column (the \pkg{broom} convention
#' for multi-equation models) and the term prefixed by its component, so that
#' \pkg{modelsummary} lays the table out without a `shape` argument. A term of
#' a zero-inflation or participation equation that the data do not bound
#' (separation; see [zi_cpb()]) has no standard error.
#'
#' The `scale` column names what each coefficient multiplies: `"log mean"`
#' for the mean-parameterized count families (`"log untruncated mean"` for a
#' zero-truncated fit or a hurdle's intensity), the binary equation's link for a
#' participation or zero-inflation row, and otherwise the family's own
#' parameter (`"log rate"` for the CPB and GEC, `"log lambda"` for the
#' rate-parameterized COM-Poisson, `"log centrality"` for the EWP, and so on;
#' see [irr.count]). Coefficients on different scales are not comparable
#' across the columns of a table; [ame()] puts every model's effect on the
#' mean.
#'
#' @param x A fitted model from this package.
#' @param conf.int If `TRUE`, add `conf.low`/`conf.high`, the intervals of
#'   [confint()] (`NA` where the fit carries no inference).
#' @param conf.level Confidence level for the interval.
#' @param ... Unused.
#' @return A data frame with one row per coefficient.
#' @exportS3Method broom::tidy
tidy.cpb <- function(x, conf.int = FALSE, conf.level = 0.95, ...)
  .tidy_finish(x, .tidy_rows(x$coefficients, x$se.beta), conf.int, conf.level)

## the participation equation's rows of a hurdle: its covariance (model-based, robust or clustered)
.tidy_part <- function(x) {
  V <- .part_vcov(x); b <- x$part_beta
  .tidy_rows(b, if (is.null(V)) NULL else sqrt(pmax(diag(V), 0))[names(b)], "participation", x$part_separation$diverging)
}

#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.cpb_fe <- function(x, conf.int = FALSE, conf.level = 0.95, ...)
  .tidy_finish(x, .tidy_rows(x$coefficients, x$se.beta[names(x$coefficients)]), conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.gec <- function(x, conf.int = FALSE, conf.level = 0.95, ...)
  .tidy_finish(x, .tidy_rows(x$coefficients, x$se.beta[names(x$coefficients)]), conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.count_reg <- function(x, conf.int = FALSE, conf.level = 0.95, ...)
  .tidy_finish(x, .tidy_rows(x$coefficients, x$se.beta), conf.int, conf.level)
.tidy_hurdle <- function(x, conf.int, conf.level)
  .tidy_finish(x, rbind(.tidy_part(x), .tidy_rows(x$intensity$coefficients, x$intensity$se.beta, "intensity")),
               conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_cpb <- function(x, conf.int = FALSE, conf.level = 0.95, ...) .tidy_hurdle(x, conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_gec <- function(x, conf.int = FALSE, conf.level = 0.95, ...) .tidy_hurdle(x, conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.hurdle_count <- function(x, conf.int = FALSE, conf.level = 0.95, ...) .tidy_hurdle(x, conf.int, conf.level)
.tidy_zi <- function(x, g, conf.int, conf.level)
  .tidy_finish(x, rbind(.tidy_rows(x$coefficients, x$se.beta, "count"),
                        .tidy_rows(g, x$se.zero, "zero", x$separation$diverging)), conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_cpb <- function(x, conf.int = FALSE, conf.level = 0.95, ...) .tidy_zi(x, x$zero_coef, conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_gec <- function(x, conf.int = FALSE, conf.level = 0.95, ...) .tidy_zi(x, x$zero_coef, conf.int, conf.level)
#' @rdname tidy.cpb
#' @exportS3Method broom::tidy
tidy.zi_count <- function(x, conf.int = FALSE, conf.level = 0.95, ...)
  .tidy_zi(x, x$zero.coefficients, conf.int, conf.level)

## ---- the dispersion parameter under one name, with its direction -------------
## The families' dispersion parameters point different ways (the CPB's alpha below 1 is underdispersion, the
## gamma-count's alpha above 1 is), so every glance() and texreg table names the parameter and the side of its
## Poisson value that is underdispersion.
.disp_rules <- list(cpb = c("alpha", "<", 1), gec = c("delta", "<", 1), negbin = c("theta", "over", NA),
                    compois = c("nu", ">", 1), mpcmp = c("nu", ">", 1), genpois = c("lambda", "<", 0),
                    gammacount = c("alpha", ">", 1), doublepois = c("theta", ">", 1), hyperpois = c("gamma", "<", 1),
                    ewp2 = c("beta", ">", 0), ewp3 = c("beta1, beta2", ">", 0))
.disp_info <- function(x) {
  key <- if (inherits(x, c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb"))) "cpb"
         else if (inherits(x, c("gec", "hurdle_gec", "zi_gec"))) "gec" else x$family
  val <- switch(key, cpb = if (!is.null(x$alpha)) x$alpha else x$intensity$alpha,
                gec = x$delta, x$theta)
  if (identical(key, "poisson")) return(list(name = NA_character_, value = NA_real_, direction = "equidispersed", rule = ""))
  r <- .disp_rules[[key]]
  dir <- if (r[2] == "over") (if (is.finite(val) && 1 / val > 1e-6) "overdispersed" else "equidispersed")
         else {
           thr <- as.numeric(r[3]); u <- if (r[2] == "<") val < thr - 1e-3 else val > thr + 1e-3
           o <- if (r[2] == "<") val > thr + 1e-3 else val < thr - 1e-3
           if (all(u)) "underdispersed" else if (all(o)) "overdispersed" else if (any(u) && any(o)) "mixed" else "equidispersed"
         }
  if (key == "cpb" && isTRUE(val > 0.999)) dir <- "equidispersed"           # the CPB's Poisson limit
  list(name = r[1], value = val, direction = dir,
       rule = if (r[2] == "over") "" else paste0(r[2], " ", r[3], " under"))
}

#' Glance at a fitted model (broom method)
#'
#' One row of fit statistics for every model class of the package:
#' `logLik`, `AIC`, `BIC`, `df` and `nobs` (the weight total of a weighted fit),
#' the dispersion parameter as `dispersion` with its name
#' (`dispersion_parameter`) and the direction its estimate points
#' (`dispersion_direction`; the families' parameters point different ways, so
#' the CPB's `alpha` below 1 and the gamma-count's `alpha` above 1 are both
#' underdispersion), `coef_scale`, the scale of the count equation's
#' coefficients (see [tidy.cpb()]), and the class's own columns (`alpha` for
#' the CPB family, `family` for the matched count families, `n.units` for the
#' fixed-effects fits).
#'
#' @param x A fitted model from this package.
#' @param ... Unused.
#' @return A one-row data frame of fit statistics.
#' @exportS3Method broom::glance
glance.cpb <- function(x, ...) .glance(x, alpha = x$alpha)
.glance <- function(x, ...) {
  ll <- as.numeric(stats::logLik(x)); df <- attr(stats::logLik(x), "df"); n <- stats::nobs(x)
  di <- .disp_info(x)
  out <- data.frame(..., logLik = ll, AIC = -2 * ll + 2 * df, BIC = -2 * ll + df * log(n), df = df, nobs = n,
                    dispersion = di$value[1L], dispersion_parameter = di$name, dispersion_direction = di$direction,
                    coef_scale = .coef_scale(x)$scale, row.names = NULL, stringsAsFactors = FALSE)
  if (length(di$value) > 1L) out$dispersion2 <- di$value[2L]               # the EWP's second tail
  out
}
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.cpb_fe <- function(x, ...) .glance(x, alpha = x$alpha, n.units = x$n_units)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_cpb <- function(x, ...) .glance(x, alpha = x$intensity$alpha)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_cpb <- function(x, ...) .glance(x, alpha = x$alpha)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.gec <- function(x, ...)
  if (inherits(x, "gec_fe")) .glance(x, n.units = x$n_units) else .glance(x)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_gec <- function(x, ...) .glance(x)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_gec <- function(x, ...) .glance(x)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.count_reg <- function(x, ...) .glance(x, family = x$family)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.zi_count <- function(x, ...) .glance(x, family = x$family)
#' @rdname glance.cpb
#' @exportS3Method broom::glance
glance.hurdle_count <- function(x, ...) .glance(x, family = x$family)

#' Augment data with CPB fitted values (broom method)
#'
#' @param x A `"cpb"` object.
#' @param data Optional: the data frame the model was fit to. Its rows the fit
#'   used (those without a missing model variable) are returned with the new
#'   columns appended.
#' @param ... Unused.
#' @return A data frame with fitted values, response residuals, and the implied
#'   per-observation ceiling (appended to `data` when it is given).
#' @exportS3Method broom::augment
augment.cpb <- function(x, data = NULL, ...) {
  out <- data.frame(.fitted = x$fitted.values, .resid = x$Y - x$fitted.values, .ceiling = x$ceiling, row.names = NULL)
  if (is.null(data)) return(out)
  rows <- if (!is.null(x$rows)) x$rows else if (nrow(data) == length(x$Y)) seq_len(nrow(data))
          else stop("'data' must be the data frame the model was fit to.", call. = FALSE)
  if (max(rows) > nrow(data)) stop("'data' must be the data frame the model was fit to.", call. = FALSE)
  cbind(data[rows, , drop = FALSE], out)
}
