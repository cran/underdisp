## ---------------------------------------------------------------------------
## Cross-family parity methods: make logLik()/nobs()/coef()/AIC() work uniformly
## across the CPB, GEC, and count families, and give the GEC family the same
## rate-ratio and first-difference quantities of interest as the others.
## ---------------------------------------------------------------------------

## --- logLik / nobs / coef for the CPB classes that lacked them -------------

#' @method logLik cpb_fe
#' @export
logLik.cpb_fe <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
#' @method nobs cpb_fe
#' @export
nobs.cpb_fe <- function(object, ...) object$n
#' @method coef cpb_fe
#' @export
coef.cpb_fe <- function(object, ...) object$coefficients

#' @method logLik hurdle_cpb
#' @export
logLik.hurdle_cpb <- function(object, ...) {
  idf <- if (!is.null(object$intensity$df)) object$intensity$df
         else length(object$intensity$coefficients) + 1L
  df <- length(stats::coef(object$participation)) + idf
  ll <- as.numeric(stats::logLik(object$participation)) + object$intensity$loglik
  structure(ll, df = df, nobs = object$n, class = "logLik")
}
#' @method nobs hurdle_cpb
#' @export
nobs.hurdle_cpb <- function(object, ...) object$n

#' @method logLik zi_cpb
#' @export
logLik.zi_cpb <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
#' @method nobs zi_cpb
#' @export
nobs.zi_cpb <- function(object, ...) object$n
## coef for the ZI classes: return BOTH blocks as a named list (count + zero),
## symmetric with coef.hurdle_* which returns list(participation, intensity).
.zi_coef <- function(object)
  list(count = object$coefficients,
       zero  = if (!is.null(object$zero.coefficients)) object$zero.coefficients else object$zero_coef)
#' @method coef zi_cpb
#' @export
coef.zi_cpb <- function(object, ...) .zi_coef(object)

## --- coef for the hurdle classes (which store split coefficients) ----------
## Returns a named list with the participation and intensity coefficient blocks,
## since a hurdle model has two coefficient vectors.
.hurdle_coef <- function(object)
  list(participation = stats::coef(object$participation),
       intensity = object$intensity$coefficients)
#' @method coef hurdle_cpb
#' @export
coef.hurdle_cpb <- function(object, ...) .hurdle_coef(object)
#' @method coef hurdle_count
#' @export
coef.hurdle_count <- function(object, ...) .hurdle_coef(object)
#' @method coef hurdle_gec
#' @export
coef.hurdle_gec <- function(object, ...) .hurdle_coef(object)

## --- residuals parity -------------------------------------------------------
## Response residuals (y - fitted) for every class that lacked them; two-part
## fits use the marginal fitted mean, matching fitted()/predict(type="response").
#' @method residuals gec
#' @export
residuals.gec <- function(object, ...) object$Y - object$fitted.values
#' @method residuals cpb_fe
#' @export
residuals.cpb_fe <- function(object, ...) object$Y - object$fitted.values
#' @method residuals hurdle_cpb
#' @export
residuals.hurdle_cpb <- function(object, ...) object$y - fitted(object)
#' @method residuals hurdle_count
#' @export
residuals.hurdle_count <- function(object, ...) object$y - fitted(object)
#' @method residuals hurdle_gec
#' @export
residuals.hurdle_gec <- function(object, ...) object$y - fitted(object)
#' @method residuals zi_cpb
#' @export
residuals.zi_cpb <- function(object, ...) object$y - fitted(object)
#' @method residuals zi_count
#' @export
residuals.zi_count <- function(object, ...) object$y - fitted(object)
#' @method residuals zi_gec
#' @export
residuals.zi_gec <- function(object, ...) object$y - fitted(object)

## --- real summaries ---------------------------------------------------------
## Every model class gets a genuine summary(): per-equation coefficient tables
## with z and p wherever the fit carries inference (bootstrap SEs, an analytic
## vcov, or an embedded glm), a dispersion line, and an explicit note when no
## inference was requested at fit time.
.sum_block <- function(b, se = NULL) {
  if (is.null(b) || !length(b)) return(NULL)
  if (is.null(se) || !all(is.finite(se))) {
    out <- cbind(Estimate = unname(b)); rownames(out) <- names(b)
    return(out)
  }
  zv <- unname(b) / se
  out <- cbind(Estimate = unname(b), `Std. Error` = se, `z value` = zv,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(zv)))
  rownames(out) <- names(b); out
}
.glm_block <- function(glmfit) {
  sm <- summary(glmfit)$coefficients
  keep <- !grepl("^factor\\(", rownames(sm))
  sm[keep, , drop = FALSE]
}
.ud_summary <- function(object, blocks, disp = NULL, note = NULL)
  structure(list(object = object, blocks = blocks, disp = disp, note = note),
            class = "summary.underdisp")
#' @method print summary.underdisp
#' @export
print.summary.underdisp <- function(x, ...) {
  if (is.null(x$blocks)) { print(x$object); return(invisible(x)) }
  cat(class(x$object)[1], "fit\n")
  if (!is.null(x$object$call)) cat("Call:  ", deparse(x$object$call), "\n", sep = "")
  for (nm in names(x$blocks)) {
    bl <- x$blocks[[nm]]
    if (is.null(bl)) next
    cat("\n", nm, ":\n", sep = "")
    if (ncol(bl) >= 4) stats::printCoefmat(bl, signif.stars = FALSE)
    else print(round(bl, 4))
  }
  if (!is.null(x$disp)) cat("\n", x$disp, "\n", sep = "")
  if (!is.null(x$note)) cat(x$note, "\n")
  invisible(x)
}
.no_inf_note <- function(has_se)
  if (!has_se) "Note: no inference was requested at fit time (se = \"none\")." else NULL
#' @method summary cpb_fe
#' @export
summary.cpb_fe <- function(object, ...) {
  se <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta)))
          object$se.beta[names(object$coefficients)] else NULL
  .ud_summary(object,
    stats::setNames(list(.sum_block(object$coefficients, se)),
                    if (is.null(se)) "Coefficients" else "Coefficients (unit-bootstrap SEs)"),
    disp = sprintf("alpha (shape) = %.4f%s", object$alpha,
                   if (identical(object$bias_correct, "jackknife"))
                     " [split-panel jackknife corrected]" else ""),
    note = .no_inf_note(!is.null(se)))
}
#' @method summary gec
#' @export
summary.gec <- function(object, ...) {
  se <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  .ud_summary(object,
    stats::setNames(list(.sum_block(object$coefficients, se)),
                    if (is.null(se)) "Coefficients" else "Coefficients (bootstrap SEs)"),
    disp = sprintf("dispersion delta (Var/Mean) = %.4f%s", object$delta,
                   if (!is.null(object$se.delta) && is.finite(object$se.delta))
                     sprintf("  (SE %.4f)", object$se.delta) else ""),
    note = .no_inf_note(!is.null(se)))
}
#' @method summary gec_fe
#' @export
summary.gec_fe <- function(object, ...) {
  se <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta)))
          object$se.beta[names(object$coefficients)] else NULL
  .ud_summary(object,
    stats::setNames(list(.sum_block(object$coefficients, se)),
                    if (is.null(se)) "Coefficients" else "Coefficients (unit-bootstrap SEs)"),
    disp = sprintf("dispersion delta (Var/Mean) = %.4f%s", object$delta,
                   if (identical(object$bias_correct, "jackknife"))
                     " [split-panel jackknife corrected]" else ""),
    note = .no_inf_note(!is.null(se)))
}
#' @method summary hurdle_cpb
#' @export
summary.hurdle_cpb <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .glm_block(object$participation),
         "Intensity (zero-truncated CPB)" = .sum_block(object$intensity$coefficients, ise)),
    disp = sprintf("intensity alpha (shape) = %.4f", object$intensity$alpha),
    note = .no_inf_note(!is.null(ise)))
}
#' @method summary hurdle_count
#' @export
summary.hurdle_count <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .glm_block(object$participation),
         "Intensity (zero-truncated)" = .sum_block(object$intensity$coefficients, ise)),
    disp = .shape_line(object$family, object$theta),
    note = .no_inf_note(!is.null(ise)))
}
#' @method summary hurdle_gec
#' @export
summary.hurdle_gec <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .glm_block(object$participation),
         "Intensity (zero-truncated GEC)" = .sum_block(object$intensity$coefficients, ise)),
    disp = sprintf("intensity dispersion delta = %.4f", object$delta),
    note = .no_inf_note(!is.null(ise)))
}
#' @method summary zi_cpb
#' @export
summary.zi_cpb <- function(object, ...) {
  cse <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  zse <- if (!is.null(object$se.zero) && any(is.finite(object$se.zero))) object$se.zero else NULL
  .ud_summary(object,
    list("Count component (CPB)" = .sum_block(object$coefficients, cse),
         "Zero-inflation (logit)" = .sum_block(object$zero_coef, zse)),
    disp = sprintf("intensity alpha (shape) = %.4f", object$alpha),
    note = .no_inf_note(!is.null(cse)))
}
#' @method summary zi_count
#' @export
summary.zi_count <- function(object, ...) {
  cse <- if (!is.null(object$vcov)) sqrt(pmax(diag(object$vcov), 0)) else NULL
  .ud_summary(object,
    list("Count component" = .sum_block(object$coefficients, cse),
         "Zero-inflation" = .sum_block(object$zero.coefficients, NULL)),
    disp = .shape_line(object$family, object$theta),
    note = .no_inf_note(!is.null(cse)))
}
#' @method summary zi_gec
#' @export
summary.zi_gec <- function(object, ...) {
  cse <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  zse <- if (!is.null(object$se.zero) && any(is.finite(object$se.zero))) object$se.zero else NULL
  .ud_summary(object,
    list("Count component (GEC)" = .sum_block(object$coefficients, cse),
         "Zero-inflation (logit)" = .sum_block(object$zero_coef, zse)),
    disp = sprintf("count dispersion delta = %.4f", object$delta),
    note = .no_inf_note(!is.null(cse)))
}

## --- rate ratios and first differences for the GEC and FE-CPB families -----
## The GEC/CPB mean is a monotone function of the rate lambda = exp(offset + fe +
## x'beta): rate ratios are exp(beta); the first difference is the change in the
## fitted mean (gec_mean_cpp for GEC, matching predict(type="response"); lambda
## for the CPB, matching predict.cpb). SEs use the S3 vcov() (bootstrap), not a
## stored $vcov slot -- gec/gec_fe expose covariance only through the method.

.safe_vcov <- function(object) tryCatch({ v <- stats::vcov(object); if (is.matrix(v)) v else NULL },
                                        error = function(e) NULL)

.irr_from_vcov <- function(b, vc, level, method = "bootstrap (normal)") {
  .ud_irr_wald(b, if (is.null(vc)) NULL else sqrt(diag(vc)), level, "count", "IRR",
               method = method)
}
#' @rdname irr.count
#' @method irr gec
#' @export
irr.gec <- function(object, level = 0.95, ...)
  .irr_from_vcov(object$coefficients, .safe_vcov(object), level)
#' @rdname irr.count
#' @method irr gec_fe
#' @export
irr.gec_fe <- function(object, level = 0.95, ...)
  .irr_from_vcov(object$coefficients, .safe_vcov(object), level)
#' @rdname irr.count
#' @method irr cpb_fe
#' @export
irr.cpb_fe <- function(object, level = 0.95, ...)
  .irr_from_vcov(object$coefficients, vcov(object), level)

## first difference on the mean scale via a supplied mean function meanfun(eta),
## eta = base + x'beta (base = 0, or the FE reference for the concentrated
## models). Returns the unified ud_fd single row.
.rate_fd <- function(b, xref, base, meanfun, vc, variable, from, to, level) {
  if (!variable %in% names(b)) stop("'variable' is not a covariate in the model.")
  mk <- function(v) { x <- xref; x[variable] <- v; x }
  mfun <- function(bb, v) meanfun(base + sum(mk(v) * bb))
  ci <- .fd_delta_ci(function(bb) mfun(bb, to) - mfun(bb, from), b, vc, level)
  .ud_fd(component = "mean", from = mfun(b, from), to = mfun(b, to),
         lower = ci[1], upper = ci[2], method = if (is.null(vc)) "none" else "delta")
}
.gec_meanfun <- function(object) function(eta) gec_mean_cpp(exp(eta), object$delta, object$max.support)
#' @rdname first_difference.count
#' @method first_difference gec
#' @export
first_difference.gec <- function(object, variable, from, to, level = 0.95, ...) {
  .fd_dots(...)
  .rate_fd(object$coefficients, colMeans(object$X), 0, .gec_meanfun(object),
           .safe_vcov(object), variable, from, to, level)
}
#' @rdname first_difference.count
#' @method first_difference gec_fe
#' @export
first_difference.gec_fe <- function(object, variable, from, to, level = 0.95, ...) {
  .fd_dots(...)
  .rate_fd(object$coefficients, object$xref, object$fe_ref, .gec_meanfun(object),
           .safe_vcov(object), variable, from, to, level)
}
#' Bootstrap covariance for a fixed-effects CPB fit
#'
#' The covariance of the covariate coefficients from the stored pairs/cluster
#' bootstrap (`se = "bootstrap"` at fit time); `NULL` when no bootstrap was run.
#' @param object A `"cpb_fe"` object.
#' @param ... Unused.
#' @return A covariance matrix, or `NULL`.
#' @method vcov cpb_fe
#' @export
vcov.cpb_fe <- function(object, ...) {
  if (is.null(object$boot)) return(NULL)
  ok <- object$boot[stats::complete.cases(object$boot), , drop = FALSE]
  if (nrow(ok) < 2L) return(NULL)
  v <- stats::cov(ok)
  dimnames(v) <- list(names(object$coefficients), names(object$coefficients))
  v
}
#' @rdname first_difference.count
#' @method first_difference cpb_fe
#' @export
first_difference.cpb_fe <- function(object, variable, from, to, level = 0.95, ...) {
  .fd_dots(...)
  .rate_fd(object$coefficients, object$xref, object$fe_ref, exp,   # CPB reports the rate lambda
           vcov(object), variable, from, to, level)
}

## --- rate/odds ratios + first differences for the two-part GEC models ------
## Mirror the CPB/count hurdle and ZI QoIs (irr.hurdle_cpb / first_difference.*).
#' @rdname irr.count
#' @method irr hurdle_gec
#' @export
irr.hurdle_gec <- function(object, level = 0.95, ...)
  rbind(.irr_from_vcov(object$int_beta, .safe_vcov(object$intensity), level),
        .or_rows(object$participation, level))
#' @rdname irr.count
#' @method irr zi_gec
#' @export
irr.zi_gec <- function(object, level = 0.95, ...)
  rbind(.ud_irr_wald(object$coefficients,
                     if (!is.null(object$se.beta) && all(is.finite(object$se.beta)))
                       object$se.beta else NULL,
                     level, "count", "IRR", method = "bootstrap (normal)"),
        .ud_irr_wald(object$zero_coef,
                     if (!is.null(object$se.zero) && all(is.finite(object$se.zero)))
                       object$se.zero else NULL,
                     level, "binary", "OR", method = "bootstrap (normal)"))
#' @rdname first_difference.count
#' @method first_difference hurdle_gec
#' @export
first_difference.hurdle_gec <- function(object, variable, from, to, level = 0.95,
                                        stage = c("both", "participation", "intensity",
                                                  "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  linkinv <- stats::make.link(if (is.null(object$link)) "logit" else object$link)$linkinv
  inx <- variable %in% names(object$int_xref); inz <- variable %in% names(object$part_xref)
  if (!inx && !inz) stop("'variable' is in neither the intensity nor the participation equation.")
  zr0 <- object$part_xref; xr0 <- object$int_xref
  pfun <- function(g, v) { zr <- zr0
    if (inz && stage %in% c("both", "participation", "zero")) zr[variable] <- v
    linkinv(sum(zr * g)) }
  mfun <- function(b, v) { xr <- xr0
    if (inx && stage %in% c("both", "intensity", "count")) xr[variable] <- v
    mu <- exp(sum(xr * b))
    p0 <- gec_lp0_cpp(c(b, log(object$delta)), matrix(xr, 1), 0L, 0, object$max.support)[1, 2]
    gec_mean_cpp(mu, object$delta, object$max.support) / (1 - p0) }           # E(Y | Y > 0)
  bvc <- tryCatch(stats::vcov(object$participation), error = function(e) NULL)
  .twopart_fd(pfun, mfun, object$part_beta, object$int_beta, bvc, NULL,
              from, to, level, c("participation", "intensity", "marginal"))
}
#' @rdname first_difference.count
#' @method first_difference zi_gec
#' @export
first_difference.zi_gec <- function(object, variable, from, to, level = 0.95,
                                    stage = c("both", "participation", "intensity",
                                              "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  xref <- colMeans(object$X); zref <- colMeans(object$Z)
  inx <- variable %in% names(xref); inz <- variable %in% names(zref)
  if (!inx && !inz) stop("'variable' is in neither the count nor the inflation equation.")
  ## rows: P(structural zero) level, count-component mean, marginal (1-pi)*m
  pfun <- function(g, v) { zr <- zref
    if (inz && stage %in% c("both", "participation", "zero")) zr[variable] <- v
    stats::plogis(sum(zr * g)) }
  mfun <- function(b, v) { xr <- xref
    if (inx && stage %in% c("both", "intensity", "count")) xr[variable] <- v
    gec_mean_cpp(exp(sum(xr * b)), object$delta, object$max.support) }
  pi_f <- pfun(object$zero_coef, from); pi_t <- pfun(object$zero_coef, to)
  m_f  <- mfun(object$coefficients, from); m_t <- mfun(object$coefficients, to)
  .ud_fd(component = c("zero", "count", "marginal"),
         from  = c(pi_f, m_f, (1 - pi_f) * m_f),
         to    = c(pi_t, m_t, (1 - pi_t) * m_t))
}
