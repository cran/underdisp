## ---------------------------------------------------------------------------
## Cross-family parity methods: make logLik()/nobs()/coef()/AIC() work uniformly
## across the CPB, GEC, and count families, and give the GEC family the same
## rate-ratio and first-difference quantities of interest as the others.
## ---------------------------------------------------------------------------

## --- logLik / nobs / coef for the CPB classes that lacked them -------------

#' @method logLik cpb_fe
#' @export
logLik.cpb_fe <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")
#' @method nobs cpb_fe
#' @export
nobs.cpb_fe <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted
#' @method coef cpb_fe
#' @export
coef.cpb_fe <- function(object, ...) object$coefficients

#' @method logLik hurdle_cpb
#' @export
logLik.hurdle_cpb <- function(object, ...) {
  ll <- as.numeric(stats::logLik(object$participation)) + object$intensity$loglik
  structure(ll, df = .ud_twopart_df(object), nobs = nobs(object), class = "logLik")
}
#' @method nobs hurdle_cpb
#' @export
nobs.hurdle_cpb <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted

#' @method logLik zi_cpb
#' @export
logLik.zi_cpb <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")
#' @method nobs zi_cpb
#' @export
nobs.zi_cpb <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted
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

## --- covariance parity ---------------------------------------------------------
## Every class answers vcov(): the bootstrap covariance where draws were stored,
## the block-diagonal of the two margins for the hurdles (which factorize), and
## NULL, never an error, when no inference was requested at fit time.
.block_diag <- function(A, B) {
  if (is.null(A) || is.null(B)) return(NULL)
  V <- matrix(0, nrow(A) + nrow(B), ncol(A) + ncol(B))
  V[seq_len(nrow(A)), seq_len(ncol(A))] <- A
  V[nrow(A) + seq_len(nrow(B)), ncol(A) + seq_len(ncol(B))] <- B
  dimnames(V) <- list(c(paste0("participation:", rownames(A)), paste0("intensity:", rownames(B))),
                      c(paste0("participation:", colnames(A)), paste0("intensity:", colnames(B))))
  V
}
## the participation block (model-based, robust, or clustered; see participation.R) and the intensity's, NA
## where the intensity has no inference
.hurdle_vcov <- function(object) {
  pv <- .part_vcov(object)
  iv <- tryCatch(stats::vcov(object$intensity), error = function(e) NULL)
  if (is.null(iv)) { nm <- names(object$int_beta); iv <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm)) }
  .block_diag(pv, iv)
}
#' @method vcov hurdle_cpb
#' @export
vcov.hurdle_cpb <- function(object, ...) .hurdle_vcov(object)
#' @method vcov hurdle_gec
#' @export
vcov.hurdle_gec <- function(object, ...) .hurdle_vcov(object)
#' @method vcov hurdle_count
#' @export
vcov.hurdle_count <- function(object, ...) .hurdle_vcov(object)
## the zero-inflated bootstraps' covariance of the count and zero-equation coefficients (the zero_fe unit dummies,
## which are not the same units across replicates, have no draws)
.zi_boot_vcov <- function(object) {
  D <- .ud_draws(object)
  if (is.null(D)) return(NULL)
  cols <- c(paste0("count:", names(object$coefficients)), paste0("zero:", names(object$zero_coef)))
  stats::cov(D[, intersect(cols, colnames(D)), drop = FALSE])
}
#' @method vcov zi_cpb
#' @export
vcov.zi_cpb <- function(object, ...) .zi_boot_vcov(object)
#' @method vcov zi_gec
#' @export
vcov.zi_gec <- function(object, ...) .zi_boot_vcov(object)

## --- confidence-interval parity ----------------------------------------------
## One rule per inference type: bootstrap percentile where draws are stored
## (cpb, gec), the normal approximation for the fixed-effects and mixture
## bootstraps (their draws share the incidental-parameters bias or only the
## standard errors are kept), and Wald intervals from the analytic covariance
## for the matched count families. A fit without inference errors informatively.
.ci_wald <- function(est, se, level) {
  a <- (1 - level) / 2; z <- stats::qnorm(1 - a)
  ci <- cbind(est - z * se, est + z * se)
  dimnames(ci) <- list(names(est), paste0(format(100 * c(a, 1 - a), trim = TRUE), " %"))
  ci
}
.ci_pick <- function(ci, parm) if (missing(parm) || is.null(parm)) ci else ci[parm, , drop = FALSE]
.no_inference <- function(how = "se = \"bootstrap\"") stop("No inference available; refit with ", how, ".", call. = FALSE)

#' Confidence intervals for every model class
#'
#' `confint()` methods for the whole family, one interval rule per inference
#' type: bootstrap percentile intervals for [cpb()] (coefficients; the
#' dispersion parameter gets its profile-likelihood interval, first-order or
#' calibrated by [calibrate_alpha()]) and [gec()]
#' (coefficients and `delta`); normal-approximation intervals from the
#' bootstrap standard errors for the fixed-effects fits ([cpb_fe()],
#' [gec_fe()]) and the zero-inflated mixtures ([zi_cpb()], [zi_gec()]); Wald
#' intervals from the analytic covariance for [count_reg()] and [zi_count()]
#' (the dispersion parameter's interval is mapped from its estimation scale to
#' the natural scale); and, for the hurdles, the participation model's Wald
#' intervals (prefixed `participation:`) stacked over the intensity model's intervals
#' (prefixed `intensity:`); the zero-inflated classes prefix `count:` and `zero:`.
#' The same names label `vcov()`, `tidy()`, and the texreg tables. A fit without inference errors informatively.
#'
#' @param object A fitted model from this package.
#' @param parm Optional subset of parameter names.
#' @param level Confidence level (default 0.95).
#' @param ... Unused.
#' @return A two-column matrix of lower and upper bounds.
#' @name confint.underdisp
#' @seealso [confint.cpb()]
NULL

#' @rdname confint.underdisp
#' @method confint gec
#' @export
confint.gec <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$boot)) .no_inference("se = \"bootstrap\"")
  a <- (1 - level) / 2; fe <- inherits(object, "gec_fe")
  est <- c(object$coefficients, if (!fe) c(delta = object$delta))
  ok <- .ud_draws(object, names(est))
  if (is.null(ok)) {                                      # fewer than two converged replicates
    .ud_boot_too_few()
    ci <- .ci_wald(est, rep(NA_real_, length(est)), level)
  } else if (fe) {                                        # FE bootstrap: normal approximation
    ci <- .ci_wald(est, apply(ok, 2, stats::sd), level)
  } else {
    ci <- t(apply(ok, 2, stats::quantile, c(a, 1 - a)))
    colnames(ci) <- paste0(format(100 * c(a, 1 - a), trim = TRUE), " %")
  }
  .ci_pick(ci, parm)
}
#' @rdname confint.underdisp
#' @method confint cpb_fe
#' @export
confint.cpb_fe <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$se.beta) || !any(is.finite(object$se.beta))) .no_inference()
  .ci_pick(.ci_wald(object$coefficients, object$se.beta[names(object$coefficients)], level), parm)
}
#' @rdname confint.underdisp
#' @method confint count_reg
#' @export
confint.count_reg <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$vcov_full)) .no_inference("se = \"analytic\", \"robust\" or \"cluster\"")
  fam <- .count_fam(object$family); p <- length(object$coefficients)
  est <- object$coefficients; se <- sqrt(pmax(diag(object$vcov_full)[seq_len(p)], 0))
  ci <- .ci_wald(est, se, level)
  if (fam$nshape) {                                        # dispersion: interval on the link scale, mapped back
    th <- fam$shape_link(object$theta); sth <- sqrt(pmax(diag(object$vcov_full)[p + seq_len(fam$nshape)], 0))
    cth <- fam$shape_inv(.ci_wald(stats::setNames(th, fam$shape_name), sth, level))
    rownames(cth) <- sub("^(log|atanh)\\((.*)\\)$", "\\2", fam$shape_name)
    ci <- rbind(ci, cth)
  }
  .ci_pick(ci, parm)
}
#' @rdname confint.underdisp
#' @method confint zi_count
#' @export
confint.zi_count <- function(object, parm, level = 0.95, ...) {
  if (is.null(object$vcov_full)) .no_inference("se = \"analytic\", \"robust\" or \"cluster\"")
  fam <- .count_fam(object$family); pc <- length(object$coefficients); pz <- length(object$zero.coefficients)
  est <- c(stats::setNames(object$coefficients, paste0("count:", names(object$coefficients))),
           stats::setNames(object$zero.coefficients, paste0("zero:", names(object$zero.coefficients))))
  se <- sqrt(pmax(diag(object$vcov_full)[seq_len(pc + pz)], 0))
  ci <- .ci_wald(est, se, level)
  if (fam$nshape) {
    th <- fam$shape_link(object$theta); sth <- sqrt(pmax(diag(object$vcov_full)[pc + pz + seq_len(fam$nshape)], 0))
    cth <- fam$shape_inv(.ci_wald(stats::setNames(th, fam$shape_name), sth, level))
    rownames(cth) <- sub("^(log|atanh)\\((.*)\\)$", "\\2", fam$shape_name)
    ci <- rbind(ci, cth)
  }
  .ci_pick(ci, parm)
}
.hurdle_confint <- function(object, parm, level) {
  pc <- object$part_beta; pv <- .part_vcov(object)
  pci <- .ci_wald(pc, sqrt(pmax(diag(pv), 0))[names(pc)], level); rownames(pci) <- paste0("participation:", names(pc))
  pci[rownames(pci) %in% paste0("participation:", object$part_separation$diverging), ] <- NA   # separated
  ici <- tryCatch(stats::confint(object$intensity, level = level), error = function(e) NULL)
  if (is.null(ici))                                        # no intensity inference was requested
    ici <- matrix(NA_real_, length(object$int_beta), 2L, dimnames = list(names(object$int_beta), NULL))
  rownames(ici) <- paste0("intensity:", rownames(ici)); colnames(ici) <- colnames(pci)
  .ci_pick(rbind(pci, ici), parm)
}
#' @rdname confint.underdisp
#' @method confint hurdle_cpb
#' @export
confint.hurdle_cpb <- function(object, parm, level = 0.95, ...) .hurdle_confint(object, parm, level)
#' @rdname confint.underdisp
#' @method confint hurdle_gec
#' @export
confint.hurdle_gec <- function(object, parm, level = 0.95, ...) .hurdle_confint(object, parm, level)
#' @rdname confint.underdisp
#' @method confint hurdle_count
#' @export
confint.hurdle_count <- function(object, parm, level = 0.95, ...) .hurdle_confint(object, parm, level)
.zi_confint <- function(object, parm, level) {
  if (is.null(object$se.beta) || !any(is.finite(object$se.beta))) .no_inference()
  est <- c(stats::setNames(object$coefficients, paste0("count:", names(object$coefficients))),
           stats::setNames(object$zero_coef, paste0("zero:", names(object$zero_coef))))
  se  <- c(object$se.beta, object$se.zero)
  .ci_pick(.ci_wald(est, se, level), parm)
}
#' @rdname confint.underdisp
#' @method confint zi_cpb
#' @export
confint.zi_cpb <- function(object, parm, level = 0.95, ...) .zi_confint(object, parm, level)
#' @rdname confint.underdisp
#' @method confint zi_gec
#' @export
confint.zi_gec <- function(object, parm, level = 0.95, ...) .zi_confint(object, parm, level)

## --- residuals parity -------------------------------------------------------
## Response residuals (y - fitted) for every class that lacked them; two-part
## fits use the marginal fitted mean, matching fitted()/predict(type="response").
## response or Pearson residuals of a single-equation fit, the variance the exact one of the fitted pmf (zero-truncated
## when the fit is), as dispersion_profile() computes it
.ud_residuals <- function(object, type) {
  r <- object$Y - stats::fitted(object)
  if (type == "pearson") r <- r / sqrt(pmax(.disp_condition(.disp_moments(object))$var, 1e-12))
  r
}
#' Residuals of a fitted model
#'
#' Response residuals (`y` minus the fitted mean) for every model class, and for
#' the single-equation classes Pearson residuals, the response residual divided
#' by the standard deviation of the fitted distribution (its exact variance;
#' the conditional one of a zero-truncated fit).
#' @param object A fitted model from this package.
#' @param type `"response"` (default) or, for `cpb`, `cpb_fe`, `gec`, `gec_fe`
#'   and `count_reg`, `"pearson"`.
#' @param ... Unused.
#' @return A numeric vector, in the rows of the data.
#' @name residuals.underdisp
NULL
#' @rdname residuals.underdisp
#' @method residuals gec
#' @export
residuals.gec <- function(object, type = c("response", "pearson"), ...) .ud_residuals(object, match.arg(type))
#' @rdname residuals.underdisp
#' @method residuals cpb_fe
#' @export
residuals.cpb_fe <- function(object, type = c("response", "pearson"), ...) .ud_residuals(object, match.arg(type))
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
  if (is.null(se) || !any(is.finite(se))) {
    out <- cbind(Estimate = unname(b)); rownames(out) <- names(b)
    return(out)
  }
  se <- unname(se); se[!is.finite(se)] <- NA; zv <- unname(b) / se
  out <- cbind(Estimate = unname(b), `Std. Error` = se, `z value` = zv,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(zv)))
  rownames(out) <- names(b); out
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
  for (nt in x$note) cat(nt, "\n")
  invisible(x)
}
.no_inf_note <- function(has_se)
  if (!has_se) "Note: no inference was requested at fit time (se = \"none\")." else NULL
## the survivor count of a bootstrap, or the no-inference note
.boot_note <- function(object, has_se) {
  if (!has_se) return(.no_inf_note(FALSE))
  nb <- if (!is.null(object$boot)) attr(object$boot, "nboot_ok") else NULL
  if (is.null(nb)) NULL else sprintf("(%d bootstrap resamples converged)", nb)
}
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
    note = .boot_note(object, !is.null(se)))
}
#' @method summary gec
#' @export
summary.gec <- function(object, ...) {
  se <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  .ud_summary(object,
    stats::setNames(list(.sum_block(object$coefficients, se)),
                    if (is.null(se)) "Coefficients" else "Coefficients (bootstrap SEs)"),
    disp = sprintf("dispersion delta (Var/Mean on an unbounded support) = %.4f%s", object$delta,
                   if (!is.null(object$se.delta) && is.finite(object$se.delta))
                     sprintf("  (SE %.4f)", object$se.delta) else ""),
    note = .boot_note(object, !is.null(se)))
}
#' @method summary gec_fe
#' @export
summary.gec_fe <- function(object, ...) {
  se <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta)))
          object$se.beta[names(object$coefficients)] else NULL
  .ud_summary(object,
    stats::setNames(list(.sum_block(object$coefficients, se)),
                    if (is.null(se)) "Coefficients" else "Coefficients (unit-bootstrap SEs)"),
    disp = sprintf("dispersion delta (Var/Mean on an unbounded support) = %.4f%s", object$delta,
                   if (identical(object$bias_correct, "jackknife"))
                     " [split-panel jackknife corrected]" else ""),
    note = .boot_note(object, !is.null(se)))
}
#' @method summary hurdle_cpb
#' @export
summary.hurdle_cpb <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .part_block(object),
         "Intensity (zero-truncated CPB)" = .sum_block(object$intensity$coefficients, ise)),
    disp = sprintf("intensity alpha (shape) = %.4f", object$intensity$alpha),
    note = c(.no_inf_note(!is.null(ise)), .part_note(object)))
}
#' @method summary hurdle_count
#' @export
summary.hurdle_count <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .part_block(object),
         "Intensity (zero-truncated)" = .sum_block(object$intensity$coefficients, ise)),
    disp = .shape_line(object$family, object$theta),
    note = c(.no_inf_note(!is.null(ise)), .part_note(object)))
}
#' @method summary hurdle_gec
#' @export
summary.hurdle_gec <- function(object, ...) {
  ise <- if (!is.null(object$intensity$se.beta) && any(is.finite(object$intensity$se.beta)))
           object$intensity$se.beta else NULL
  .ud_summary(object,
    list("Participation (binary stage)" = .part_block(object),
         "Intensity (zero-truncated GEC)" = .sum_block(object$intensity$coefficients, ise)),
    disp = sprintf("intensity dispersion delta = %.4f", object$delta),
    note = c(.no_inf_note(!is.null(ise)), .part_note(object)))
}
#' @method summary zi_cpb
#' @export
summary.zi_cpb <- function(object, ...) {
  cse <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  zse <- if (!is.null(object$se.zero) && any(is.finite(object$se.zero))) object$se.zero else NULL
  .ud_summary(object,
    list("Count component (CPB)" = .sum_block(object$coefficients, cse),
         "Zero-inflation (logit)" = .sum_mask(.sum_block(object$zero_coef, zse), object$separation$diverging)),
    disp = sprintf("intensity alpha (shape) = %.4f", object$alpha),
    note = c(.no_inf_note(!is.null(cse)), .sep_note(object$separation, "zero-inflation"), .zi_zero_note(object)))
}
#' @method summary zi_count
#' @export
summary.zi_count <- function(object, ...) {
  cse <- if (!is.null(object$vcov)) object$se.beta else NULL
  .ud_summary(object,
    list("Count component" = .sum_block(object$coefficients, cse),
         "Zero-inflation" = .sum_mask(.sum_block(object$zero.coefficients, if (is.null(cse)) NULL else object$se.zero),
                                      object$separation$diverging)),
    disp = .shape_line(object$family, object$theta),
    note = c(.no_inf_note(!is.null(cse)), .sep_note(object$separation, "zero-inflation")))
}
#' @method summary zi_gec
#' @export
summary.zi_gec <- function(object, ...) {
  cse <- if (!is.null(object$se.beta) && any(is.finite(object$se.beta))) object$se.beta else NULL
  zse <- if (!is.null(object$se.zero) && any(is.finite(object$se.zero))) object$se.zero else NULL
  .ud_summary(object,
    list("Count component (GEC)" = .sum_block(object$coefficients, cse),
         "Zero-inflation (logit)" = .sum_mask(.sum_block(object$zero_coef, zse), object$separation$diverging)),
    disp = sprintf("count dispersion delta = %.4f", object$delta),
    note = c(.no_inf_note(!is.null(cse)), .sep_note(object$separation, "zero-inflation")))
}

## --- rate ratios and first differences for the GEC and FE-CPB families -----
## The GEC/CPB mean is a monotone function of the rate lambda = exp(offset + fe +
## x'beta): rate ratios are exp(beta); the first difference is the change in the
## fitted mean (gec_mean_cpp for GEC, .cpb_mean for the CPB, each matching
## predict(type = "response")), with the interval of .qoi_ci() (qoi_unified.R).

.safe_vcov <- function(object) tryCatch({ v <- stats::vcov(object); if (is.matrix(v)) v else NULL },
                                        error = function(e) NULL)

.irr_from_vcov <- function(b, vc, level, method = "bootstrap (normal)") {
  .ud_irr_wald(b, if (is.null(vc)) NULL else sqrt(diag(vc)), level, "count", "IRR",
               method = method)
}
#' @rdname irr.count
#' @method irr gec
#' @export
irr.gec <- function(object, level = 0.95, ...) .ud_irr_draws(object, object$coefficients, level)
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

## offset held at its mean in every first difference
.mean_offset <- function(object) if (is.null(object$offset)) 0 else mean(object$offset)
#' @rdname first_difference.count
#' @method first_difference gec
#' @export
first_difference.gec <- function(object, variable, from, to, level = 0.95, at = NULL, ...) {
  .fd_dots(...)
  .fd_build(.fd_spec(object, variable, from, to, at = at), level)
}
#' @rdname first_difference.count
#' @method first_difference gec_fe
#' @export
first_difference.gec_fe <- function(object, variable, from, to, level = 0.95, at = NULL, ...) {
  .fd_dots(...)
  .fd_build(.fd_spec(object, variable, from, to, at = at), level)
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
  ok <- .ud_draws(object, names(object$coefficients))
  if (is.null(ok)) NULL else stats::cov(ok)
}
#' @rdname first_difference.count
#' @method first_difference cpb_fe
#' @export
first_difference.cpb_fe <- function(object, variable, from, to, level = 0.95, at = NULL, ...) {
  .fd_dots(...)
  .fd_build(.fd_spec(object, variable, from, to, at = at), level)
}

## --- rate/odds ratios + first differences for the two-part GEC models ------
## Mirror the CPB/count hurdle and ZI QoIs (irr.hurdle_cpb / first_difference.*).
#' @rdname irr.count
#' @method irr hurdle_gec
#' @export
irr.hurdle_gec <- function(object, level = 0.95, ...)
  rbind(irr(object$intensity, level = level), .or_rows(object, level))
#' @rdname irr.count
#' @method irr zi_gec
#' @export
irr.zi_gec <- function(object, level = 0.95, ...) .zi_irr(object, level)
#' @rdname first_difference.count
#' @method first_difference hurdle_gec
#' @export
first_difference.hurdle_gec <- function(object, variable, from, to, level = 0.95,
                                        stage = c("both", "participation", "intensity",
                                                  "zero", "count"), at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  .fd_build(.fd_spec(object, variable, from, to, at = at, stage = stage), level)
}
#' @rdname first_difference.count
#' @method first_difference zi_gec
#' @export
first_difference.zi_gec <- function(object, variable, from, to, level = 0.95,
                                    stage = c("both", "participation", "intensity",
                                              "zero", "count"), at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  .fd_build(.fd_spec(object, variable, from, to, at = at, stage = stage), level)
}
