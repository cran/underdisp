## ---------------------------------------------------------------------------
## Quantities of interest for the hurdle-CPB and ZI-CPB, and a boundary-corrected
## test for zero-inflation. The headline is the extensive/intensive decomposition
## of a first difference: how much of a covariate's effect runs through whether a
## unit participates versus how much it does.
## ---------------------------------------------------------------------------

## mean of the (possibly truncated) CPB at a single lambda
.cpb_mean1 <- function(lambda, alpha, truncated = FALSE) {
  kmax <- floor(lambda / (1 - alpha)) + 1L
  pm <- .cpb_pmf1(lambda, alpha, kmax, truncated = truncated)
  sum((0:kmax) * pm)
}

## --- implied_ceiling ------------------------------------------------------

#' @method implied_ceiling hurdle_cpb
#' @export
implied_ceiling.hurdle_cpb <- function(object, newdata = NULL, level = 0.95, ...)
  implied_ceiling(object$intensity, newdata = newdata, level = level)

#' @method implied_ceiling zi_cpb
#' @export
implied_ceiling.zi_cpb <- function(object, newdata = NULL, ...) {
  lam <- predict(object, newdata = newdata, type = "intensity")
  data.frame(lambda = lam, ceiling = lam / (1 - object$alpha), row.names = NULL)
}

## --- irr / odds ratios ----------------------------------------------------

#' @method irr hurdle_cpb
#' @export
irr.hurdle_cpb <- function(object, level = 0.95, ...) {
  intensity <- if (!is.null(object$intensity$boot)) irr(object$intensity, level = level)
               else .ud_irr_wald(object$intensity$coefficients, NULL, level, "count", "IRR")
  rbind(intensity, .or_rows(object$participation, level))
}

#' @method irr zi_cpb
#' @export
irr.zi_cpb <- function(object, level = 0.95, ...)
  rbind(.ud_irr_wald(object$coefficients, NULL, level, "count", "IRR"),
        .ud_irr_wald(object$zero_coef,   NULL, level, "binary", "OR"))

## --- first-difference decompositions --------------------------------------

.hurdle_decomp <- function(fit, variable, from, to, stage = "both") {
  a <- fit$intensity$alpha
  inx <- variable %in% names(fit$int_xref); inz <- variable %in% names(fit$part_xref)
  if (!inx && !inz) stop("'variable' is in neither the intensity nor the participation equation.")
  build <- function(v) {
    xr <- fit$int_xref; zr <- fit$part_xref
    if (inx && stage %in% c("both", "intensity", "count"))        xr[variable] <- v
    if (inz && stage %in% c("both", "participation", "zero"))     zr[variable] <- v
    lam <- exp(fit$int_fe_ref + sum(xr * fit$int_beta))
    p   <- stats::plogis(sum(zr * fit$part_beta))
    ey1 <- .cpb_mean1(lam, a, truncated = TRUE)
    c(p, ey1, p * ey1)
  }
  .ud_fd(component = c("participation", "intensity", "marginal"),
         from = build(from), to = build(to))
}

.zi_decomp <- function(fit, variable, from, to, stage = "both") {
  a <- fit$alpha
  inx <- variable %in% names(fit$int_xref); inz <- variable %in% names(fit$zero_xref)
  if (!inx && !inz) stop("'variable' is in neither the count nor the inflation equation.")
  build <- function(v) {
    xr <- fit$int_xref; zr <- fit$zero_xref
    if (inx && stage %in% c("both", "intensity", "count"))    xr[variable] <- v
    if (inz && stage %in% c("both", "participation", "zero")) zr[variable] <- v
    lam <- exp(fit$int_fe_ref + sum(xr * fit$int_beta))
    pi_ <- stats::plogis(sum(zr * fit$zero_beta))
    eyc <- .cpb_mean1(lam, a, truncated = FALSE)
    c(pi_, eyc, (1 - pi_) * eyc)
  }
  .ud_fd(component = c("zero", "count", "marginal"),
         from = build(from), to = build(to))
}

.boot_diff <- function(refit, decomp, data, variable, from, to, B, level, stage = "both") {
  n <- nrow(data); a <- (1 - level) / 2
  D <- vapply(seq_len(B), function(b) {
    f2 <- tryCatch(refit(data[sample(n, replace = TRUE), , drop = FALSE]), error = function(e) NULL)
    if (is.null(f2)) rep(NA_real_, 3) else decomp(f2, variable, from, to, stage)$diff
  }, numeric(3))
  list(lower = apply(D, 1, stats::quantile, a, na.rm = TRUE),
       upper = apply(D, 1, stats::quantile, 1 - a, na.rm = TRUE))
}

#' Extensive/intensive first-difference decomposition for a hurdle- or ZI-CPB
#'
#' The effect of moving `variable` `from` one value `to` another, decomposed into
#' the change in the participation (extensive) margin, the intensity (intensive)
#' margin, and the marginal expected count, holding other covariates at their
#' means. This is the quantity the single-equation defaults cannot separate.
#'
#' @param object A `"hurdle_cpb"` or `"zi_cpb"` object.
#' @param variable Name of a model-matrix column to vary (in either margin).
#' @param from,to The two values to contrast.
#' @param B Bootstrap replicates for percentile intervals; `0` (default) returns
#'   point estimates only. For `zi_cpb` the bootstrap is slow (each replicate
#'   refits the mixture).
#' @param level Confidence level.
#' @param data The original data frame; required when `B > 0`.
#' @param stage Which equation(s) the covariate moves in: `"both"` (default),
#'   the binary stage only (`"participation"`/`"zero"`), or the count stage
#'   only (`"intensity"`/`"count"`), holding the covariate at its reference in
#'   the other equation. All three component rows are always returned.
#' @param ... Unused; unknown arguments error.
#' @return A `"ud_fd"` data frame (the package-wide first-difference contract):
#'   columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`,
#'   one row per margin level.
#' @method first_difference hurdle_cpb
#' @export
first_difference.hurdle_cpb <- function(object, variable, from, to, B = 0, level = 0.95,
                                        data = NULL, stage = c("both", "participation",
                                        "intensity", "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  pt <- .hurdle_decomp(object, variable, from, to, stage)
  if (B > 0) {
    if (is.null(data)) stop("Bootstrap (B > 0) requires the original 'data'.")
    refit <- function(d) hurdle_cpb(object$formula, data = d,
                                    participation = object$part_formula, se = "none")
    ci <- .boot_diff(refit, .hurdle_decomp, data, variable, from, to, B, level, stage)
    pt$lower <- ci$lower; pt$upper <- ci$upper
    pt$method <- sprintf("bootstrap (refit, B=%d)", B)
  }
  pt
}

#' @rdname first_difference.hurdle_cpb
#' @method first_difference zi_cpb
#' @export
first_difference.zi_cpb <- function(object, variable, from, to, B = 0, level = 0.95,
                                    data = NULL, stage = c("both", "participation",
                                    "intensity", "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  pt <- .zi_decomp(object, variable, from, to, stage)
  if (B > 0) {
    if (is.null(data)) stop("Bootstrap (B > 0) requires the original 'data'.")
    refit <- function(d) zi_cpb(object$formula, data = d, zero = object$zero_formula)
    ci <- .boot_diff(refit, .zi_decomp, data, variable, from, to, B, level, stage)
    pt$lower <- ci$lower; pt$upper <- ci$upper
    pt$method <- sprintf("bootstrap (refit, B=%d)", B)
  }
  pt
}

## --- boundary-corrected test for zero-inflation ---------------------------

#' Boundary-corrected test for zero-inflation (CPB vs ZI-CPB)
#'
#' A likelihood-ratio test of whether a count needs a structural-zero component.
#' The (untruncated) CPB is nested in the zero-inflated CPB at a structural-zero
#' probability of zero. Because that value lies on the boundary of the parameter
#' space, the LR statistic follows a \eqn{\tfrac12\chi^2_0 + \tfrac12\chi^2_1}
#' mixture (Self and Liang 1987), which halves the naive \eqn{\chi^2_1} p-value.
#' The test uses an intercept-only structural-zero probability, so it is a clean
#' single-parameter boundary test; a covariate-dependent structural-zero model is
#' better compared with information criteria and proper scores via
#' [compare_dispersion()].
#'
#' @param object Either a model formula --- then `data` is required and the
#'   nested pair is fit internally --- or a fitted model (`cpb`, `cpb_fe`, or
#'   `zi_cpb`).
#' @param object2 The second argument: the data frame when `object` is a formula
#'   (so `zi_test(y ~ x, mydata)` works), or the second fitted model when
#'   `object` is a fit. In the model interface exactly one of the two models must
#'   be a zero-inflated `zi_cpb` and the other its non-inflated nest (`cpb` or
#'   `cpb_fe`); the two may carry any combination of fixed effects and
#'   robust/clustered standard errors.
#' @param data A data frame; an alternative to passing it as `object2`.
#' @param ... Passed to [cpb()] and [zi_cpb()] in the formula interface.
#' @return An object of class `"zi_test"` with the two log-likelihoods, the LR
#'   statistic, and the boundary-corrected p-value.
#' @examples
#' \donttest{
#' set.seed(1); x <- rnorm(500)
#' y <- rzicpb(500, lambda = exp(1.2 + 0.4 * x), alpha = 0.5, pi = 0.3)
#' zi_test(y ~ x, data = data.frame(y = y, x = x))
#' }
#' @export
zi_test <- function(object, object2 = NULL, data = NULL, ...) {
  if (inherits(object, "formula")) {
    dd <- if (!is.null(data)) data else object2
    if (is.null(dd) || !is.data.frame(dd))
      stop("zi_test(formula, data): supply the data frame as the second argument or via 'data ='.")
    cpb_fit <- cpb(object, data = dd, truncated = FALSE, se = "none", ...)
    zi_fit  <- zi_cpb(object, data = dd, zero = ~ 1, ...)
  } else {
    if (is.null(object2))
      stop("zi_test(model1, model2): supply two fitted models; one must be a zero-inflated 'zi_cpb'.")
    models <- list(object, object2)
    is_zi  <- vapply(models, function(m) inherits(m, c("zi_cpb", "zi_gec")), logical(1))
    if (!any(is_zi))
      stop("zi_test compares a zero-inflated model against its non-inflated nest, but neither ",
           "argument is zero-inflated. At least one model must be a 'zi_cpb' or 'zi_gec'.")
    if (all(is_zi))
      stop("Both arguments are zero-inflated; zi_test needs the non-inflated nest ",
           "(a 'cpb'/'cpb_fe' for zi_cpb, or 'gec'/'gec_fe' for zi_gec) as the other model.")
    zi_fit  <- models[[which(is_zi)]]
    cpb_fit <- models[[which(!is_zi)]]
    if (inherits(zi_fit, "zi_gec")) {
      if (!inherits(cpb_fit, c("gec", "gec_fe")))
        stop("The non-inflated model must be a 'gec' or 'gec_fe' (the nest of the zero-inflated GEC); got '",
             class(cpb_fit)[1L], "'.")
    } else if (!inherits(cpb_fit, c("cpb", "cpb_fe"))) {
      stop("The non-inflated model must be a 'cpb' or 'cpb_fe' (the nest of the zero-inflated CPB); got '",
           class(cpb_fit)[1L], "'.")
    }
    if (!is.null(cpb_fit$n) && !is.null(zi_fit$n) && cpb_fit$n != zi_fit$n)
      warning("The two models were fit on different numbers of observations (", cpb_fit$n,
              " vs ", zi_fit$n, "); the likelihood-ratio test assumes the same data.")
    if (length(zi_fit$zero_coef) > 1L)
      warning("The zero-inflated model has a covariate-dependent structural-zero probability; the ",
              "boundary p-value assumes intercept-only inflation. For covariate-dependent inflation, ",
              "compare by information criteria via compare_dispersion().")
  }
  LR <- max(2 * (zi_fit$loglik - cpb_fit$loglik), 0)
  p  <- 0.5 * stats::pchisq(LR, 1, lower.tail = FALSE)
  fam <- if (inherits(zi_fit, "zi_gec")) "GEC" else "CPB"
  structure(list(loglik_cpb = cpb_fit$loglik, loglik_zicpb = zi_fit$loglik,
                 LR = LR, p.value = p, family = fam), class = "zi_test")
}

#' @method print zi_test
#' @export
print.zi_test <- function(x, ...) {
  fam <- if (is.null(x$family)) "CPB" else x$family
  cat(sprintf("Boundary-corrected LR test for zero-inflation (%s vs ZI-%s)\n", fam, fam))
  cat(sprintf("  logLik: %s = %.2f, ZI-%s = %.2f\n", fam, x$loglik_cpb, fam, x$loglik_zicpb))
  cat(sprintf("  LR = %.2f,  p = %.4g   (0.5 chi^2_0 + 0.5 chi^2_1 mixture)\n", x$LR, x$p.value))
  cat(if (x$p.value < 0.05) "  Reject: a structural-zero component improves fit.\n"
      else "  No evidence of structural zero-inflation beyond the CPB.\n")
  invisible(x)
}
