## ---------------------------------------------------------------------------
## Quantities of interest for the hurdle-CPB and ZI-CPB, and a boundary-corrected
## test for zero-inflation. The headline is the extensive/intensive decomposition
## of a first difference: how much of a covariate's effect runs through whether a
## unit participates versus how much it does.
## ---------------------------------------------------------------------------

## --- implied_ceiling ------------------------------------------------------

#' @method implied_ceiling hurdle_cpb
#' @export
implied_ceiling.hurdle_cpb <- function(object, newdata = NULL, level = 0.95, ...) {
  ## the intensity's ceiling applies to every observation (the stored intensity
  ## rate is full-length), not only to the positives the intensity was fit on
  int <- object$intensity
  lam <- if (is.null(newdata)) object$lambda_full else predict(int, newdata = newdata, type = "rate")
  out <- data.frame(lambda = lam, ceiling = lam / (1 - int$alpha),
                    max_count = .cpb_top(lam / (1 - int$alpha), .cpb_support(object)), row.names = NULL)
  if (inherits(int, "cpb")) {                         # a profile interval exists for the pooled intensity
    aci <- .cpb_alpha_ci(int, level = level)
    out$ceiling_at_alpha_lower <- lam / (1 - aci[["lower"]]); out$ceiling_at_alpha_upper <- lam / (1 - aci[["upper"]])
  }
  out
}

#' @method implied_ceiling zi_cpb
#' @export
implied_ceiling.zi_cpb <- function(object, newdata = NULL, ...) {
  lam <- predict(object, newdata = newdata, type = "rate")
  data.frame(lambda = lam, ceiling = lam / (1 - object$alpha),
             max_count = .cpb_top(lam / (1 - object$alpha), .cpb_support(object)), row.names = NULL)
}

## --- irr / odds ratios ----------------------------------------------------

#' @method irr hurdle_cpb
#' @export
irr.hurdle_cpb <- function(object, level = 0.95, ...)
  rbind(irr(object$intensity, level = level), .or_rows(object, level))

#' @method irr zi_cpb
#' @export
irr.zi_cpb <- function(object, level = 0.95, ...) .zi_irr(object, level)
## a zero-inflated CPB or GEC fit's ratios: bootstrap-normal, as confint() gives them, the zero_fe dummies left out
## and the separated inflation terms NA
.zi_irr <- function(object, level) {
  zfe <- if (!is.null(object$zero_fe)) object$zero_fe else object$refit_args$zero_fe
  g <- object$zero_coef; keep <- !.is_fe_col(names(g), zfe)
  rbind(.ud_irr_wald(object$coefficients, object$se.beta, level, "count", "IRR", method = "bootstrap (normal)"),
        .irr_mask(.ud_irr_wald(g[keep], if (!is.null(object$se.zero)) object$se.zero[names(g)][keep], level, "inflation",
                               "OR", method = "bootstrap (normal)"), object$separation$diverging, "inflation"))
}


## --- first-difference decompositions --------------------------------------

#' Extensive/intensive first-difference decomposition for a hurdle- or ZI-CPB
#'
#' The effect of moving `variable` `from` one value `to` another on the
#' participation (extensive) margin, the intensity (intensive) margin, and the
#' marginal expected count, holding other covariates at their means (or at the
#' values `at` gives them). The two margins are reported side by side; their
#' changes do not add up to the marginal change (see [first_difference.count]).
#' With `B = 0` the intervals follow the package's one rule (see
#' [first_difference()]): the stored bootstrap draws of a `zi_cpb` fit, and for
#' a hurdle the delta method over its two margins' covariances (the
#' participation model's, and the intensity's stored bootstrap draws). `B > 0`
#' refits the whole model on resampled data instead.
#'
#' @param object A `"hurdle_cpb"` or `"zi_cpb"` object.
#' @param variable Name of a model-matrix column to vary (in either margin).
#' @param from,to The two values to contrast.
#' @param B Refit-bootstrap replicates for percentile intervals; `0` (default)
#'   uses the fit's own inference instead (see Details). For `zi_cpb` the refit
#'   bootstrap is slow (each replicate refits the mixture).
#' @param level Confidence level.
#' @param data The original data frame; required when `B > 0`. Each replicate
#'   refits the model with the fit's own settings (fixed effects, link, offset,
#'   weights, support, and the like), resampling whole fixed-effects units (or
#'   the fit's clusters) when it has them; an offset or weights given to the fit
#'   as a vector must then be a column of `data` instead.
#' @param stage Which equation(s) the covariate moves in: `"both"` (default),
#'   the binary stage only (`"participation"`/`"zero"`), or the count stage
#'   only (`"intensity"`/`"count"`), holding the covariate at its reference in
#'   the other equation. All three component rows are always returned.
#' @param cores Worker processes for the bootstrap refits (default 1); see [cpb()].
#' @param at Optional named numeric vector fixing other model-matrix columns at
#'   given values instead of their means (see [first_difference()]).
#' @param ... Unused; unknown arguments error.
#' @return A `"ud_fd"` data frame (the package-wide first-difference contract):
#'   columns `component`, `from`, `to`, `diff`, `lower`, `upper`, `method`,
#'   one row per margin level.
#' @method first_difference hurdle_cpb
#' @export
first_difference.hurdle_cpb <- function(object, variable, from, to, B = 0, level = 0.95,
                                        data = NULL, stage = c("both", "participation",
                                        "intensity", "zero", "count"), cores = 1L, at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  spec_of <- function(f) .fd_spec(f, variable, from, to, at = at, stage = stage)
  pt <- .fd_build(spec_of(object), level)
  if (B > 0) {
    if (is.null(data)) stop("Bootstrap (B > 0) requires the original 'data'.")
    ra <- .fd_refit_args(object, c("fe", "part_fe", "link"))
    ra$dots$support <- .cpb_support(object)
    ## a refit replicate follows the fit's own jackknife outcome, as ame(B > 0) does: no correction when the gate
    ## refused it on the data (on King's support a replicate skips the gate)
    if (!is.null(ra$dots$bias_correct) && inherits(object$intensity, "cpb_fe"))
      ra$dots$bias_correct <- object$intensity$bias_correct
    refit <- function(d) do.call(hurdle_cpb, c(list(object$formula, data = d, participation = object$part_formula,
                                                    se = "none", fe = ra$fe, part_fe = ra$part_fe, link = ra$link,
                                                    offset = ra$offset, weights = ra$weights), ra$dots))
    ## the fit's cluster when it has one, else its fixed-effects units (relabelled within each drawn block), as the
    ## fits' own bootstraps and ame(B > 0) resample
    pt <- .fd_refit_boot(pt, refit, spec_of, data, B, level, cores, cl = if (!is.null(ra$cluster)) ra$cluster else ra$fe,
                         relabel = if (!is.null(ra$fe)) unique(c(ra$fe, if (identical(ra$part_fe, ra$fe)) ra$part_fe)))
  }
  pt
}

#' @rdname first_difference.hurdle_cpb
#' @method first_difference zi_cpb
#' @export
first_difference.zi_cpb <- function(object, variable, from, to, B = 0, level = 0.95,
                                    data = NULL, stage = c("both", "participation",
                                    "intensity", "zero", "count"), cores = 1L, at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  spec_of <- function(f) .fd_spec(f, variable, from, to, at = at, stage = stage)
  pt <- .fd_build(spec_of(object), level)
  if (B > 0) {
    if (is.null(data)) stop("Bootstrap (B > 0) requires the original 'data'.")
    ra <- .fd_refit_args(object, c("fe"))
    ra$support <- .cpb_support(object)
    args <- ra[intersect(names(ra), c("fe", "zero_fe", "method", "max.support", "maxit", "tol", "offset", "weights",
                                       "support"))]
    refit <- function(d) do.call(zi_cpb, c(list(object$formula, data = d, zero = object$zero_formula, se = "none"), args))
    pt <- .fd_refit_boot(pt, refit, spec_of, data, B, level, cores, cl = if (!is.null(ra$cluster)) ra$cluster else ra$fe,
                         relabel = if (!is.null(ra$fe)) unique(c(ra$fe, if (identical(ra$zero_fe, ra$fe)) ra$zero_fe)))
  }
  pt
}

## The settings a first-difference bootstrap refits with: those the fit stored
## (`refit_args`, from 0.2.0), else the fields an earlier object carries. An
## offset or weights given as a vector cannot follow resampled rows of `data`.
.fd_refit_args <- function(object, fields) {
  ra <- object$refit_args
  if (is.null(ra)) ra <- c(object[intersect(fields, names(object))], list(dots = list()))
  for (nm in c("offset", "weights")) {
    v <- ra[[nm]]
    if (!is.null(v) && !(is.character(v) && length(v) == 1L))
      stop("first_difference(B > 0) refits the model on resampled rows of 'data': give '", nm,
           "' to the fit as a column name of 'data'.", call. = FALSE)
  }
  if (is.null(ra$dots)) ra$dots <- list()
  ra
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
#' single-parameter boundary test. Given two fitted models it refuses the cases
#' where that reference distribution does not hold: a covariate-dependent
#' structural-zero probability (several parameters on their boundary), a
#' fixed-effects `zi_cpb()` (one count effect per unit, whose number grows with
#' the panel, so the statistic has no fixed-dimension boundary reference), a
#' separated zero equation, models fit to different data, a zero-truncated
#' nest, and CPB models on different supports. It is not a test of a hurdle
#' against a mixture, which are not nested. Compare those with
#' [compare_models()].
#'
#' @param object Either a model formula --- then `data` is required and the
#'   nested pair is fit internally --- or a fitted model (`cpb`, `gec`,
#'   `zi_cpb`, or `zi_gec`).
#' @param object2 The second argument: the data frame when `object` is a formula
#'   (so `zi_test(y ~ x, mydata)` works), or the second fitted model when
#'   `object` is a fit. In the model interface exactly one of the two models must
#'   be a zero-inflated `zi_cpb` (or `zi_gec`) with intercept-only inflation and
#'   the other its untruncated non-inflated nest (`cpb`, or `gec`) on the same
#'   data.
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
    ## the nest's support-guard warning is not the point of this test: the
    ## untruncated CPB of zero-inflated counts leans toward its Poisson limit
    cpb_fit <- .ud_quiet_guard(cpb(object, data = dd, truncated = FALSE, se = "none", ...))
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
    ## the likelihood-ratio test needs the nest of the same model on the same data: the same rows, the
    ## untruncated count model, and (for the CPB) the same support
    if (!is.null(cpb_fit$n) && !is.null(zi_fit$n) && cpb_fit$n != zi_fit$n)
      stop("The two models were fit on different numbers of observations (", cpb_fit$n, " vs ", zi_fit$n,
           "); the likelihood-ratio test needs both fitted to the same data.", call. = FALSE)
    if (isTRUE(cpb_fit$truncated))
      stop("The non-inflated model is zero-truncated; the nest of the zero-inflated model is the untruncated ",
           "fit (truncated = FALSE) on the same data.", call. = FALSE)
    if (inherits(zi_fit, "zi_cpb") && !identical(.cpb_support(cpb_fit), .cpb_support(zi_fit)))
      stop("The two models use different CPB supports (", .cpb_support(cpb_fit), " vs ", .cpb_support(zi_fit),
           "; an object saved before version 0.2.0 reads as the floor support); refit them on one support.",
           call. = FALSE)
    ## the boundary reference distribution holds for one parameter on its boundary at the two models' joint
    ## maxima: not with several inflation parameters (every structural-zero probability at zero is several
    ## parameters at their bounds), not for a fixed-effects fit (its unit effects are a first-stage plug-in, so
    ## its log-likelihood is not the joint maximum, or one per unit, so the parameters grow with the panel), and
    ## not on a divergent ridge (separation)
    if (length(zi_fit$zero_coef) > 1L)
      stop("The zero-inflated model has a covariate-dependent structural-zero probability: under no inflation each ",
           "of its zero-equation parameters is unidentified, and the likelihood-ratio statistic has no ",
           "0.5 chi^2_0 + 0.5 chi^2_1 reference distribution. Test with intercept-only inflation (zero = ~ 1), or ",
           "compare the models with compare_models().", call. = FALSE)
    if (!is.null(zi_fit$fe))
      stop("zi_test() needs both models at their joint maximum, with a fixed number of parameters. A fixed-effects ",
           "zi_cpb() estimates one count effect per unit, whose number grows with the panel (on the floor support, ",
           "method = \"ml\" holds them at a first-stage plug-in), so the likelihood-ratio statistic has no boundary ",
           "reference distribution; compare the models with compare_models().", call. = FALSE)
    if (length(.ud_separated(zi_fit)))
      stop("The zero-inflated model's zero equation is separated (no useful bound on ",
           paste(.ud_separated(zi_fit), collapse = ", "), "); its log-likelihood is a supremum along a divergent ",
           "ridge, not a maximum. See ?zi_cpb, Separation.", call. = FALSE)
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
