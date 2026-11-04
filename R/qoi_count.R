## ---------------------------------------------------------------------------
## Equation-specific quantities of interest for the matched Poisson/NB families:
## predict() (marginal E(Y), or a specific equation) and irr() (rate ratios for
## the count component, odds ratios for the binary component), mirroring the CPB
## hurdle/ZI. For discrete-change effects use first_difference(); marginaleffects
## compatibility is discussed in ?count_reg (marginaleffects gates on a fixed
## class whitelist, so avg_slopes() does not yet recognize these classes).
## ---------------------------------------------------------------------------

## --- equation-specific predict for the two-part Poisson/NB models ----------

#' Predict from a hurdle Poisson/NB fit
#'
#' @param object A `"hurdle_count"` object.
#' @param newdata Optional data frame of covariate profiles. A fixed-effects
#'   unit an equation has not seen, or every row when its column is absent,
#'   enters that equation at the reference unit, and a unit at a limit at that
#'   limit (see `fe,part_fe` in [hurdle_count()]).
#' @param type `"response"` (marginal E(Y), zeros included), `"participation"`
#'   (P(Y > 0)), or `"intensity"` (E(Y | Y > 0)).
#' @param offset Optional offset (log scale) for the intensity when predicting on
#'   `newdata`; a numeric vector or a column name in `newdata`.
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict hurdle_count
#' @export
predict.hurdle_count <- function(object, newdata = NULL,
                                 type = c("response", "participation", "intensity"), offset = NULL,
                                 se.fit = FALSE, interval = c("none", "confidence"),
                                 level = 0.95, ...) {
  type <- match.arg(type); fam <- .count_fam(object$family)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  if (is.null(newdata)) { p <- object$p_full; mu <- object$lambda_full }
  else {
    p  <- .part_at_limit(.ud_part_prob(object, newdata), object, newdata)
    mu <- exp(.count_reg_eta(object$intensity, newdata, offset, xref = object$int_xref))   # natural parameter
  }
  bad <- is.na(mu) | is.na(p); mu[bad] <- 1                    # a missing covariate predicts NA
  if (type == "participation") return(.ud_mask(as.numeric(p), bad))
  ztmean <- .ud_mask(fam$meanfun(mu, object$theta) / pmax(1 - fam$p0(mu, object$theta), 1e-12), bad)   # E(Y | Y > 0)
  if (type == "intensity") return(ztmean)
  as.numeric(p) * ztmean                                 # marginal E(Y)
}

#' Predict from a zero-inflated Poisson/NB fit
#'
#' @param object A `"zi_count"` object.
#' @param newdata Optional data frame of covariate profiles. A fixed-effects
#'   unit an equation has not seen, or every row when its column is absent,
#'   enters that equation at the reference unit. A `zero_fe` unit with no
#'   positive count gets the structural-zero probability one and one whose every
#'   count is positive zero, and an `fe` unit with no positive count a count
#'   mean of zero (see [compare_models()], Degrees of freedom).
#' @param type `"response"` (marginal E(Y)), `"count"` (equivalently
#'   `"intensity"`; the count component's mean E(Y) of that component), or
#'   `"zero"` (structural-zero probability). `"intensity"` is accepted as a
#'   synonym for `"count"` for consistency with [predict.zi_cpb()].
#' @param offset Optional offset (log scale) of the count component for the
#'   `newdata` rows; by default the fit's own offset column is read from
#'   `newdata` (see [predict.cpb()]).
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict zi_count
#' @export
predict.zi_count <- function(object, newdata = NULL, type = c("response", "count", "zero", "intensity"),
                             offset = NULL, se.fit = FALSE, interval = c("none", "confidence"),
                             level = 0.95, ...) {
  type <- match.arg(type); if (type == "intensity") type <- "count"
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  fam <- .count_fam(object$family); at0 <- FALSE
  if (is.null(newdata)) { pistar <- object$pi_full; lam <- object$lambda_full }
  else {
    linkinv <- stats::make.link(object$link)$linkinv
    ez <- .zi_count_etas(object, newdata, offset)
    lam    <- exp(ez$count)
    pistar <- .zero_at_limit(as.numeric(linkinv(ez$zero)), object, newdata)
    at0 <- .count_at_limit(object, newdata)                  # a unit with no positive count: its mean at zero
  }
  bad <- is.na(lam) | is.na(pistar); lam[bad] <- 1             # a missing covariate or offset predicts NA
  if (type == "zero") return(.ud_mask(pistar, bad))
  cmean <- fam$meanfun(lam, object$theta)                # E(Y) of the count component
  cmean[at0] <- 0
  cmean <- .ud_mask(cmean, bad)
  if (type == "count") return(cmean)
  (1 - pistar) * cmean
}

## The count and zero equations' linear predictors of a zi_count() fit in new data (the count one with its offset).
## A unit of an `fe` or `zero_fe` column the fit has not seen (or every row, when `newdata` lacks the column) enters
## at that equation's reference unit (`xref`, `zref`), and with `all_factors` (held-out scoring) an unseen level of
## any factor; the rows the fit can evaluate are computed as before.
.zi_count_etas <- function(object, newdata, offset = NULL, all_factors = FALSE) {
  Tc <- stats::terms(object$formula)
  Tz <- if (!is.null(object$zero_terms)) object$zero_terms else stats::terms(object$zero.formula)
  bc <- object$coefficients; gz <- object$zero.coefficients
  eq <- function(Tm, coefs, contr, xref, fe) {
    ## the fit's levels are shared by its two equations; each equation's model frame takes its own variables'
    ## (another variable's levels make model.frame() warn that it is not a factor)
    vars <- vapply(as.list(attr(stats::delete.response(Tm), "variables"))[-1L], function(v) paste(deparse(v), collapse = ""), "")
    xl <- object$levels[intersect(names(object$levels), vars)]
    if (!length(xl)) xl <- NULL
    if (!.ud_has_unseen(xl, newdata, fe, all_factors)) {
      X <- .ud_newdata_matrix(Tm, newdata, xl, contr, names(coefs))
      return(as.numeric(X[, names(coefs), drop = FALSE] %*% coefs))
    }
    .ud_dummy_eta(stats::delete.response(Tm), xl, coefs, newdata, contr, xref, fe, all_factors)
  }
  ec <- eq(Tc, bc, object$contrasts, object$xref, object$fe)
  list(count = .ud_newdata_offset(object, newdata, offset, nrow(newdata)) + ec,
       zero = eq(Tz, gz, object$zero_contrasts, object$zref, object$zero_fe))
}

## --- equation-specific rate/odds ratios -----------------------------------

#' Rate and odds ratios for the matched count families
#'
#' Rate ratios (`exp` of the count/intensity coefficients, column `IRR`) and, for
#' the two-part models, odds ratios of the binary stage (column `OR`), by
#' equation. Element names match the CPB family (`intensity` for the count
#' component; `participation` for the hurdle stage, `zero` for the inflation
#' stage).
#'
#' @param object A fitted model from this package.
#' @param level Confidence level.
#' @param ... Unused.
#' @return A `"ud_irr"` data frame -- the package-wide ratio contract shared by
#'   every `irr()` method: columns `term`, `equation` (`"count"`,
#'   `"participation"` for a hurdle's binary stage, or `"inflation"` for a
#'   zero-inflated model's structural-zero stage; the two point in opposite
#'   directions, an odds ratio above one raising P(Y > 0) in the first and
#'   P(structural zero) in the second), `ratio` (`"IRR"`, `"OR"`, or the
#'   parameter's ratio named below), `estimate`, `lower`, `upper`, and `method`
#'   (the interval source; `"none"` with `NA` bounds when no covariance is
#'   available). Two-part models stack both equations. The count rows are
#'   \eqn{\exp(\beta)}, the ratio of the family's natural parameter. That is a
#'   ratio of means for the Poisson, the negative binomial, the
#'   mean-parameterized COM-Poisson (`"mpcmp"`), the hyper-Poisson, and the
#'   generalized Poisson with \eqn{\lambda \ge 0} (`"IRR"`). It is not for the
#'   others, whose rows name what it is, with a message: the rate-parameterized
#'   COM-Poisson (`"compois"`, `"lambda ratio"`; its mean moves roughly as the
#'   rate to the power \eqn{1/\nu}, so that \eqn{\exp(\beta)} is roughly the
#'   mean ratio to the power \eqn{\nu}), the EWP (`"centrality ratio"`), the
#'   gamma count (`"long-run rate ratio"`), the double Poisson
#'   (`"location ratio"`) and the generalized Poisson below zero
#'   (`"mu ratio"`). The CPB and the underdispersed GEC keep `"IRR"`, a ratio of
#'   their rate, which their mean follows only approximately (the rate
#'   renormalized over a finite support). For all of these, [first_difference()]
#'   and [ame()] report the exact mean contrast.
#' @name irr.count
NULL

## OR rows of a hurdle's participation equation: its covariance (participation.R), its part_fe dummies left out,
## its separated terms NA
.or_rows <- function(object, level) {
  b <- object$part_beta; V <- .part_vcov(object)
  keep <- !.is_fe_col(names(b), object$part_fe) & !is.na(b)
  se <- if (is.null(V)) NULL else sqrt(pmax(diag(V), 0))[names(b)][keep]
  .irr_mask(.binary_ratio_rows(b[keep], se, level, "participation", object$link),
            object$part_separation$diverging, "participation")
}

## exp(coefficient) of a binary equation is an odds ratio only under the logit link; under probit or cloglog
## no ratio is reported for that equation (first_difference() gives its effects on the probability scale)
.binary_ratio_rows <- function(b, se, level, equation, link, method = "wald") {
  if (is.null(link) || identical(link, "logit")) return(.ud_irr_wald(b, se, level, equation, "OR", method = method))
  message("irr(): the ", equation, " equation uses the ", link, " link, under which exp(coefficient) is not an ",
          "odds ratio, so no ratio is reported for it; first_difference() gives its effects on the probability scale.")
  NULL
}

#' @rdname irr.count
#' @method irr hurdle_count
#' @export
irr.hurdle_count <- function(object, level = 0.95, ...)
  rbind(irr(object$intensity, level = level), .or_rows(object, level))

#' @rdname irr.count
#' @method irr zi_count
#' @export
irr.zi_count <- function(object, level = 0.95, ...) {
  g <- object$zero.coefficients; keep <- !.is_fe_col(names(g), object$zero_fe)
  se <- if (!is.null(object$se.zero)) unname(object$se.zero)[keep]        # named count:/zero: by the covariance
  rbind(.irr_drop_fe(.count_irr_rows(object, object$coefficients,
                                     if (!is.null(object$vcov)) sqrt(pmax(diag(object$vcov), 0)), level), object$fe),
        .irr_mask(.binary_ratio_rows(g[keep], se, level, "inflation", object$link),
                  object$separation$diverging, "inflation"))
}


## --- first-difference (discrete change), equation-decomposed for two-part ---

#' First differences for the matched count families
#'
#' The discrete-change effect on the expected count E(Y) as `variable` moves from
#' `from` to `to`, holding the other covariates at their sample means (or at
#' the values `at` gives them). The interval follows the package's one rule (see
#' [first_difference()]): percentile intervals from the fit's stored bootstrap
#' draws for the GEC, fixed-effects CPB and zero-inflated GEC fits, the delta
#' method over the full covariance, the dispersion parameter included, for the
#' matched count families and the hurdles' margins. For
#' `count_reg` this is a single number; for the
#' two-part models it is reported on each margin: the binary stage's
#' probability (participation, or the structural-zero probability), the count
#' stage's mean, and the marginal mean, their product, each at `from` and `to`.
#' The first two are different quantities and their changes do not add up to
#' the marginal change, which is
#' \eqn{\Delta(pm) = \Delta p\, m_0 + p_0 \Delta m + \Delta p\, \Delta m}.
#'
#' @param object A fitted `count_reg`, `gec`, `gec_fe`, `cpb_fe`,
#'   `hurdle_count`, `zi_count`, `hurdle_gec`, or `zi_gec` model.
#' @param variable Name of the covariate to change (must be in the model).
#' @param from,to The two values of `variable`.
#' @param level Confidence level of the intervals.
#' @param at Optional named numeric vector fixing other model-matrix columns at
#'   given values instead of their means (see [first_difference()]); in a
#'   two-part model a column of both equations is fixed in both.
#' @param stage For the two-part models, which equation(s) the change is applied
#'   to. `"both"` (default) moves the variable wherever it appears;
#'   `"intensity"`/`"count"` or `"participation"`/`"zero"` moves it in only
#'   that equation, holding it at the reference in the other. For a variable in
#'   only one equation all options coincide; for a variable in both, `"both"`
#'   gives the total effect and the single-stage options the partial effect
#'   through that margin. All component rows are always returned.
#' @param ... Unused; unknown arguments error.
#' @return A `"ud_fd"` data frame -- the package-wide first-difference contract
#'   shared by every `first_difference()` method: columns `component`, `from`,
#'   `to`, `diff`, `lower`, `upper`, `method`. Single-equation fits return one
#'   row (`component = "mean"`); two-part fits return one row per margin level
#'   (binary-stage probability, count-stage mean, marginal mean). `method`
#'   records the uncertainty source per row (`"delta"`, a bootstrap label, or
#'   `"none"` with `NA` bounds).
#' @name first_difference.count
NULL

#' @rdname first_difference.count
#' @method first_difference count_reg
#' @export
first_difference.count_reg <- function(object, variable, from, to, level = 0.95, at = NULL, ...) {
  .fd_dots(...)
  .fd_build(.fd_spec(object, variable, from, to, at = at), level)
}

#' @rdname first_difference.count
#' @method first_difference hurdle_count
#' @export
first_difference.hurdle_count <- function(object, variable, from, to, level = 0.95,
                                          stage = c("both", "participation", "intensity",
                                                    "zero", "count"), at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  .fd_build(.fd_spec(object, variable, from, to, at = at, stage = stage), level)
}

#' @rdname first_difference.count
#' @method first_difference zi_count
#' @export
first_difference.zi_count <- function(object, variable, from, to, level = 0.95,
                                      stage = c("both", "participation", "intensity",
                                                "zero", "count"), at = NULL, ...) {
  .fd_dots(...); stage <- match.arg(stage)
  .fd_build(.fd_spec(object, variable, from, to, at = at, stage = stage), level)
}
