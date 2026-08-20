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
#' @param newdata Optional data frame of covariate profiles.
#' @param type `"response"` (marginal E(Y), zeros included), `"participation"`
#'   (P(Y > 0)), or `"intensity"` (E(Y | Y > 0)).
#' @param offset Optional offset (log scale) for the intensity when predicting on
#'   `newdata`; a numeric vector or a column name in `newdata`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict hurdle_count
#' @export
predict.hurdle_count <- function(object, newdata = NULL,
                                 type = c("response", "participation", "intensity"), offset = NULL, ...) {
  type <- match.arg(type); fam <- .count_fam(object$family)
  if (is.null(newdata)) { p <- object$p_full; mu <- object$lambda_full }
  else {
    p  <- as.numeric(stats::predict(object$participation, newdata = newdata, type = "response"))
    mu <- predict(object$intensity, newdata = newdata, type = "response", offset = offset)
  }
  if (type == "participation") return(as.numeric(p))
  ztmean <- fam$meanfun(mu, object$theta) / (1 - fam$p0(mu, object$theta))   # E(Y | Y > 0)
  if (type == "intensity") return(ztmean)
  as.numeric(p) * ztmean                                 # marginal E(Y)
}

#' Predict from a zero-inflated Poisson/NB fit
#'
#' @param object A `"zi_count"` object.
#' @param newdata Optional data frame of covariate profiles.
#' @param type `"response"` (marginal E(Y)), `"count"` (equivalently
#'   `"intensity"`; the count component's mean E(Y) of that component), or
#'   `"zero"` (structural-zero probability). `"intensity"` is accepted as a
#'   synonym for `"count"` for consistency with [predict.zi_cpb()].
#' @param offset Optional offset (log scale) for the count component when
#'   predicting on `newdata`; a numeric vector or a column name in `newdata`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict zi_count
#' @export
predict.zi_count <- function(object, newdata = NULL, type = c("response", "count", "zero", "intensity"),
                             offset = NULL, ...) {
  type <- match.arg(type); if (type == "intensity") type <- "count"
  fam <- .count_fam(object$family)
  if (is.null(newdata)) { pistar <- object$pi_full; lam <- object$lambda_full }
  else {
    linkinv <- stats::make.link(object$link)$linkinv
    Xc <- stats::model.matrix(stats::delete.response(stats::terms(object$formula)), newdata)
    Zz <- stats::model.matrix(object$zero.formula, newdata)
    eta <- as.numeric(Xc[, names(object$coefficients), drop = FALSE] %*% object$coefficients)
    if (!is.null(offset)) {
      ov <- if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset
      eta <- eta + as.numeric(ov)
    }
    lam    <- exp(eta)
    pistar <- as.numeric(linkinv(Zz[, names(object$zero.coefficients), drop = FALSE] %*% object$zero.coefficients))
  }
  if (type == "zero") return(pistar)
  cmean <- fam$meanfun(lam, object$theta)                # E(Y) of the count component
  if (type == "count") return(cmean)
  (1 - pistar) * cmean
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
#'   every `irr()` method: columns `term`, `equation` (`"count"` or
#'   `"binary"`), `ratio` (`"IRR"` or `"OR"`), `estimate`, `lower`, `upper`,
#'   and `method` (the interval source; `"none"` with `NA` bounds when no
#'   covariance is available). Two-part models stack both equations.
#' @name irr.count
NULL

## OR rows for an embedded binary glm (factor-dummy terms dropped from display)
.or_rows <- function(glmfit, level) {
  sm <- summary(glmfit)$coefficients
  keep <- !grepl("^factor\\(", rownames(sm))
  b <- stats::setNames(sm[keep, 1], rownames(sm)[keep])
  .ud_irr_wald(b, sm[keep, 2], level, "binary", "OR")
}

#' @rdname irr.count
#' @method irr hurdle_count
#' @export
irr.hurdle_count <- function(object, level = 0.95, ...)
  rbind(irr(object$intensity, level = level), .or_rows(object$participation, level))

#' @rdname irr.count
#' @method irr zi_count
#' @export
irr.zi_count <- function(object, level = 0.95, ...)
  rbind(.ud_irr_wald(object$coefficients,
                     if (!is.null(object$vcov)) sqrt(diag(object$vcov)) else NULL,
                     level, "count", "IRR"),
        .ud_irr_wald(object$zero.coefficients, NULL, level, "binary", "OR"))

## --- first-difference (discrete change), equation-decomposed for two-part ---

#' First differences for the matched count families
#'
#' The discrete-change effect on the expected count E(Y) as `variable` moves from
#' `from` to `to`, holding the other covariates at their sample means. For
#' `count_reg` this is a single number with a delta-method interval; for the
#' two-part models it is decomposed exactly into extensive (participation /
#' non-structural-zero) and intensive (count) channels that sum to the total.
#'
#' @param object A fitted `count_reg`, `gec`, `gec_fe`, `cpb_fe`,
#'   `hurdle_count`, `zi_count`, `hurdle_gec`, or `zi_gec` model.
#' @param variable Name of the covariate to change (must be in the model).
#' @param from,to The two values of `variable`.
#' @param level Confidence level for the delta-method intervals (all methods on
#'   this page that can compute one).
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
first_difference.count_reg <- function(object, variable, from, to, level = 0.95, ...) {
  .fd_dots(...)
  fam <- .count_fam(object$family); xr <- colMeans(object$X); b <- object$coefficients
  if (!variable %in% names(b)) stop("'variable' is not a covariate in the model.")
  mk <- function(v) { x <- xr; if (variable %in% names(x)) x[variable] <- v; x }
  mfun <- function(bb, v) fam$meanfun(exp(sum(mk(v) * bb)), object$theta)
  ci <- .fd_delta_ci(function(bb) mfun(bb, to) - mfun(bb, from), b, object$vcov, level)
  .ud_fd(component = "mean", from = mfun(b, from), to = mfun(b, to),
         lower = ci[1], upper = ci[2],
         method = if (is.null(object$vcov)) "none" else "delta")
}

## LEVEL rows for a two-part fit: binary-stage probability p, count-stage mean m,
## and marginal p*m, each at `from` and `to`. Delta intervals per row where a
## coherent covariance exists for that row's parameters (the hurdle factorizes,
## so the marginal combines independent blocks).
.twopart_fd <- function(pfun, mfun, bpar, mpar, bvc, mvc,
                        from, to, level, labels) {
  pf <- pfun(bpar, from); pt_ <- pfun(bpar, to)
  mf <- mfun(mpar, from); mt <- mfun(mpar, to)
  ci_p <- .fd_delta_ci(function(bb) pfun(bb, to) - pfun(bb, from), bpar, bvc, level)
  ci_m <- .fd_delta_ci(function(bb) mfun(bb, to) - mfun(bb, from), mpar, mvc, level)
  ci_g <- if (!is.null(bvc) && !is.null(mvc)) {
    par <- c(bpar, mpar); nb <- length(bpar)
    vc  <- matrix(0, length(par), length(par))
    vc[seq_len(nb), seq_len(nb)] <- bvc
    vc[nb + seq_along(mpar), nb + seq_along(mpar)] <- mvc
    .fd_delta_ci(function(pp) pfun(pp[seq_len(nb)], to)  * mfun(pp[-seq_len(nb)], to) -
                              pfun(pp[seq_len(nb)], from) * mfun(pp[-seq_len(nb)], from),
                 par, vc, level)
  } else c(NA_real_, NA_real_)
  .ud_fd(component = labels,
         from  = c(pf, mf, pf * mf),
         to    = c(pt_, mt, pt_ * mt),
         lower = c(ci_p[1], ci_m[1], ci_g[1]),
         upper = c(ci_p[2], ci_m[2], ci_g[2]),
         method = c(if (is.null(bvc)) "none" else "delta",
                    if (is.null(mvc)) "none" else "delta",
                    if (is.null(bvc) || is.null(mvc)) "none" else "delta"))
}

#' @rdname first_difference.count
#' @method first_difference hurdle_count
#' @export
first_difference.hurdle_count <- function(object, variable, from, to, level = 0.95,
                                          stage = c("both", "participation", "intensity",
                                                    "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  fam <- .count_fam(object$family); linkinv <- stats::make.link(object$link)$linkinv
  inx <- variable %in% names(object$int_xref); inz <- variable %in% names(object$part_xref)
  if (!inx && !inz) stop("'variable' is in neither the intensity nor the participation equation.")
  zr0 <- object$part_xref; xr0 <- object$int_xref
  pfun <- function(g, v) { zr <- zr0
    if (inz && stage %in% c("both", "participation", "zero")) zr[variable] <- v
    linkinv(sum(zr * g)) }
  mfun <- function(b, v) { xr <- xr0
    if (inx && stage %in% c("both", "intensity", "count")) xr[variable] <- v
    lam <- exp(sum(xr * b))
    fam$meanfun(lam, object$theta) / (1 - fam$p0(lam, object$theta)) }
  bvc <- tryCatch(stats::vcov(object$participation), error = function(e) NULL)
  mvc <- object$intensity$vcov
  .twopart_fd(pfun, mfun, object$part_beta, object$int_beta, bvc, mvc,
              from, to, level, c("participation", "intensity", "marginal"))
}

#' @rdname first_difference.count
#' @method first_difference zi_count
#' @export
first_difference.zi_count <- function(object, variable, from, to, level = 0.95,
                                      stage = c("both", "participation", "intensity",
                                                "zero", "count"), ...) {
  .fd_dots(...); stage <- match.arg(stage)
  linkinv <- stats::make.link(object$link)$linkinv; fam <- .count_fam(object$family)
  inx <- variable %in% names(object$xref); inz <- variable %in% names(object$zref)
  if (!inx && !inz) stop("'variable' is in neither the count nor the inflation equation.")
  zr0 <- object$zref; xr0 <- object$xref
  ## rows: P(structural zero) level, count-component mean, marginal (1-pi)*m
  pfun <- function(g, v) { zr <- zr0
    if (inz && stage %in% c("both", "participation", "zero")) zr[variable] <- v
    linkinv(sum(zr * g)) }
  mfun <- function(b, v) { xr <- xr0
    if (inx && stage %in% c("both", "intensity", "count")) xr[variable] <- v
    fam$meanfun(exp(sum(xr * b)), object$theta) }
  pi_f <- pfun(object$zero.coefficients, from); pi_t <- pfun(object$zero.coefficients, to)
  m_f  <- mfun(object$coefficients, from);      m_t  <- mfun(object$coefficients, to)
  ci_m <- .fd_delta_ci(function(bb) mfun(bb, to) - mfun(bb, from),
                       object$coefficients, object$vcov, level)
  .ud_fd(component = c("zero", "count", "marginal"),
         from  = c(pi_f, m_f, (1 - pi_f) * m_f),
         to    = c(pi_t, m_t, (1 - pi_t) * m_t),
         lower = c(NA_real_, ci_m[1], NA_real_),
         upper = c(NA_real_, ci_m[2], NA_real_),
         method = c("none", if (is.null(object$vcov)) "none" else "delta", "none"))
}
