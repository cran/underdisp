## ---------------------------------------------------------------------------
## The parameters of each class's quantities of interest, and the intervals of
## predictions and average marginal effects, by the rule of .qoi_ci()
## (qoi_unified.R). A quantity is predict() on the fit with its parameters
## replaced by a draw or a perturbation, so its value at the estimate is the
## fit's own prediction.
## ---------------------------------------------------------------------------

## A fit's parameters `theta` (named as the stored draws' columns or the covariance's dimnames), the source of their
## uncertainty (`draws` or `vc`, with `recenter`, `hold` and `bounds` as in .qoi_ci()), `set(th)`, the fit with its parameters
## replaced by `th`, and `own(newdata)`, the rows that enter at a parameter neither the draws nor the covariance
## carry: a concentrated fixed-effects unit's own effect (the bootstrap resamples units, so a unit's effect is not a
## column of the draws) or a zero_fe unit dummy (not the same units across replicates).
.qoi_par <- function(object) UseMethod(".qoi_par")
## the uncertainty fields of .qoi_par(), as a first difference's specification takes them
.qoi_unc <- function(p) p[intersect(c("theta", "draws", "vc", "recenter", "hold", "bounds"), names(p))]
.qoi_none <- function(newdata) rep(FALSE, nrow(newdata))
## the rows of `newdata` whose column `v` names one of `units`
.qoi_rows_in <- function(newdata, v, units)
  if (is.null(v) || !v %in% names(newdata)) .qoi_none(newdata) else as.character(newdata[[v]]) %in% as.character(units)
## the levels of the zero_fe dummies `zfe` of a zero-inflated fit
.qoi_zero_levels <- function(object, zfe) if (is.null(zfe)) character(0) else object$zero_xlev[[paste0("factor(", zfe, ")")]]

.qoi_par.cpb <- function(object) {
  tt <- .fd_theta1(object, object$alpha); b <- names(object$coefficients)
  c(tt, list(set = function(th) { object$coefficients[b] <- th[b]; object$alpha <- th[[".disp"]]; object },
             own = .qoi_none))
}

.qoi_par.cpb_fe <- function(object) {
  tt <- .fd_theta1(object, object$alpha, .fe_ref(object)); b <- names(object$coefficients)
  c(tt, list(recenter = TRUE, bounds = list(.disp = c(0, 1)),
             set = function(th) { object$coefficients[b] <- th[b]; object$alpha <- th[[".disp"]]
                                  object$fe_ref <- th[[".fe_ref"]]; object },
             own = function(nd) .qoi_rows_in(nd, object$fe_var, names(object$fe))))
}

## gec() and gec_fe() (whose class is c("gec_fe", "gec"))
.qoi_par.gec <- function(object) {
  fe <- inherits(object, "gec_fe"); b <- names(object$coefficients)
  tt <- .fd_theta1(object, object$delta, if (fe) .fe_ref(object))
  c(tt, list(recenter = fe, bounds = list(.disp = c(0, Inf)),
             set = function(th) { object$coefficients[b] <- th[b]; object$delta <- th[[".disp"]]
                                  if (fe) object$fe_ref <- th[[".fe_ref"]]
                                  object },
             own = if (fe) function(nd) .qoi_rows_in(nd, object$fe_var, names(object$fe)) else .qoi_none))
}

.qoi_par.count_reg <- function(object) {
  fam <- .count_fam(object$family); b <- names(object$coefficients)
  sh <- if (fam$nshape) fam$shape_name else character(0)
  list(theta = c(object$coefficients, if (fam$nshape) stats::setNames(fam$shape_link(object$theta), sh)),
       vc = object$vcov_full,
       set = function(th) { object$coefficients[b] <- th[b]
                            if (fam$nshape) object$theta <- unname(fam$shape_inv(th[sh]))
                            object },
       own = .qoi_none)
}

## a two-part fit with its participation glm's coefficients replaced by the "participation:" entries of `th` (an
## aliased coefficient stays NA)
.qoi_set_part <- function(object, th) {
  g <- object$part_beta; nm <- names(g)[!is.na(g)]
  v <- th[paste0("participation:", nm)]
  object$participation$coefficients[nm] <- v; object$part_beta[nm] <- v
  object
}

.qoi_par.hurdle_count <- function(object) {
  fam <- .count_fam(object$family); pp <- .fd_part(object); b <- object$int_beta
  sh <- if (fam$nshape) paste0("intensity:", fam$shape_name) else character(0)
  theta <- c(pp$theta, stats::setNames(b, paste0("intensity:", names(b))),
             if (fam$nshape) stats::setNames(fam$shape_link(object$theta), sh))
  iv <- object$intensity$vcov_full
  if (!is.null(iv)) { nm <- paste0("intensity:", rownames(iv)); dimnames(iv) <- list(nm, nm) }
  list(theta = theta, vc = .qoi_vc(names(theta), pp$vc, iv),
       set = function(th) {
         object <- .qoi_set_part(object, th); v <- th[paste0("intensity:", names(b))]
         object$intensity$coefficients[names(b)] <- v; object$int_beta[names(b)] <- v
         if (fam$nshape) object$theta <- object$intensity$theta <- unname(fam$shape_inv(th[sh]))
         object },
       own = .qoi_none)
}

.qoi_par.zi_count <- function(object) {
  fam <- .count_fam(object$family); b <- object$coefficients; g <- object$zero.coefficients
  sh <- if (fam$nshape) fam$shape_name else character(0)
  theta <- c(stats::setNames(b, paste0("count:", names(b))), stats::setNames(g, paste0("zero:", names(g))),
             if (fam$nshape) stats::setNames(fam$shape_link(object$theta), sh))
  list(theta = theta, vc = object$vcov_full,
       set = function(th) {
         object$coefficients[names(b)] <- th[paste0("count:", names(b))]
         object$zero.coefficients[names(g)] <- th[paste0("zero:", names(g))]
         if (fam$nshape) object$theta <- unname(fam$shape_inv(th[sh]))
         object },
       own = .qoi_none)
}

.qoi_par.hurdle_cpb <- function(object) {
  int <- object$intensity; fe <- inherits(int, "cpb_fe"); pp <- .fd_part(object); b <- object$int_beta
  theta <- c(pp$theta, stats::setNames(b, paste0("intensity:", names(b))), alpha = int$alpha,
             if (fe) c(fe_ref = .fe_ref(int)))
  iv <- .fd_draws_vc(int, c(names(b), "alpha", if (fe) "fe_ref"), names(b), "intensity:")
  list(theta = theta, vc = .qoi_vc(names(theta), pp$vc, iv),
       set = function(th) {
         object <- .qoi_set_part(object, th); v <- th[paste0("intensity:", names(b))]
         object$intensity$coefficients[names(b)] <- v; object$int_beta[names(b)] <- v
         object$intensity$alpha <- th[["alpha"]]
         if (fe) object$intensity$fe_ref <- th[["fe_ref"]]
         object },
       own = if (fe) function(nd) .qoi_rows_in(nd, object$fe, names(int$fe)) else .qoi_none)
}

.qoi_par.hurdle_gec <- function(object) {
  pp <- .fd_part(object); b <- object$int_beta
  theta <- c(pp$theta, stats::setNames(b, paste0("intensity:", names(b))), delta = object$delta)
  iv <- .fd_draws_vc(object$intensity, c(names(b), "delta"), names(b), "intensity:")
  list(theta = theta, vc = .qoi_vc(names(theta), pp$vc, iv),
       set = function(th) {
         object <- .qoi_set_part(object, th); v <- th[paste0("intensity:", names(b))]
         object$int_beta[names(b)] <- v; object$intensity$coefficients[names(b)] <- v
         object$delta <- object$intensity$delta <- th[["delta"]]
         object },
       own = .qoi_none)
}

.qoi_par.zi_cpb <- function(object) {
  fe <- !is.null(object$fe); b <- object$int_beta; g <- object$zero_beta
  theta <- c(stats::setNames(b, paste0("count:", names(b))), stats::setNames(g, paste0("zero:", names(g))),
             alpha = object$alpha, if (fe) c(fe_ref = object$int_fe_ref))
  draws <- .ud_draws(object); zfe <- object$refit_args$zero_fe
  list(theta = theta, draws = draws, recenter = fe, bounds = list(alpha = c(0, 1)), hold = .fd_zero_hold(names(theta), draws),
       set = function(th) {
         v <- th[paste0("count:", names(b))]; w <- th[paste0("zero:", names(g))]
         object$coefficients[names(b)] <- v; object$int_beta[names(b)] <- v
         object$zero_coef[names(g)] <- w; object$zero_beta[names(g)] <- w; object$zero.coefficients[names(g)] <- w
         object$alpha <- th[["alpha"]]
         if (fe) object$int_fe_ref <- th[["fe_ref"]]
         object },
       own = function(nd) .qoi_rows_in(nd, object$fe, names(object$fe_hat)) |
                          .qoi_rows_in(nd, zfe, .qoi_zero_levels(object, zfe)))
}

.qoi_par.zi_gec <- function(object) {
  b <- object$coefficients; g <- object$zero_coef
  theta <- c(stats::setNames(b, paste0("count:", names(b))), stats::setNames(g, paste0("zero:", names(g))),
             delta = object$delta)
  draws <- .ud_draws(object)
  list(theta = theta, draws = draws, hold = .fd_zero_hold(names(theta), draws),
       set = function(th) {
         object$coefficients[names(b)] <- th[paste0("count:", names(b))]
         w <- th[paste0("zero:", names(g))]; object$zero_coef[names(g)] <- w; object$zero.coefficients[names(g)] <- w
         object$delta <- th[["delta"]]
         object },
       own = function(nd) .qoi_rows_in(nd, object$zero_fe, .qoi_zero_levels(object, object$zero_fe)))
}

.qoi_check_level <- function(level)
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("'level' must be a single number between 0 and 1.", call. = FALSE)

## the arguments a predict() call passes on, reduced to the rows `ok` of newdata where they have one value per row
.qoi_dots_rows <- function(dots, ok) {
  n <- length(ok)
  lapply(dots, function(v) if (n > 1L && !is.character(v) && length(v) == n) v[ok] else v)
}

.qoi_own_message <- function(fun, k)
  message(fun, ": ", k, if (k == 1L) " row enters" else " rows enter", " at a unit's own fixed effect or its own ",
          "zero-equation dummy, which the fit's bootstrap draws do not carry, so ", if (k == 1L) "it has" else "they have",
          " no interval. A row of a unit the fit has not seen, or newdata without the unit column, enters at the ",
          "reference unit.", if (fun == "ame()") " ame(B = ) gives the average an interval by refitting the model on resampled units.")

## predict() with standard errors or a confidence interval (see predict.cpb(), `se.fit`): the quantity predict()
## gives at each row of `newdata`, as a function of the fit's parameters, through .qoi_ci() row by row
.ud_predict_ci <- function(object, newdata, type, se.fit, interval, level, ...) {
  if (is.null(newdata))
    stop("Standard errors and intervals need 'newdata' (the fit does not keep its data): pass the rows to predict, ",
         "for example the data the model was fit to.", call. = FALSE)
  .qoi_check_level(level)
  dots <- list(...); dots <- dots[!vapply(dots, is.null, logical(1))]
  fit <- as.numeric(do.call(stats::predict, c(list(object, newdata = newdata, type = type), dots)))
  n <- length(fit); ok <- is.finite(fit)
  lo <- hi <- se <- rep(NA_real_, n); m <- rep("none", n)
  if (any(ok)) {
    p <- .qoi_par(object); nd <- newdata[ok, , drop = FALSE]; dk <- .qoi_dots_rows(dots, ok)
    q <- function(th) as.numeric(do.call(stats::predict, c(list(p$set(th), newdata = nd, type = type), dk)))
    ci <- .qoi_ci(q, p$theta, p$draws, p$vc, level, isTRUE(p$recenter), if (is.null(p$hold)) character(0) else p$hold,
                  bounds = p$bounds, each = TRUE)
    own <- p$own(nd)
    if (any(own)) {
      ci$lower[own] <- ci$upper[own] <- ci$se[own] <- NA_real_; ci$method[own] <- "none"
      .qoi_own_message("predict()", sum(own))
    }
    lo[ok] <- ci$lower; hi[ok] <- ci$upper; se[ok] <- ci$se; m[ok] <- ci$method
  }
  out <- if (interval == "confidence") cbind(fit = fit, lwr = lo, upr = hi) else fit
  if (isTRUE(se.fit)) return(list(fit = out, se.fit = se, method = m))
  attr(out, "method") <- m
  out
}

## the variables of a fit's model formulas (the responses and the fixed-effects unit columns left out)
.ud_model_vars <- function(object) {
  fs <- list(object$formula, object$terms, object$int_terms, object$zero_terms, object$zero_formula,
             object$zero.formula, object$part_formula,
             if (inherits(object$participation, "glm")) stats::formula(object$participation))
  fs <- Filter(function(f) inherits(f, "formula"), fs)
  v <- unique(unlist(lapply(fs, function(f) all.vars(stats::delete.response(stats::terms(f))))))
  setdiff(v, .ud_unit_vars(object))
}
## the columns a refit of the fit needs: its formulas' variables, the responses included, and its unit columns
.ud_data_vars <- function(object) {
  fs <- list(object$formula, object$terms, object$int_terms, object$zero_terms, object$zero_formula,
             object$zero.formula, object$part_formula,
             if (inherits(object$participation, "glm")) stats::formula(object$participation))
  fs <- Filter(function(f) inherits(f, "formula"), fs)
  setdiff(unique(c(unlist(lapply(fs, all.vars)), .ud_unit_vars(object))), ".d")
}
.ud_unit_vars <- function(object)
  unique(unlist(Filter(is.character, list(object$fe, object$fe_var, object$zero_fe, object$part_fe,
                                          object$refit_args$zero_fe, object$refit_args$part_fe))))

## `newdata` with column `v` set to `value` in every row, in the column's own type
.ud_set_value <- function(newdata, v, value) {
  x <- newdata[[v]]
  if (length(value) != 1L || is.na(value)) stop("'from' and 'to' must be single values of '", v, "'.", call. = FALSE)
  newdata[[v]] <- if (is.factor(x)) {
    if (!as.character(value) %in% levels(x))
      stop("'", value, "' is not a level of '", v, "' (", paste(levels(x), collapse = ", "), ").", call. = FALSE)
    f <- factor(rep(as.character(value), nrow(newdata)), levels = levels(x)); attr(f, "contrasts") <- attr(x, "contrasts"); f
  } else if (is.character(x)) rep(as.character(value), nrow(newdata))
  else if (is.logical(x)) rep(as.logical(value), nrow(newdata))
  else {
    if (!is.numeric(value)) stop("'from' and 'to' must be numbers for the numeric variable '", v, "'.", call. = FALSE)
    rep(as.numeric(value), nrow(newdata))
  }
  newdata
}

#' Average marginal effects
#'
#' The average, over the rows of `newdata`, of the change in a predicted
#' quantity when one variable changes: its derivative (a numeric variable, with
#' `from` and `to` left out) or the contrast between two of its values (`from`
#' and `to`, for a variable of any type). Every row keeps its own values of the
#' other variables, and the variable enters through the model's formula, so
#' interactions and transformations of it (`x:z`, `I(x^2)`, `log(x)`) are
#' handled. [first_difference()] instead moves one model-matrix column at a
#' single reference profile.
#'
#' The interval follows the rule of [first_difference()]. With stored bootstrap
#' draws, the average is computed at every draw, jointly over the coefficients,
#' the dispersion parameter and a fixed-effects fit's reference unit (the draws
#' of a fixed-effects fit recentered on its estimate), and the percentile
#' interval is taken. Otherwise it uses the delta method over the fit's
#' covariance, and otherwise there is none. A row at a fixed-effects unit's own
#' effect, or at its own zero-equation dummy, depends on a parameter the draws
#' do not carry, so an average over such rows, the sample's average partial
#' effect, has no interval of this kind (`B` below gives it one). Without the unit column in `newdata` every row
#' enters at the reference unit, the mean of the identified units' effects on
#' the log scale (a typical, geometric-mean unit), whose draws the fit carries:
#' that is the average effect for such a unit over the rows' covariates, not
#' the sample's average partial effect, and the two can differ substantially
#' when the unit effects are dispersed.
#'
#' A fixed-effects count family ([count_reg()] with `fe`) carries its unit
#' dummies in its covariance, so its average partial effect has a delta-method
#' interval. The slope of a fit on the floor support (`support = "floor"`) is
#' refused, since its mean jumps at integer ceilings; a contrast is available.
#'
#' `B > 0` gives the average a refit-bootstrap interval instead, for any class:
#' the model is refit on `B` resampled copies of `newdata`, which must then be
#' the data the model was fit to, resampling whole clusters of the fit's
#' `cluster`, else whole units of its first fixed-effects column, else rows (a
#' unit drawn twice becomes two units). Each replicate averages over its
#' own rows at its own unit effects, and the percentile interval is taken.
#' Resampling the units treats them as a sample, so this is an interval for a
#' fixed-effects fit's average partial effect in the population of units they
#' represent, which carries the variability of which units are drawn: with few
#' units, or a few with large effects, it is wide and asymmetric. Each refit is the fit's own call with its
#' data replaced and its standard errors off (on King's support, every
#' replicate of a jackknife-corrected fit is corrected without re-running the
#' validity gate, as in [cpb_fe()]'s own bootstrap, and none is when the gate
#' refused the fit's correction); the call's other arguments
#' are evaluated where `ame()` is called, as [update()] does, and an offset or
#' weights given to the fit as one value per row follow the resampled rows.
#' The run takes about `B` times one fit (see `cores`).
#'
#' @param object A fitted model of any class in the package.
#' @param variable The name of a variable of the model, a column of `newdata`.
#' @param newdata The rows to average over, usually the data the model was fit
#'   to.
#' @param from,to Two values of `variable` for a contrast. Leave both out for the
#'   derivative of a numeric variable.
#' @param type The quantity, as the `type` of the class's `predict()` method.
#'   The default `"response"` is the expected count (for a zero-truncated fit,
#'   such as [cpb()]'s default, the conditional mean E(Y | Y >= 1)); the
#'   two-part models also give their equation-specific quantities.
#' @param eps The step of the central difference for the derivative; by default
#'   1e-4 times the range of `variable` in `newdata`.
#' @param level Confidence level.
#' @param offset Passed to `predict()`.
#' @param weights Optional weights for the average over the rows: a numeric
#'   vector or the name of a column of `newdata`, such as the frequency weights
#'   of a fit to tabulated data (the average then equals the one over the
#'   expanded rows). By default the rows count equally.
#' @param B Refit-bootstrap replicates for a percentile interval; `0`
#'   (default) uses the fit's own inference (see Details).
#' @param cores Worker processes for the refits; see [cpb()].
#' @return A data frame of class `"ud_ame"` with `term`, `contrast` (`"dY/dx"`,
#'   or `"to - from"`), `type`, `estimate`, `lower`, `upper`, `method` (the
#'   interval's source) and `n`, the rows averaged (a row with a missing value
#'   is left out).
#' @seealso [first_difference()], [predict.cpb()] (`se.fit`, `interval`).
#' @examples
#' \donttest{
#' set.seed(1); n <- 400; x <- rnorm(n); g <- factor(sample(c("a", "b"), n, TRUE))
#' d <- data.frame(y = rcpb(n, exp(1 + 0.4 * x + 0.3 * (g == "b")), 0.5), x = x, g = g)
#' fit <- cpb(y ~ x * g, data = d, truncated = FALSE, se = "bootstrap", B = 100)
#' ame(fit, "x", d)
#' ame(fit, "g", d, from = "a", to = "b")
#'
#' ## a fixed-effects fit's average partial effect over its own units
#' pd <- do.call(rbind, lapply(1:30, function(i) {
#'   x <- rnorm(8); data.frame(unit = i, x = x, y = rcpb(8, exp(rnorm(1, 0, 0.5) + 0.4 * x), 0.5))
#' }))
#' ffit <- cpb_fe(y ~ x, data = pd, fe = "unit")
#' ame(ffit, "x", pd, B = 100)
#' }
#' @export
ame <- function(object, variable, newdata, from = NULL, to = NULL, type = "response", eps = NULL, level = 0.95,
                offset = NULL, weights = NULL, B = 0, cores = 1L) {
  .qoi_check_level(level)
  if (!is.numeric(B) || length(B) != 1L || !is.finite(B) || B < 0 || B != floor(B))
    stop("'B' must be a whole number, 0 or more.", call. = FALSE)
  if (B == 1) stop("'B' must be 0 or at least 2.", call. = FALSE)
  if (!is.data.frame(newdata) || !nrow(newdata)) stop("'newdata' must be a data frame with at least one row.",
                                                      call. = FALSE)
  wt <- if (is.null(weights)) rep(1, nrow(newdata))
        else if (is.character(weights) && length(weights) == 1L) {
          if (!weights %in% names(newdata)) stop("'weights' names no column of 'newdata'.", call. = FALSE)
          as.numeric(newdata[[weights]])
        } else as.numeric(weights)
  if (length(wt) != nrow(newdata) || anyNA(wt) || any(wt < 0) || !any(wt > 0))
    stop("'weights' must give one non-negative weight per row of 'newdata', not all zero.", call. = FALSE)
  if (is.null(weights) && !is.null(object$weights) && length(object$weights) == nrow(newdata) &&
      any(object$weights != 1))
    message("ame(): the fit has frequency weights; the rows of 'newdata' count equally unless 'weights' is given.")
  if (!is.character(variable) || length(variable) != 1L || !variable %in% names(newdata))
    stop("'variable' must name one column of 'newdata'.", call. = FALSE)
  if (variable %in% .ud_unit_vars(object))
    stop("'", variable, "' is a fixed-effects unit column, not a covariate.", call. = FALSE)
  vars <- .ud_model_vars(object)
  if (!variable %in% vars)
    stop("'", variable, "' is not a variable of the model (", paste(vars, collapse = ", "), ").", call. = FALSE)
  if (xor(is.null(from), is.null(to)))
    stop("Give both 'from' and 'to' for a contrast, or neither for the derivative.", call. = FALSE)
  slope <- is.null(from)
  if (slope && inherits(object, c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb")) && .cpb_fl(.cpb_support(object)))
    stop("The floor-support CPB's mean jumps at integer ceilings, so its derivative depends on the step; give ",
         "'from' and 'to' for a contrast, or fit King's CPB (support = \"ceiling\").", call. = FALSE)
  if (slope) {
    x <- newdata[[variable]]
    if (!is.numeric(x))
      stop("The derivative needs a numeric variable; give 'from' and 'to' for a contrast of '", variable, "'.",
           call. = FALSE)
    h <- if (is.null(eps)) 1e-4 * diff(range(x, na.rm = TRUE)) else eps
    if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) h <- 1e-4
    nlo <- newdata; nlo[[variable]] <- x - h / 2; nhi <- newdata; nhi[[variable]] <- x + h / 2; div <- h
  } else {
    nlo <- .ud_set_value(newdata, variable, from); nhi <- .ud_set_value(newdata, variable, to); div <- 1
  }
  dots <- if (is.null(offset)) list() else list(offset = offset)
  pr <- function(o, nd, dk) as.numeric(do.call(stats::predict, c(list(o, newdata = nd, type = type), dk)))
  d0 <- (pr(object, nhi, dots) - pr(object, nlo, dots)) / div
  ok <- is.finite(d0) & wt > 0
  if (!any(ok)) stop("No row of 'newdata' gives a finite prediction.", call. = FALSE)
  nlo <- nlo[ok, , drop = FALSE]; nhi <- nhi[ok, , drop = FALSE]; dk <- .qoi_dots_rows(dots, ok)
  wk <- wt[ok] / sum(wt[ok])
  p <- .qoi_par(object)
  own <- sum(p$own(newdata[ok, , drop = FALSE]))
  ci <- if (B > 0) {
    shift <- function(nd) {
      if (!slope) return(list(.ud_set_value(nd, variable, from), .ud_set_value(nd, variable, to)))
      lo <- nd; lo[[variable]] <- nd[[variable]] - h / 2; hi <- nd; hi[[variable]] <- nd[[variable]] + h / 2
      list(lo, hi)
    }
    est <- function(o, nd, w, dk) {
      s <- shift(nd); v <- (pr(o, s[[2L]], dk) - pr(o, s[[1L]], dk)) / div
      k <- is.finite(v) & w > 0
      if (any(k)) sum(w[k] * v[k]) / sum(w[k]) else NA_real_
    }
    .ame_refit_boot(object, newdata, est, wt, offset, B, level, cores, parent.frame())
  } else if (own) {
    .qoi_own_message("ame()", own); list(lower = NA_real_, upper = NA_real_, method = "none")
  } else {
    q <- function(th) { o <- p$set(th); sum(wk * (pr(o, nhi, dk) - pr(o, nlo, dk))) / div }
    .qoi_ci(q, p$theta, p$draws, p$vc, level, isTRUE(p$recenter), if (is.null(p$hold)) character(0) else p$hold,
            bounds = p$bounds)
  }
  out <- data.frame(term = variable, contrast = if (slope) "dY/dx" else paste(to, "-", from), type = type,
                    estimate = sum(wk * d0[ok]), lower = ci$lower, upper = ci$upper, method = ci$method, n = sum(ok),
                    row.names = NULL, stringsAsFactors = FALSE)
  class(out) <- c("ud_ame", "data.frame")
  out
}

## ame(B > 0): the model refit on B resampled copies of `newdata` (the data it was fit to): whole units of its first
## fixed-effects column (a unit drawn twice relabeled as two units), else whole clusters of its call's `cluster`, else
## rows. Each refit is the fit's own call with `data` replaced and its standard errors off, its other arguments
## evaluated in `env` (as update() does), run as a bootstrap replicate (.ud_replicate(): a jackknife fit is corrected
## without its gate); an offset or weights given as one value per row follow the rows. `est(fit, rows, weights,
## predict_args)` is the replicate's average; the interval is the percentile one.
.ame_refit_boot <- function(object, newdata, est, wt, offset, B, level, cores, env) {
  cl <- object$call
  if (is.null(cl)) stop("ame(B > 0) refits the model from its call, which this fit does not store.", call. = FALSE)
  miss <- setdiff(.ud_data_vars(object), names(newdata))
  if (length(miss))
    stop("ame(B > 0) refits the model on resampled rows of 'newdata', which must be the data the model was fit to; ",
         "it lacks ", paste0("'", miss, "'", collapse = ", "), ".", call. = FALSE)
  fun <- get(intersect(class(object), c("cpb_fe", "gec_fe", "cpb", "gec", "count_reg", "hurdle_cpb", "hurdle_gec",
                                        "hurdle_count", "zi_cpb", "zi_gec", "zi_count"))[1L],
             envir = asNamespace("underdisp"), mode = "function")
  args <- as.list(cl)[-1L]
  ccl <- args$cluster
  args <- args[setdiff(names(args), c("data", "se", "B", "cores", "cluster"))]
  ev <- function(a) tryCatch(eval(a, env), error = function(e)
    stop("ame(B > 0) refits the model from its call, whose arguments must still be available where ame() is called: ",
         conditionMessage(e), call. = FALSE))
  args <- lapply(args, ev); ccl <- if (!is.null(ccl)) ev(ccl)
  ## the fit's own CPB support (an object saved before 0.2.0 is a floor-support fit whose call has no `support`)
  if (inherits(object, c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb"))) args$support <- .cpb_support(object)
  ## a replicate skips the jackknife's gate, so it follows the fit's own outcome (no correction when the gate refused)
  bc <- if (!is.null(object$bias_correct)) object$bias_correct else object$intensity$bias_correct
  if (!is.null(args$bias_correct) && !is.null(bc)) args$bias_correct <- bc
  cols <- unlist(Filter(function(v) is.character(v) && length(v) == 1L,
                        c(args[intersect(c("time", "offset", "weights"), names(args))], list(ccl))))
  miss <- setdiff(cols, names(newdata))
  if (length(miss))
    stop("ame(B > 0) refits the model on resampled rows of 'newdata', which must be the data the model was fit to; ",
         "it lacks ", paste0("'", miss, "'", collapse = ", "), ".", call. = FALSE)
  n <- nrow(newdata)
  rowvec <- character(0)
  for (v in intersect(c("offset", "weights"), names(args))) {
    a <- args[[v]]
    if (is.null(a) || (is.character(a) && length(a) == 1L)) next
    if (length(a) != n)
      stop("ame(B > 0) refits the model on resampled rows of 'newdata': the fit's '", v, "' must be a column of ",
           "the data or one value per row of 'newdata'.", call. = FALSE)
    rowvec <- c(rowvec, v)
  }
  ## the fit's cluster when it has one (its units relabelled within each drawn cluster, as cpb_fe()'s own bootstrap
  ## does), else its first fixed-effects column, else rows
  units <- intersect(.ud_unit_vars(object), names(newdata))
  cv <- if (length(units)) units[1L]; clv <- NULL; lab <- "obs"
  if (is.character(ccl) && length(ccl) == 1L) { clv <- newdata[[ccl]]; lab <- ccl }
  else if (!is.null(ccl) && length(ccl) == n) { clv <- ccl; lab <- "custom" }
  else if (!is.null(cv)) { clv <- newdata[[cv]]; lab <- cv }
  if (is.factor(clv)) clv <- droplevels(clv)                 # an unused level would be drawn as an empty cluster
  off_rows <- !is.null(offset) && !is.character(offset) && length(offset) == n && n > 1L
  D <- .ud_boot(n, clv, function(idx, blk) {
    d <- newdata[idx, , drop = FALSE]
    if (!is.null(cv)) d[[cv]] <- paste0(blk, "_", as.character(d[[cv]]))
    a <- args; for (v in rowvec) a[[v]] <- a[[v]][idx]
    f2 <- .ud_quiet_guard(.ud_replicate(do.call(fun, c(list(data = d, se = "none"), a))))
    if (isFALSE(f2$converged)) return(NULL)
    dk <- if (is.null(offset)) list() else list(offset = if (off_rows) offset[idx] else offset)
    v <- est(f2, d, wt[idx], dk)
    if (is.finite(v)) v else NULL
  }, "ame", B, cores, label = lab)
  ok <- D[is.finite(D[, 1L]), 1L]
  if (length(ok) < 2L) return(list(lower = NA_real_, upper = NA_real_, method = "none"))
  a <- (1 - level) / 2
  list(lower = stats::quantile(ok, a, names = FALSE), upper = stats::quantile(ok, 1 - a, names = FALSE),
       method = sprintf("bootstrap (refit, B=%d)", B))
}

#' @export
print.ud_ame <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1)) & names(y) != "n"
  y[num] <- lapply(y[num], round, digits)
  if (all(is.na(y$lower))) y$lower <- y$upper <- NULL
  if (all(y$method == "none")) y$method <- NULL
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}
