## ---------------------------------------------------------------------------
## Zero-inflated continuous parameter binomial (mixture variant), fit by EM.
## Use when a zero is ambiguous between structural nonparticipation and a
## sampling zero from a low-rate count; for a pure participation process where
## participants have Y >= 1 by construction, prefer hurdle_cpb().
## ---------------------------------------------------------------------------

#' Simulate from the zero-inflated continuous parameter binomial
#'
#' @param n Number of units.
#' @param lambda Intensity rate (the CPB's `lambda`); scalar or length-`n`.
#' @param alpha CPB shape parameter in (0, 1).
#' @param pi Structural-zero probability; scalar or length-`n`.
#' @param support The CPB's support, `"ceiling"` (King's CPB, the default) or
#'   `"floor"`; see [cpb()].
#' @return An integer vector: a structural zero with probability `pi`, otherwise
#'   an (untruncated) CPB draw, which may itself be zero.
#' @examples
#' set.seed(1)
#' y <- rzicpb(1000, lambda = 3, alpha = 0.5, pi = 0.3)
#' mean(y == 0)
#' @export
rzicpb <- function(n, lambda, alpha, pi, support = c("ceiling", "floor")) {
  support <- match.arg(support)
  s <- stats::rbinom(n, 1, pi)
  y <- integer(n); idx <- which(s == 0)
  if (length(idx)) y[idx] <- rcpb(length(idx), rep_len(lambda, n)[idx], alpha, truncated = FALSE, support = support)
  y
}

#' Zero-inflated continuous parameter binomial regression
#'
#' Fits the zero-inflated CPB. The model mixes a structural-zero process
#' (a logistic model for the probability `pi` that a unit is a structural zero)
#' with an untruncated CPB for the count, so a zero can arise either structurally
#' or as a sampling zero from the CPB. Unlike [hurdle_cpb()], the likelihood does
#' not factorize, so the parameters are estimated by direct joint maximum
#' likelihood. The mixture likelihood's optima differ with the structural share
#' of zeros, so the search starts from points that span it: separate fits (a
#' zero-truncated CPB on the positives and a logit of the zero indicator), the
#' zero-inflated Poisson's coefficients, and a low-inflation start (a constant
#' structural-zero probability of 0.02), followed by a BFGS pass; it is
#' cross-checked against an expectation-maximization climb of the same
#' observed-data likelihood, keeping the better optimum. The separation check
#' (Separation) restarts the final polish from any better point it finds. An
#' EM-only path (`method = "em"`) is retained for comparison.
#'
#' The inflation equation uses a logit link. For a probit or cloglog inflation
#' link, use [zi_count()] (Poisson/NB/COM-Poisson count), whose `link` argument
#' covers the binary stage.
#'
#' @param formula Intensity (count) model formula (`y ~ x`).
#' @param data A data frame.
#' @param zero Optional one-sided formula (`~ z`) for the zero-inflation model;
#'   defaults to the intensity model's right-hand side.
#' @param fe Optional column name for unit fixed effects in the intensity (count)
#'   component. Under `method = "ml"` the unit effects enter as plug-in offsets
#'   from a first-stage fit to the positive counts, and on King's support the
#'   fit is then polished on the observed-data likelihood with every unit's
#'   effect free, so that it is a joint maximum of the mixture likelihood (on
#'   the floor support the fit keeps the plug-in effects, which are not a joint
#'   maximum: a unit whose positive counts are all 1 takes an arbitrary point of
#'   the first stage's flat set, so the inflation equation and the
#'   log-likelihood can depend on the order of the rows); under
#'   `method = "em"` a classification (hard-assignment) step handles the
#'   structural zeros before the same polish. A unit with no
#'   positive count is left out of the estimation, as fixest leaves such units
#'   out of a Poisson fixed-effects fit: at the maximum its count effect is minus
#'   infinity, so its zeros are certain whatever its structural-zero
#'   probability, and its rows add nothing to the likelihood or the degrees of
#'   freedom. Its rows stay in the fit at that limit (the effect at the
#'   [cpb_fe()] sentinel, 30 log units below the unit's largest linear
#'   predictor), with fitted means of zero, and the fit records
#'   the units in `zero_units` and their number in `n_zero_units`. Their
#'   structural-zero probability is `NA` where the zero equation would give
#'   them their own dummy (`zero_fe`), which their zeros cannot identify. On the
#'   floor support such units stay in the estimation with a filled-in count
#'   effect (the smallest estimated effect minus 2 under `"ml"`, the mean effect
#'   under `"em"`).
#' @param zero_fe Optional column name for fixed effects in the zero-inflation
#'   equation, entered as factor dummies (opt-in). On King's support the dummy
#'   of a unit with no positive count is not counted in the degrees of freedom
#'   (see [compare_models()], Degrees of freedom), and in new data such a unit's
#'   rows get the structural-zero probability one and those of a unit whose
#'   every count is positive (`zero_fe_one_units`) zero; a unit the fit has not
#'   seen, or every row when `newdata` lacks the column, enters at the zero
#'   equation's reference unit.
#' @param method Estimation method: `"ml"` (default) is direct joint maximum
#'   likelihood from the starts described above; `"em"` is
#'   expectation-maximization, retained for comparison but prone to a degenerate
#'   structural-zero-probability collapse under heavy zero-inflation.
#' @param se Coefficient inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Optional cluster for the bootstrap (a column name or vector);
#'   under fixed effects the units are resampled by default. The mixture EM makes
#'   the bootstrap costly, so keep `B` modest.
#' @param max.support Guard on the evaluated support of the count component; see
#'   [cpb()]. The default `NULL` sets `max(500, 10 * max(y))` (500 on the floor
#'   support). It is kept with the fit for [simulate()].
#' @param maxit Optimizer iteration budget (scaled internally for the joint
#'   maximization; also the EM iteration cap under `method = "em"`).
#' @param offset Optional exposure offset (log scale) for the count component:
#'   a numeric vector or the name of a column in `data`, making the count a
#'   rate model. Supported for `method = "ml"` without `fe`; with `fe`, an
#'   offset is refused with an error. The bootstrap carries the offset through
#'   every replicate.
#' @param tol Relative convergence tolerance on the observed-data log-likelihood.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [cpb()]. Supported for `method = "ml"`.
#' @param cores Worker processes for the bootstrap; see [cpb()].
#' @param support The CPB's support, `"ceiling"` (King's CPB, the default) or
#'   `"floor"`; see [cpb()].
#' @section Separation:
#' A zero contributes `log(pi + (1 - pi) f0)` to the likelihood, which rises
#' with `pi`, and a positive count contributes `log(1 - pi) + log f(y)`, which
#' falls with it. So a coefficient of the zero-inflation equation can grow without
#' bound, and the likelihood with it, under three conditions. The zeros it pushes
#' toward `pi = 1` can go there freely, and so can the positives it pushes toward
#' `pi = 0`. The count component's probability of a zero, `f0`, covers the zeros
#' it pushes toward `pi = 0`. Unlike separation in a logistic regression, this
#' needs no covariate that splits the zeros from the positives. So a
#' [hurdle_cpb()] whose participation logit uses the same covariates can have a
#' finite maximum while the mixture has none. The optimizer then stops somewhere
#' along the divergent path. The zero equation's coefficients there are not
#' estimates, and the log-likelihood approximates a supremum, so likelihood
#' comparisons with the fit (likelihood ratio, AIC, Vuong) lack their usual
#' reference distributions.
#'
#' After the fit, the profile likelihood of each covariate coefficient in the
#' zero equation is evaluated at the estimate. It is evaluated again at 10 and at
#' 100 standard deviations of its covariate further from zero. Each time, the
#' other zero-equation coefficients are re-maximized with the count component
#' held at its estimate. The fit warns when a coefficient's profile 100 standard
#' deviations out is less than `qchisq(0.95, 1) / 2` below its value at the
#' estimate. That coefficient's 95% profile interval then reaches past an effect
#' of 100 on the log-odds per standard deviation, and the data put no useful
#' bound on it. This covers a coefficient with no finite estimate and one whose
#' finite maximum lies on a near-flat ridge, which a fit cannot tell apart.
#' Holding the count component fixed can only lower the profile, so the full
#' model's interval is at least as wide. Every profile can also reach the
#' no-inflation boundary, where each structural-zero probability is zero and the
#' zero equation drops out. A coefficient whose profile at the estimate lies
#' within `qchisq(0.95, 1) / 2` of that limit is unbounded because the data show
#' little inflation, not by separation, so it is not flagged, and the fit is
#' reported as near the boundary instead ([zi_test()] tests it). The
#' intercept and the `zero_fe` dummies are not checked. Their limits are the
#' no-inflation boundary ([zi_test()]) and the boundaries of single units (a unit
#' with no positive count, whose structural-zero probability goes to one, as in a
#' fixed-effects logit; with `fe`, such a unit is left out of the estimation,
#' see `fe`). The component `separation` of the fit holds:
#' - `diverging`: the terms the data do not bound;
#' - `profile`: the profile steps;
#' - `short_of_max`: how far the best point the check found lies above the fit,
#'   with the count component unchanged. That point can come from the
#'   re-maximization at the estimate, when the fit is not a stationary point of
#'   the zero equation, or from a profile that rises above its value at the
#'   estimate, since the mixture likelihood is not concave and a profile can rise
#'   before it falls. A positive value means the optimizer stopped at least that
#'   far below the maximum;
#' - `best`: the zero-equation coefficients at that point (`NULL` when there is
#'   none);
#' - `restarts`, `restart_gain`: when the check finds a point at least 0.01
#'   higher, the fit restarts its final polish there, with every parameter it
#'   estimates free, keeps the result only if it gains, and checks again, up to
#'   three times. These record the restarts kept and their total gain in
#'   log-likelihood. The other components describe the final fit;
#' - `no_inflation`, `inflation_gain`: whether the best point lies within
#'   `qchisq(0.95, 1) / 2` of the no-inflation limit (the count component held),
#'   and its distance from it;
#' - `n_extreme`: the number of observations whose structural-zero probability is
#'   numerically 0 or 1, outside units that lie wholly at one edge.
#'
#' [zi_gec()] and [zi_count()] run the same check.
#' @return An object of class `"zi_cpb"`. Its component `separation` holds the
#'   separation check (see Separation).
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 800; x <- rnorm(n); z <- rnorm(n)
#' y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
#' zi_cpb(y ~ x, data = data.frame(y = y, x = x, z = z), zero = ~ z)
#' }
#' @export
zi_cpb <- function(formula, data, zero = NULL, fe = NULL, zero_fe = NULL, method = c("ml", "em"),
                   se = c("none", "bootstrap"), B = 500, cluster = NULL, max.support = NULL,
                   maxit = 200, tol = 1e-6, offset = NULL, weights = NULL, cores = 1L,
                   support = c("ceiling", "floor")) {
  .ud_no_formula_offset(formula, zero)
  support <- match.arg(support)
  refit_args <- list(fe = fe, zero_fe = zero_fe, method = match.arg(method), max.support = max.support, maxit = maxit,
                     tol = tol, offset = offset, weights = weights, support = support,
                     cluster = if (is.character(cluster) && length(cluster) == 1L) cluster)
  data <- .ud_drop_na_fe(data, list(fe, zero_fe))
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  se <- match.arg(se); method <- match.arg(method)
  .ud_check_whole(max.support, "max.support", null_ok = TRUE)
  if (se == "bootstrap") .ud_check_whole(B, "B", lower = 2)
  if (!is.null(weights)) {
    if (identical(method, "em")) stop("'weights' is supported for method = \"ml\" only.")
    if (!(is.character(weights) && length(weights) == 1L)) {
      if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
      data[[".zi_w"]] <- as.numeric(weights); weights <- ".zi_w"   # vector -> column, subsets with data
    }
  }
  offset_name <- .ud_offset_name(offset)                        # before a vector becomes a column
  if (!is.null(offset)) {
    if (!is.null(fe))
      stop("'offset' with 'fe' is not yet supported for zi_cpb (the fixed effects enter as ",
           "plug-in offsets of their own); fit without 'fe' or fold the exposure into the data.")
    if (identical(method, "em"))
      stop("'offset' is supported for method = \"ml\" only.")
    if (!(is.character(offset) && length(offset) == 1L)) {
      if (length(offset) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
      data[[".zi_off"]] <- as.numeric(offset); offset <- ".zi_off"   # vector -> column, so it subsets with data
    }
  }
  if (!is.null(cluster) && se != "bootstrap")
    warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
  ## reduce to complete cases on every model variable so the count design, the
  ## inflation design, the offset, and the weights stay aligned
  zrhs0 <- if (is.null(zero)) formula[-2L] else zero
  mv <- unique(c(all.vars(formula), all.vars(zrhs0), fe, zero_fe, offset, weights))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  cl_vals <- .ud_cluster_values(cluster, data, keep)           # checked before the rows are reduced
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) cluster <- cl_vals
  data <- data[keep, , drop = FALSE]
  if (is.null(max.support))                                     # the replicates and refits keep the fit's guard
    max.support <- if (.cpb_fl(support)) 500L
                   else max(500L, 10L * max(stats::model.response(stats::model.frame(formula, data)), 0L))
  refit_args$max.support <- max.support
  res <- .zi_cpb_full(formula, data, zero, fe, zero_fe, max.support, maxit, tol, method, offset, weights, support)
  if (isTRUE(res$n_zero_units > 0L))
    message("zi_cpb: ", .zi_zero_units_text(res), " left out of the estimation: their count effect is at its limit, ",
            "minus infinity, so their zeros are certain and add nothing to the likelihood.")
  res$call <- match.call(); res$refit_args <- refit_args; res$offset_name <- offset_name
  .cpb_warn_limit(res$alpha, support)
  .zi_warn_separation(res$separation, "zi_cpb")
  res$se.type <- se
  res$se.beta <- stats::setNames(rep(NA_real_, length(res$coefficients)), names(res$coefficients))
  res$se.zero <- stats::setNames(rep(NA_real_, length(res$zero_coef)), names(res$zero_coef))
  if (se == "bootstrap")
    res <- .zi_cpb_boot(res, formula, data, zero, fe, zero_fe, cluster, B, max.support, maxit, tol,
                        method, offset, weights, cores, support = support)
  res
}

## The fit, its separation check, and the restart of its final polish from any better point the check finds (see
## Separation); the bootstrap replicates run the same. With `fe` on King's support, the units with no positive count
## are left out of the estimation and their rows restored at the limit afterwards (.zi_add_zero_units).
.zi_cpb_full <- function(formula, data, zero, fe, zero_fe, max.support, maxit, tol, method, offset, weights, support) {
  zu <- .zi_zero_units(formula, data, fe, support)
  if (length(zu)) {
    est <- !(as.character(data[[fe]]) %in% zu)
    res <- .zi_cpb_full(formula, data[est, , drop = FALSE], zero, fe, zero_fe, max.support, maxit, tol, method,
                        offset, weights, support)
    return(.zi_add_zero_units(res, data, est, zu, zero_fe, weights))
  }
  res <- .zi_cpb_fit(formula, data, zero, fe, zero_fe, max.support, maxit, tol, method, offset, weights,
                     support = support)
  chk <- .zi_restart(res,
    check = function(f) .zi_separation(f$Z, f$zero_coef, f$y, f$p0_full, f$weights, f$zero_group,
                                       .zi_fe_cols(f$Z, zero_fe)),
    refit = function(f, g) .zi_cpb_fit(formula, data, zero, fe, zero_fe, max.support, maxit, tol, method, offset,
                                       weights, support = support,
                                       restart = list(b = f$coefficients, alpha = f$alpha, g = g, fe_hat = f$fe_hat,
                                                      estimated = f$fe_estimated, iterations = f$iterations)),
    loglik = function(f) f$loglik)
  res <- chk$fit; res$separation <- chk$separation
  ## on King's support, a zero_fe level with no positive count among the rows fitted has its dummy at the limit: not
  ## counted (.ud_zero_levels); the units left out above have no rows here, so none is counted twice
  zl <- if (is.null(zero_fe) || .cpb_fl(support)) character(0)
        else if (length(zero_fe) == 1L) .ud_zero_levels(res$y, res$zero_group)
        else .ud_fe_zero_levels(res$y, data, zero_fe)              # res rows are in the data's order
  res$df <- res$df - length(zl); res$zero_fe_zero_units <- zl
  ## and one whose every count is positive has it at the other limit (counted; predicted there in new data)
  res$zero_fe_one_units <- if (is.null(zero_fe) || .cpb_fl(support) || length(zero_fe) != 1L) character(0)
                           else .ud_one_levels(res$y, res$zero_group)
  res
}

## The units of a fixed-effects fit on King's support whose counts are all zero (none when every unit's are: the fit
## then warns that every count is zero). The floor support keeps the fill values of versions before 0.2.0.
.zi_zero_units <- function(formula, data, fe, support) {
  if (is.null(fe) || .cpb_fl(support)) return(character(0))
  y <- stats::model.response(stats::model.frame(formula, data))
  pos <- tapply(y > 0, as.character(data[[fe]]), any)
  if (!any(pos)) character(0) else names(pos)[!pos]
}

## The rows of the units .zi_zero_units() left out, restored to a fit made without them (`est` marks the rows it
## used): their count effect at the cpb_fe() sentinel, 30 log units below the unit's largest x'b, so every rate of
## the unit is at most exp(-30) and P(0) is 1 to within about 1e-13; they add nothing to the log-likelihood or the
## degrees of freedom. Their structural-zero probability comes from
## the zero equation where it covers them, and is NA where they would need their own zero_fe dummy.
.zi_add_zero_units <- function(res, data, est, zu, zero_fe, weights) {
  n <- nrow(data); dz <- data[!est, , drop = FALSE]; m <- nrow(dz)
  put <- function(v, z) { out <- rep(NA, n); out[est] <- v; out[!est] <- z; out }
  b <- res$coefficients
  Xz <- stats::model.matrix(stats::delete.response(stats::terms(res$formula)), dz)
  xb <- if (all(names(b) %in% colnames(Xz))) as.numeric(Xz[, names(b), drop = FALSE] %*% b) else rep(0, m)
  ue <- -30 - tapply(xb, as.character(dz[[res$fe]]), max)        # each unit's sentinel effect
  lam <- exp(ue[as.character(dz[[res$fe]])] + xb)
  p0 <- cpb_lp0_cpp(c(1, stats::qlogis(res$alpha)), matrix(1, m, 1), integer(m), log(lam) - 1, res$max.support,
                    FALSE, .cpb_fl(res$support))[, 2]
  seen <- if (is.null(zero_fe)) rep(TRUE, m)
          else if (length(zero_fe) == 1L) as.character(dz[[zero_fe]]) %in% as.character(res$zero_group)
          else Reduce(`&`, lapply(zero_fe, function(f) as.character(dz[[f]]) %in% as.character(data[[f]][est])))
  Zz <- matrix(NA_real_, m, ncol(res$Z), dimnames = list(NULL, colnames(res$Z)))
  if (any(seen)) {
    Zs <- tryCatch(.ud_newdata_matrix(res$zero_terms, dz[seen, , drop = FALSE], res$zero_xlev, res$zero_contrasts,
                                      names(res$zero_coef)), error = function(e) NULL)
    if (!is.null(Zs)) Zz[seen, ] <- Zs
  }
  Zf <- matrix(NA_real_, n, ncol(res$Z), dimnames = list(NULL, colnames(res$Z))); Zf[est, ] <- res$Z; Zf[!est, ] <- Zz
  lv <- levels(factor(data[[res$fe]]))
  fe_all <- c(res$fe_hat, ue[zu])[lv]
  if (!is.null(res$fe_estimated))
    res$fe_estimated <- unname(c(stats::setNames(res$fe_estimated, names(res$fe_hat)),
                                 stats::setNames(rep(FALSE, length(zu)), zu))[lv])
  res$fe_hat <- fe_all
  res$y <- as.integer(put(res$y, 0L))
  if (!is.null(res$weights)) res$weights <- put(res$weights, as.numeric(dz[[weights]]))
  res$n <- n; res$nobs_weighted <- if (is.null(res$weights)) n else sum(res$weights)
  res$lambda_full <- put(res$lambda_full, lam); res$p0_full <- put(res$p0_full, p0)
  res$pi_full <- put(res$pi_full, as.numeric(stats::plogis(Zz %*% res$zero_coef)))
  res$ceiling <- put(res$ceiling, lam / (1 - res$alpha))
  res$Z <- Zf
  if (!is.null(zero_fe)) res$zero_group <- .ud_fe_group(data, zero_fe)
  res$zero_units <- zu; res$n_zero_units <- length(zu); res$n_zero_rows <- m; res$fit_rows <- est
  res
}
.zi_zero_units_text <- function(x)
  paste0(x$n_zero_units, if (x$n_zero_units == 1L) " unit" else " units", " (", x$n_zero_rows,
         if (x$n_zero_rows == 1L) " row" else " rows", ") with no positive count")
## A structural-zero probability that is NA for a left-out unit (see .zi_add_zero_units) is set to 0 where a
## distribution or a mean is computed: the unit's count component puts all its mass on zero, so any value gives the
## same distribution.
.zi_pi0 <- function(p) { p[is.na(p)] <- 0; p }

## Pairs or cluster bootstrap for the zero-inflated CPB (.ud_boot): rows, the clusters, or by default the
## fixed-effects units, a unit drawn twice relabeled as two units (in the zero equation's unit dummies too, when
## they are on the same variable). Each replicate refits the whole model, the separation restart included. The
## draws hold the count and zero-equation coefficients (a zero_fe fit keeps its covariate coefficients only; its
## unit dummies are not the same units across replicates), alpha, and for a fixed-effects fit the reference unit
## effect.
.zi_cpb_boot <- function(res, formula, data, zero, fe, zero_fe, cluster, B, max.support, maxit, tol,
                         method = "ml", offset = NULL, weights = NULL, cores = 1L, support) {
  use_fe <- !is.null(fe)
  cvar <- if (!is.null(cluster)) cluster else if (use_fe) fe else NULL
  cval <- if (is.null(cvar)) NULL else if (is.character(cvar) && length(cvar) == 1L) data[[cvar]] else cvar
  fev  <- if (use_fe) as.character(data[[fe]]) else NULL
  zcov <- setdiff(names(res$zero_coef), .zi_fe_cols(res$Z, zero_fe))
  cols <- c(paste0("count:", names(res$coefficients)), paste0("zero:", zcov), "alpha", if (use_fe) "fe_ref")
  relabel_zero <- use_fe && identical(zero_fe, fe)
  boot <- .ud_boot(nrow(data), cval, function(idx, blk) {
    bd <- data[idx, , drop = FALSE]; ffe <- fe; zfe <- zero_fe
    if (use_fe && !is.null(blk)) {
      bd[[".bootunit"]] <- paste0(blk, "_", fev[idx]); ffe <- ".bootunit"
      if (relabel_zero) zfe <- ".bootunit"
    }
    fb <- .ud_quiet_guard(.zi_cpb_full(formula, bd, zero, ffe, zfe, max.support, maxit, tol, method, offset,
                                       weights, support))
    if (!isTRUE(fb$converged)) return(NULL)
    c(fb$coefficients, fb$zero_coef[zcov], fb$alpha, if (use_fe) fb$int_fe_ref)
  }, cols, B, cores, label = if (is.null(cvar)) "obs" else if (is.null(cluster)) fe else .ud_cluster_label(cluster))
  res$boot <- boot
  bc <- .ud_boot_stats(boot, paste0("count:", names(res$coefficients)), res$coefficients, normal = TRUE)
  bz <- .ud_boot_stats(boot, paste0("zero:", zcov), res$zero_coef[zcov], normal = TRUE)
  res$se.beta <- bc$se; res$ci.beta <- bc$ci
  res$se.zero[zcov] <- bz$se
  res$se.type <- "bootstrap"
  res$cluster <- attr(boot, "cluster")
  res
}

## Seeds from the zero-inflated Poisson, the alpha -> 1 limit of the ZI-CPB, fit by
## EM with exact glm M-steps on the optimizer's (scaled) designs. The mixture
## likelihood can hold two optima at the same dispersion, and the separate-stage
## seeds can all lie in the basin of the lower one (on real panels, tens of
## log-likelihood points below the optimum the ZIP's coefficients lead to), so
## the multistart also runs from the ZIP's coefficients, at dispersions
## spread between the smallest alpha whose supports cover every positive count at
## the ZIP's rates (floor support: y <= lambda / (1 - alpha); King's: y - 1 <
## lambda / (1 - alpha)) and 1.
.zi_zip_seeds <- function(Xs, Zs, y, off, w, maxit = 200L, tol = 1e-8, support) {
  is0 <- y == 0
  b <- stats::glm.fit(Xs, y, weights = w, offset = off, family = stats::poisson())$coefficients
  b[!is.finite(b)] <- 0
  g <- suppressWarnings(stats::glm.fit(Zs, as.numeric(is0), weights = w, family = stats::binomial())$coefficients)
  g[!is.finite(g)] <- 0
  ll_old <- -Inf
  for (it in seq_len(maxit)) {
    lam <- as.numeric(exp(off + Xs %*% b)); pit <- as.numeric(stats::plogis(Zs %*% g))
    wr <- ifelse(is0, pit / (pit + (1 - pit) * exp(-lam)), 0)
    b <- suppressWarnings(stats::glm.fit(Xs, y, weights = (1 - wr) * w, offset = off, start = b,
                                         family = stats::poisson())$coefficients)
    g <- suppressWarnings(stats::glm.fit(Zs, wr, weights = w, start = g, family = stats::quasibinomial())$coefficients)
    if (!all(is.finite(b)) || !all(is.finite(g))) return(NULL)
    lam <- as.numeric(exp(off + Xs %*% b)); pit <- as.numeric(stats::plogis(Zs %*% g))
    ll <- sum(w * ifelse(is0, log(pmax(pit + (1 - pit) * exp(-lam), 1e-300)),
                         log(pmax(1 - pit, 1e-300)) + stats::dpois(y, lam, log = TRUE)))
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) break
    ll_old <- ll
  }
  lam <- as.numeric(exp(off + Xs %*% b))
  af <- if (.cpb_fl(support)) max(0, (1 - lam / y)[y > 0]) else max(0, (1 - lam / (y - 1))[y > 1])
  list(b = unname(b), g = unname(g), alpha = af + (1 - af) * c(0.25, 0.5, 0.75))
}

## The zero-inflated Poisson by direct maximum likelihood (BFGS with the analytic
## gradient) from the EM's ZIP (`zs`), from the Poisson regression with a
## constant structural-zero probability of 0.02, and from the Poisson regression
## with the logit of the zero indicator (pscl's start). King's zero-inflated CPB
## contains the ZIP as its limit alpha -> 1, and on some data the EM's ZIP sits
## in the high-inflation basin while the ZIP's maximum, and with it the CPB
## mixture's supremum, lies in the low-inflation one. Returns the best solution
## in the scaled basis, or NULL.
.zi_zip_ml <- function(Xs, Zs, y, off, w, zs = NULL) {
  is0 <- y == 0; pb <- ncol(Xs); pg <- ncol(Zs)
  nll <- function(p) {
    lam <- as.numeric(exp(off + Xs %*% p[seq_len(pb)])); pit <- as.numeric(stats::plogis(Zs %*% p[pb + seq_len(pg)]))
    v <- -sum(w * ifelse(is0, log(pit + (1 - pit) * exp(-lam)), log1p(-pit) + stats::dpois(y, lam, log = TRUE)))
    if (is.finite(v)) v else 1e10
  }
  gr <- function(p) {
    lam <- as.numeric(exp(off + Xs %*% p[seq_len(pb)])); pit <- as.numeric(stats::plogis(Zs %*% p[pb + seq_len(pg)]))
    L0 <- pmax(pit + (1 - pit) * exp(-lam), 1e-300)
    dc <- ifelse(is0, -(1 - pit) * exp(-lam) * lam / L0, y - lam)
    dz <- ifelse(is0, pit * (1 - pit) * (1 - exp(-lam)) / L0, -pit)
    -c(as.numeric(crossprod(Xs, w * dc)), as.numeric(crossprod(Zs, w * dz)))
  }
  b0 <- tryCatch({ v <- stats::glm.fit(Xs, y, weights = w, offset = off, family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, pb))
  icz <- match("(Intercept)", colnames(Zs))
  glow <- rep(0, pg); if (!is.na(icz)) glow[icz] <- stats::qlogis(0.02)   # the intercept column is unscaled
  gl <- suppressWarnings(stats::glm.fit(Zs, as.numeric(is0), weights = w, family = stats::binomial())$coefficients)
  gl[!is.finite(gl)] <- 0
  starts <- list(c(b0, glow), c(b0, gl))
  if (!is.null(zs)) starts <- c(list(c(zs$b, zs$g)), starts)
  best <- NULL
  for (st in starts) {
    o <- tryCatch(stats::optim(st, nll, gr, method = "BFGS", control = list(maxit = 1000, reltol = 1e-12)),
                  error = function(e) NULL)
    if (!is.null(o) && o$value < 1e9 && (is.null(best) || o$value < best$value)) best <- o
  }
  if (is.null(best)) return(NULL)
  list(b = best$par[seq_len(pb)], g = best$par[pb + seq_len(pg)], value = best$value)
}

## A low-inflation start for the mixture: the single-equation (untruncated) CPB's
## coefficients and dispersion, and a constant structural-zero probability of
## 0.02, on the original scale (NULL when the inflation design has no intercept
## or the fit fails). The mixture likelihood can peak near pi = 0 with a count
## model that alone explains most zeros, where the separate-stage and ZIP starts
## can all end lower.
.zi_low_start <- function(X, Z, y, off, w, max.support, support) {
  icz <- match("(Intercept)", colnames(Z)); if (is.na(icz)) return(NULL)
  sf <- tryCatch(.ud_quiet_guard(.cpb_fit_ms(X, y, off, max.support, FALSE, c(0.5, 0.9), 4000L, 1e-8, w = w,
                                             restarts = 0L, support = support)), error = function(e) NULL)
  if (is.null(sf) || !is.finite(sf$value) || sf$value >= 1e9) return(NULL)
  g <- rep(0, ncol(Z)); g[icz] <- stats::qlogis(0.02)
  list(b = sf$par[seq_len(ncol(X))], alpha = stats::plogis(sf$par[ncol(X) + 1L]), g = g)
}

## Observed-data polish of the fixed-effects classification EM: block coordinate
## ascent on the mixture log-likelihood (each branch's probability floored only at
## the smallest positive double, against underflow; a positive count the support
## excludes stays impossible).
## Each round maximizes, for every estimated unit, its count intercept jointly
## with its own inflation dummy when the inflation equation has one (zero_fe;
## Nelder-Mead in two dimensions, since the two compete for the unit's zeros and
## alternating between them creeps), or its count intercept alone (bounded below
## by the feasibility of its counts); then the inflation coefficients (BFGS,
## analytic gradient); then the covariate coefficients with alpha (Nelder-Mead
## with restarts, and a BFGS pass on King's support). Each block's move is kept
## only when it gains; rounds stop once one gains less than 1e-8 relative. Units
## the EM never estimated (no positive count) keep its fill value, as in its own
## likelihood. Inputs are sorted by unit; ui indexes the units; `w` holds the
## frequency weights of the rows (the direct ML path's; the EM path takes none).
.zi_fe_polish <- function(y, Xcov, Z, ui, fe_hat, b, a, g, ms, fl, estimated, maxit = 50L, w = NULL) {
  n <- length(y); is0 <- y == 0; one <- matrix(1, n, 1); fe_hat <- as.numeric(fe_hat); w <- .ud_w1(w, n)
  b <- as.numeric(b); g <- as.numeric(g); k <- length(b)
  rows <- split(seq_len(n), factor(ui, levels = seq_along(fe_hat)))
  ## the inflation column, if any, that is unit u's own dummy (one on its rows, zero elsewhere)
  own <- rep(NA_integer_, length(fe_hat))
  for (j in seq_len(ncol(Z))) {
    nz <- Z[, j] != 0; uu <- unique(ui[nz])
    if (length(uu) == 1L && all(Z[nz, j] == 1) && sum(nz) == length(rows[[uu]])) own[uu] <- j
  }
  rowll <- function(eta_c, a, eta_z, idx = seq_len(n)) {
    m <- cpb_lp0_cpp(c(1, stats::qlogis(a)), one[idx, , drop = FALSE], y[idx], eta_c - 1, ms, FALSE, fl)
    pit <- stats::plogis(eta_z)
    ## a positive count outside the support (the kernel's -1e300) stays impossible
    tiny <- .Machine$double.xmin
    ifelse(is0[idx], log(pmax(pit + (1 - pit) * m[, 2], tiny)),
           log(pmax(1 - pit, tiny)) + ifelse(m[, 1] < -1e299, -Inf, m[, 1]))
  }
  total <- function(fe_hat, b, a, g) sum(w * rowll(fe_hat[ui] + as.numeric(Xcov %*% b), a, as.numeric(Z %*% g)))
  ll <- total(fe_hat, b, a, g); it <- 0L
  for (it in seq_len(maxit)) {
    ll_old <- ll
    ## (1) unit count intercepts, each with the unit's own inflation dummy where there is one
    xb <- as.numeric(Xcov %*% b); zg <- as.numeric(Z %*% g)
    for (u in which(estimated)) {
      r <- rows[[u]]; yr <- y[r]; j <- own[u]
      if (!is.na(j)) {
        zr <- zg[r] - g[j]
        f2 <- function(p) { v <- -sum(w[r] * rowll(p[1] + xb[r], a, zr + p[2], r)); if (is.finite(v)) v else 1e10 }
        o <- stats::optim(c(fe_hat[u], g[j]), f2, method = "Nelder-Mead", control = list(reltol = 1e-12, maxit = 2000))
        if (o$value < f2(c(fe_hat[u], g[j]))) { fe_hat[u] <- o$par[1]; zg[r] <- zr + o$par[2]; g[j] <- o$par[2] }
        next
      }
      fu <- function(au) { v <- -sum(w[r] * rowll(au + xb[r], a, zg[r], r)); if (is.finite(v)) v else 1e10 }
      need <- if (fl) yr > 0 else yr > 1
      lo <- if (any(need)) max(log((if (fl) yr[need] else yr[need] - 1) * (1 - a)) - xb[r][need]) + 1e-9
            else fe_hat[u] - 6
      lo <- min(lo, fe_hat[u]); hi <- max(fe_hat[u], lo) + 3
      o <- stats::optimize(fu, c(lo, hi), tol = 1e-8)
      if (o$objective < fu(fe_hat[u])) fe_hat[u] <- o$minimum
    }
    ## (2) inflation coefficients, given the count component's P(0)
    eta_c <- fe_hat[ui] + xb
    p0 <- cpb_lp0_cpp(c(1, stats::qlogis(a)), one, y, eta_c - 1, ms, FALSE, fl)[, 2]
    fg <- function(gg) { pit <- stats::plogis(as.numeric(Z %*% gg)); tiny <- .Machine$double.xmin
      -sum(w * ifelse(is0, log(pmax(pit + (1 - pit) * p0, tiny)), log(pmax(1 - pit, tiny)))) }
    gr <- function(gg) { pit <- stats::plogis(as.numeric(Z %*% gg))
      d <- ifelse(is0, pit * (1 - pit) * (1 - p0) / pmax(pit + (1 - pit) * p0, .Machine$double.xmin), -pit)
      -as.numeric(crossprod(Z, w * d)) }
    og <- tryCatch(stats::optim(g, fg, gr, method = "BFGS", control = list(maxit = 500, reltol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(og) && og$value < fg(g)) g <- og$par
    ## (3) covariate coefficients and alpha
    zg <- as.numeric(Z %*% g)
    fb <- function(p) { v <- -sum(w * rowll(fe_hat[ui] + as.numeric(Xcov %*% p[seq_len(k)]), stats::plogis(p[k + 1L]), zg))
      if (is.finite(v)) v else 1e10 }
    p0b <- c(b, stats::qlogis(a)); v0 <- fb(p0b)
    o <- stats::optim(p0b, fb, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
    for (rs in 1:3) { o2 <- stats::optim(o$par, fb, method = "Nelder-Mead", control = list(maxit = 4000, reltol = 1e-12))
      if (o2$value > o$value - 1e-8) break; o <- o2 }
    if (!fl) o <- .ud_bfgs_pass(o, fb)
    if (o$value < v0) { b <- o$par[seq_len(k)]; a <- stats::plogis(o$par[k + 1L]) }
    ll <- total(fe_hat, b, a, g)
    if (ll - ll_old < 1e-8 * (abs(ll_old) + 1)) break
  }
  list(fe_hat = fe_hat, b = b, a = a, g = g, ll = ll, rounds = it)
}

## `restart`: a fit's coefficients (original units) with its zero-equation coefficients replaced, from which only
## the final polish runs (the restart from the separation check, .zi_restart): the observed-data polish on the
## fixed-effects EM path, the Nelder-Mead restarts and BFGS pass of the direct ML otherwise.
.zi_cpb_fit <- function(formula, data, zero = NULL, fe = NULL, zero_fe = NULL, max.support = 500, maxit = 200,
                        tol = 1e-6, method = "ml", offset = NULL, weights = NULL, support, restart = NULL) {
  fl <- .cpb_fl(support)
  mf <- stats::model.frame(formula, data, drop.unused.levels = TRUE)
  y  <- stats::model.response(mf)
  if (!is.numeric(y)) stop("Response must be a numeric count; got ", class(y)[1L], ".")
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  .ud_warn_all_zero(y)
  y  <- as.integer(y)
  X  <- stats::model.matrix(formula, data); .ud_rank_check(X, "count design")
  w  <- .ud_w1(.ud_weights(weights, data, seq_len(nrow(data))), nrow(X))
  zrhs <- if (is.null(zero)) formula[-2L] else zero
  if (!is.null(zero_fe))                                          # fixed effects in the zero equation
    zrhs <- stats::reformulate(c(labels(stats::terms(zrhs)), paste0("factor(", zero_fe, ")")))
  Zt <- stats::terms(stats::update(zrhs, ~ .))
  Z  <- stats::model.matrix(Zt, data); .ud_rank_check(Z, "zero design")
  int_terms0 <- stats::delete.response(stats::terms(formula))            # stored factor levels so
  int_xlev  <- stats::.getXlevels(int_terms0, mf)                        # predict(newdata=) does not
  zero_xlev <- stats::.getXlevels(Zt, stats::model.frame(Zt, data, drop.unused.levels = TRUE))      # drop levels absent from newdata
  int_contrasts <- attr(X, "contrasts"); zero_contrasts <- attr(Z, "contrasts")   # the fit's factor coding
  n  <- length(y); kmax <- max(y); is0 <- y == 0

  ## fixed-effects path. The default ("ml") starts as a TWO-STEP estimator: the unit
  ## count intercepts are taken from a first-stage zero-truncated FE-CPB fit to the
  ## positive counts and held fixed as an offset, after which the observed-data
  ## mixture likelihood is maximized over (beta, alpha, gamma), the inflation
  ## equation carrying the terms in `zero` (and `zero_fe`'s unit dummies, which the
  ## outer optimizer fits directly). On King's support that plug-in fit then starts
  ## the observed-data polish (.zi_fe_polish), which frees every unit's count effect
  ## with the other parameters, so the fit is a joint maximum that uses the zeros;
  ## the floor support keeps the plug-in, as before 0.2.0.
  ## method = "em" is a seeded classification-EM that re-estimates the unit count
  ## effects from the FULL data (positives and sampling zeros) and fits the
  ## inflation fixed effects (`zero_fe`) cheaply via glm at each step, then polishes
  ## on the observed-data likelihood; seeding it from the separate-stage fits avoids
  ## the degenerate pi -> 0 basin that an unseeded EM can fall into on heavily-zero
  ## panels.
  if (!is.null(fe)) {
    ord  <- order(data[[fe]]); data <- data[ord, , drop = FALSE]
    y    <- y[ord]; is0 <- is0[ord]; Z <- Z[ord, , drop = FALSE]; w <- w[ord]
    Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
    Xcov <- Xall[, setdiff(colnames(Xall), "(Intercept)"), drop = FALSE]
    uf   <- factor(data[[fe]]); nu <- nlevels(uf)
    ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
    pcov <- ncol(Xcov); pg <- ncol(Z); ms <- as.integer(max.support)
    cf <- if (is.null(restart))
            tryCatch(.ud_quiet_guard(cpb_fe(formula, data = data[y > 0, , drop = FALSE], fe = fe, truncated = TRUE,
                                            max.support = max.support, weights = weights, support = support)),
                     error = function(e) stop("zi_cpb (fe): the first-stage fixed-effects fit to the positive ",
                                              "counts failed, so the unit effects cannot be estimated: ",
                                              conditionMessage(e), call. = FALSE))
    b0 <- if (!is.null(cf)) cf$coefficients[colnames(Xcov)] else stats::setNames(rep(0, pcov), colnames(Xcov))
    b0[!is.finite(b0)] <- 0
    a0 <- if (!is.null(cf)) cf$alpha else 0.5
    g0 <- suppressWarnings(stats::glm.fit(Z, as.numeric(is0), family = stats::binomial())$coefficients)
    g0[!is.finite(g0)] <- 0

    if (method == "em" && !is.null(restart)) {
      b <- restart$b; a <- restart$alpha; g <- restart$g; fe_hat <- restart$fe_hat; estimated <- restart$estimated
      fe_df <- sum(estimated); iters <- restart$iterations
      pol <- .zi_fe_polish(y, Xcov, Z, as.integer(uf), fe_hat, b, a, g, as.integer(max.support), fl, estimated)
      fe_hat[] <- pol$fe_hat; b[] <- pol$b; a <- pol$a; g[] <- pol$g; ll <- pol$ll
      conv_fe <- fl || pol$rounds < 50L                            # the polish stopped before its round cap
    } else if (method == "em") {
      b <- b0; a <- a0; g <- g0; fe_hat <- if (!is.null(cf)) cf$fe else stats::setNames(rep(0, nu), levels(uf))
      lam_fun <- function() { fi <- fe_hat[as.character(data[[fe]])]; fi[is.na(fi)] <- mean(fe_hat)
        as.numeric(exp(fi + Xcov %*% b)) }
      p0_fun <- function(lam) vapply(lam, function(l) .cpb_pmf1(l, a, kmax, FALSE, support = support)[1L], numeric(1))
      ll_old <- -Inf; s_old <- rep(NA, length(y)); best <- NULL; iters <- 0L; em_done <- FALSE
      for (it in seq_len(maxit)) {
        iters <- it; lam <- lam_fun(); pit <- as.numeric(stats::plogis(Z %*% g))
        w <- ifelse(is0, pit / (pit + (1 - pit) * p0_fun(lam)), 0)
        g <- suppressWarnings(stats::glm.fit(Z, w, family = stats::quasibinomial()))$coefficients
        g[!is.finite(g)] <- 0
        s <- is0 & (w > 0.5); if (it > 1L && all(s == s_old)) { em_done <- TRUE; break }; s_old <- s
        cf2 <- tryCatch(.ud_quiet_guard(cpb_fe(formula, data = data[!s, , drop = FALSE], fe = fe, truncated = FALSE,
                                               max.support = max.support, maxit = 800L, inner_it = 20L,
                                               support = support)),
                        error = function(e) NULL)
        if (!is.null(cf2)) { b <- cf2$coefficients[colnames(Xcov)]; b[!is.finite(b)] <- 0; a <- cf2$alpha; fe_hat <- cf2$fe }
        lam <- lam_fun(); pit <- as.numeric(stats::plogis(Z %*% g))
        py <- vapply(seq_along(y), function(i) .cpb_pmf1(lam[i], a, kmax, FALSE, support = support)[y[i] + 1L], numeric(1))
        ## Both mixture branches floored at 1e-12 (the score() convention): an early
        ## iterate can put pit ~ 0 and P(0) ~ 0 on a zero observation (feasibility
        ## leaves the count component no zero mass), and an unfloored log() then
        ## sends ll to -Inf, which turns the convergence test into if (NaN).
        ll <- sum(ifelse(is0, log(pmax(pit + (1 - pit) * p0_fun(lam), 1e-12)),
                         log(pmax(1 - pit, 1e-12)) + log(pmax(py, 1e-12))))
        if (is.finite(ll) && (is.null(best) || ll > best$ll)) best <- list(ll = ll, b = b, a = a, fe_hat = fe_hat, g = g)
        if (is.finite(ll) && is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) { em_done <- TRUE; break }
        ll_old <- ll
      }
      if (is.null(best)) stop("zi_cpb (EM): the log-likelihood was non-finite at every iteration; ",
                              "the count component may not cover the data (try method = \"ml\").")
      b <- best$b; a <- best$a; fe_hat <- best$fe_hat; g <- best$g; ll <- best$ll; conv_fe <- TRUE
      ## Units whose rows are all classified structural (e.g. units with no
      ## positive counts) never enter the count fit, so fe_hat covers only the
      ## units it was estimated on. Expand to all levels with the mean-effect
      ## fallback lam_fun used throughout the iterations, and count only the
      ## genuinely estimated unit effects toward the degrees of freedom.
      fe_df <- length(fe_hat)
      fe_full <- stats::setNames(rep(NA_real_, nu), levels(uf))
      fe_full[names(fe_hat)] <- fe_hat
      fe_full[is.na(fe_full)] <- mean(fe_hat)
      fe_hat <- fe_full
      ## the classification EM maximizes the classification likelihood, not the
      ## mixture likelihood the fit reports (Celeux and Govaert 1992): its
      ## solution is polished on the observed-data likelihood (.zi_fe_polish)
      estimated <- levels(uf) %in% names(best$fe_hat)
      pol <- .zi_fe_polish(y, Xcov, Z, as.integer(uf), fe_hat, b, a, g, as.integer(max.support), fl, estimated)
      fe_hat[] <- pol$fe_hat; b[] <- pol$b; a <- pol$a; g[] <- pol$g; ll <- pol$ll
      conv_fe <- fl || (em_done && pol$rounds < 50L)
    } else {
      ## Two-step robust ML. PLUG IN the unit count intercepts from a first-stage
      ## zero-truncated FE-CPB fit to the positives (held fixed as an offset -- this
      ## is a plug-in, not a within-mixture concentration, so the offset does not use
      ## the zeros), then jointly maximize the observed-data mixture likelihood over
      ## (beta, alpha, gamma) with that offset fixed. The likelihood is the vectorized
      ## non-FE mixture nll (cpb_lp0_cpp with the offset folded in as a unit-coefficient
      ## column), so this is fast. On King's support the polish below then frees the
      ## unit effects, so they use the zeros; on the floor support they stay at the
      ## plug-in (method = "em" re-estimates them there).
      ## On King's support the second stage also frees a common shift of the plug-in
      ## effects (k0 = 1): the first stage's effects absorb x-bar times ITS slopes, so
      ## without the shift a second-stage slope that differs from the first-stage one
      ## makes the fit depend on where the covariates are centred. The floor support
      ## keeps the fit of versions before 0.2.0 (k0 = 0).
      k0 <- if (fl) 0L else 1L
      off_lvl <- if (!is.null(restart)) restart$fe_hat
                 else if (!is.null(cf)) cf$fe[levels(uf)] else stats::setNames(rep(0, nu), levels(uf))
      off_lvl[!is.finite(off_lvl)] <- NA
      off_lvl[is.na(off_lvl)] <- min(off_lvl, na.rm = TRUE) - 2   # all-zero / unfit units -> ~zero rate
      ## the optimizer works on covariates scaled to unit (weighted) standard deviation
      sx <- .ud_colscale(Xcov, w); sz <- .ud_colscale(Z, w); Zs <- sweep(Z, 2, sz, "/")
      Xoff <- cbind(offset = off_lvl[as.integer(uf)], shift = if (k0) 1, sweep(Xcov, 2, sx, "/"))   # offset column, coefficient fixed at 1
      kb <- k0 + pcov                                           # the shift and the slopes
      nll <- function(par) {
        m   <- cpb_lp0_cpp(c(1, par[seq_len(kb)], par[kb + 1]), Xoff, y, rep(0, nrow(Xoff)), ms, FALSE, fl)
        pit <- as.numeric(stats::plogis(Zs %*% par[(kb + 2):(kb + 1 + pg)]))
        v <- -sum(w * ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
        if (is.finite(v)) v else 1e10
      }
      best <- NULL
      zs <- if (is.null(restart))
              tryCatch(.zi_zip_seeds(Xoff[, -1L, drop = FALSE], Zs, y, Xoff[, 1L], w, support = support), error = function(e) NULL)
      starts <- if (!is.null(restart))
                  list(c(if (k0) 0, restart$b * sx, stats::qlogis(min(max(restart$alpha, 1e-6), 1 - 1e-6)), restart$g * sz))
                else c(lapply(unique(c(a0, 0.3, 0.6)), function(as0) c(if (k0) 0, b0 * sx, stats::qlogis(min(max(as0, 0.05), 0.95)), g0 * sz)),
                       if (!is.null(zs)) lapply(zs$alpha, function(as0) c(zs$b, stats::qlogis(as0), zs$g)))
      for (st in starts) {
        op <- tryCatch(stats::optim(st, nll, method = "Nelder-Mead", control = list(maxit = 20L * maxit, reltol = tol)),
                       error = function(e) NULL)
        if (!is.null(op) && op$value < 1e9 && (is.null(best) || op$value < best$value)) best <- op
      }
      if (is.null(best)) stop("zi_cpb (fe): all optimization starts failed.")
      best <- .ud_bfgs_pass(best, nll)
      b <- best$par[k0 + seq_len(pcov)] / sx; a <- stats::plogis(best$par[kb + 1]); g <- best$par[(kb + 2):(kb + 1 + pg)] / sz
      ll <- -best$value; iters <- as.integer(best$counts[1]); conv_fe <- isTRUE(best$convergence == 0)
      fe_hat <- off_lvl + (if (k0) best$par[1] else 0)      # the plug-in effects with the common shift
      ## every level of the fit: on King's support the units with no positive count were left out before
      ## (.zi_zero_units), and the floor support counts its fill levels, as versions before 0.2.0 did
      fe_df <- length(fe_hat)
      ## On King's support the plug-in fit then starts the observed-data polish of the EM path, which frees every
      ## unit's count effect with the other parameters, so the fit is a joint maximum of the mixture likelihood.
      ## The plug-in alone is not: the first stage sees only the positive counts, and a unit whose positive
      ## counts are all 1 gets an arbitrary point of the zero-truncated CPB's flat set; its log-likelihood, which
      ## compare_models() ranks, could lie tens of points below the maximum. The floor support keeps the plug-in.
      if (!fl) {
        pol <- .zi_fe_polish(y, Xcov, Z, as.integer(uf), fe_hat, b, a, g, as.integer(max.support), fl,
                             as.logical(tapply(y > 0, uf, any)), w = w)
        fe_hat[] <- pol$fe_hat; b[] <- pol$b; a <- pol$a; g[] <- pol$g; ll <- pol$ll
        conv_fe <- conv_fe && pol$rounds < 50L
      }
    }
    names(b) <- colnames(Xcov); names(g) <- colnames(Z); names(fe_hat) <- levels(uf)
    fi <- fe_hat[as.character(data[[fe]])]; fi[is.na(fi)] <- mean(fe_hat)
    lambda_full <- as.numeric(exp(fi + Xcov %*% b))
    p0_full <- cpb_lp0_cpp(c(1, stats::qlogis(a)), matrix(1, length(y), 1), y, log(lambda_full) - 1, ms, FALSE, fl)[, 2]
    return(.ud_unsort(structure(list(coefficients = b, alpha = a, zero_coef = g, zero.coefficients = g, loglik = ll, iterations = iters,
        df = length(b) + fe_df + 1L + length(g), n = length(y), y = y, converged = conv_fe,
        weights = if (is.null(weights)) NULL else w, nobs_weighted = if (is.null(weights)) length(y) else sum(w),
        pi_full = as.numeric(stats::plogis(Z %*% g)), lambda_full = lambda_full, p0_full = p0_full,
        zero_group = .ud_fe_group(data, zero_fe), fe_estimated = if (method == "em") estimated,
        ceiling = lambda_full / (1 - a), fe = fe, fe_hat = fe_hat, Z = Z, support = support,
        int_beta = b, int_xref = colMeans(Xcov),
        int_fe_ref = .fe_mean_identified(fe_hat, if (method == "em") estimated else tapply(y, uf, max) > 0),
        max.support = ms,
        zero_beta = g, zero_xref = .fe_dummy_ref(colMeans(Z), Z, y, zero_fe, .ident_binary),
        int_terms = int_terms0, zero_terms = Zt, int_xlev = int_xlev, zero_xlev = zero_xlev,
                 int_contrasts = int_contrasts, zero_contrasts = zero_contrasts,
        formula = formula, zero_formula = zero, call = match.call()), class = "zi_cpb"),
      ord, c("y", "weights", "pi_full", "lambda_full", "p0_full", "zero_group", "ceiling", "Z")))
  }

  ## per-observation log P(Y = y) and P(Y = 0) for the untruncated CPB (C++)
  ms <- as.integer(max.support)
  cpb_lp <- function(beta, alpha) {
    ## `off` resolves from the enclosing frame at call time (defined below,
    ## before any call site)
    m <- cpb_lp0_cpp(c(beta, stats::qlogis(alpha)), Xs, y, off, ms, FALSE, fl)   # beta in the scaled basis
    list(lp = pmax(m[, 1], -700), p0 = m[, 2])
  }

  ## Robust estimation: good starting values from separate fits (a zero-truncated
  ## CPB on the positive counts for (beta, alpha), and a logit of the zero
  ## indicator for gamma), then DIRECT joint maximum likelihood over all three,
  ## from several dispersion starts. Starting the count model at its positives-only
  ## fit keeps lambda sized to the positives and the structural-zero probability
  ## away from the degenerate pi -> 0 optimum that coordinate-wise EM falls into
  ## under heavy zero-inflation. (This mirrors how pscl::zeroinfl seeds ZIP/ZINB.)
  pb <- ncol(X); pg <- ncol(Z)
  ## exposure offset (log scale), column-name form only by the time we get here
  off <- if (is.null(offset)) rep(0, nrow(X)) else as.numeric(data[[offset]])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  ## the optimizer works on covariates scaled to unit (weighted) standard
  ## deviation in both equations and maps the coefficients back at the end
  sx <- .ud_colscale(X, w); sz <- .ud_colscale(Z, w); Xs <- sweep(X, 2, sx, "/"); Zs <- sweep(Z, 2, sz, "/")
  pf <- if (is.null(restart))
          tryCatch(.ud_quiet_guard(cpb(formula, data[y > 0, , drop = FALSE], truncated = TRUE, se = "none",
                                       max.support = max.support, support = support,
                                       offset = if (is.null(offset)) NULL else offset, weights = weights)),
                   error = function(e) NULL)
  b0 <- if (!is.null(pf)) as.numeric(pf$coefficients)
        else { v <- stats::glm.fit(X, y, offset = off, family = stats::poisson())$coefficients
               v[!is.finite(v)] <- 0; v }
  b0 <- b0 * sx                                                          # in the scaled basis
  a0 <- if (!is.null(pf)) pf$alpha else 0.5
  g0 <- suppressWarnings(stats::glm.fit(Zs, as.numeric(is0), family = stats::binomial())$coefficients)
  g0[!is.finite(g0)] <- 0
  nll <- function(par) {
    m   <- cpb_lp0_cpp(par[1:(pb + 1)], Xs, y, off, ms, FALSE, fl)
    pit <- as.numeric(stats::plogis(Zs %*% par[(pb + 2):(pb + 1 + pg)]))
    -sum(w * ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
  }
  run_em <- function() {
    ## Expectation-maximization from the same seeds. The default "ml" branch also
    ## calls this as a cross-check: both climb the same observed-data likelihood
    ## by different routes, and on some data each escapes a basin that strands
    ## the other. Frequency weights multiply the responsibilities in both M-steps.
    b <- b0; a <- a0; g <- g0; ll_old <- -Inf; bestem <- NULL
    for (it in seq_len(maxit)) {
      pit <- as.numeric(stats::plogis(Zs %*% g)); cl <- cpb_lp(b, a)
      wr <- ifelse(is0, pit / (pit + (1 - pit) * cl$p0), 0)               # E-step: responsibilities
      g <- suppressWarnings(stats::glm.fit(Zs, wr, weights = w, family = stats::quasibinomial()))$coefficients
      g[!is.finite(g)] <- 0                                              # M-step: zero-inflation
      op <- stats::optim(c(b, stats::qlogis(a)),                         # M-step: weighted CPB
                         function(par) cpb_wnll_cpp(par, Xs, y, (1 - wr) * w, off, ms, FALSE, fl),
                         method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-7))
      b <- op$par[-length(op$par)]; a <- stats::plogis(op$par[length(op$par)])
      pit <- as.numeric(stats::plogis(Zs %*% g)); cl <- cpb_lp(b, a)
      ## floored as in the FE branch: an unfloored zero-branch log() can hit -Inf
      ## and turn the convergence test into if (NaN)
      ll <- sum(w * ifelse(is0, log(pmax(pit + (1 - pit) * cl$p0, 1e-12)),
                           log(pmax(1 - pit, 1e-12)) + cl$lp))
      if (is.finite(ll) && (is.null(bestem) || ll > bestem$ll)) bestem <- list(ll = ll, b = b, a = a, g = g, it = it)
      if (is.finite(ll) && is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) {
        if (!is.null(bestem)) bestem$done <- TRUE
        break
      }
      ll_old <- ll
    }
    bestem
  }
  if (!is.null(restart)) {
    ## the final polish only, from the restart point, on either method
    st <- c(restart$b * sx, stats::qlogis(min(max(restart$alpha, 1e-6), 1 - 1e-6)), restart$g * sz)
    v0 <- nll(st); best <- list(par = st, value = if (is.finite(v0)) v0 else Inf, counts = c(0L, 0L), convergence = 0L)
    for (rs in 1:4) {
      op <- tryCatch(stats::optim(best$par, nll, method = "Nelder-Mead", control = list(maxit = 20L * maxit, reltol = tol)),
                     error = function(e) NULL)
      if (is.null(op) || !is.finite(op$value) || op$value > best$value - 1e-4) {
        if (!is.null(op) && is.finite(op$value) && op$value < best$value) best <- op
        break
      }
      best <- op
    }
    best <- .ud_bfgs_pass(best, nll)
    b <- best$par[1:pb]; a <- stats::plogis(best$par[pb + 1]); g <- best$par[(pb + 2):(pb + 1 + pg)]
    ll <- -best$value; it <- restart$iterations; conv <- isTRUE(best$convergence == 0)
  } else if (method == "em") {
    bestem <- run_em()
    if (is.null(bestem)) stop("zi_cpb (EM): the log-likelihood was non-finite at every iteration; ",
                              "the count component may not cover the data (try method = \"ml\").")
    b <- bestem$b; a <- bestem$a; g <- bestem$g; ll <- bestem$ll; it <- bestem$it; conv <- fl || isTRUE(bestem$done)
    ll_obs <- -nll(c(b, stats::qlogis(min(max(a, 1e-12), 1 - 1e-12)), g))   # unfloored, as the ML path reports
    if (is.finite(ll_obs)) ll <- ll_obs
  } else {
    best <- NULL                                                          # direct joint ML (default)
    ## Starts that span the structural share of zeros, which is where the mixture
    ## likelihood's optima differ: the separate-stage seeds (nearly every zero
    ## structural) at three dispersions, the zero-inflated Poisson's coefficients
    ## at dispersions between the smallest alpha its rates allow and 1
    ## (.zi_zip_seeds), a low-inflation start, the single-equation CPB with a
    ## constant structural-zero probability of 0.02 (.zi_low_start), and on King's
    ## support the ZIP's direct-ML solution near alpha = 1 (.zi_zip_ml).
    zs <- tryCatch(.zi_zip_seeds(Xs, Zs, y, off, w, support = support), error = function(e) NULL)
    lw <- .zi_low_start(X, Z, y, off, w, max.support, support)
    ## on King's support the ZIP is the limit alpha -> 1: its direct-ML solution enters near that limit
    zm <- if (fl) NULL else tryCatch(.zi_zip_ml(Xs, Zs, y, off, w, zs), error = function(e) NULL)
    starts <- c(lapply(unique(c(a0, 0.3, 0.7)), function(as0) c(b0, stats::qlogis(as0), g0)),
                if (!is.null(zs)) lapply(zs$alpha, function(as0) c(zs$b, stats::qlogis(as0), zs$g)),
                if (!is.null(lw)) lapply(unique(c(lw$alpha, 0.9)), function(as0) c(lw$b * sx, stats::qlogis(as0), lw$g * sz)),
                if (!is.null(zm)) lapply(c(0.999, 1 - 1e-6), function(as0) c(zm$b, stats::qlogis(as0), zm$g)))
    for (st in starts) {
      op <- tryCatch(stats::optim(st, nll, method = "Nelder-Mead",
                                  control = list(maxit = 20L * maxit, reltol = tol)),
                     error = function(e) NULL)
      ## Restart each start from its own solution until the polish stops paying:
      ## a single Nelder-Mead run can stall on a collapsed simplex far from the
      ## optimum, and the stalled value can then win the multi-start by accident.
      for (rs in 1:3) {
        if (is.null(op)) break
        op2 <- tryCatch(stats::optim(op$par, nll, method = "Nelder-Mead",
                                     control = list(maxit = 20L * maxit, reltol = tol)),
                        error = function(e) NULL)
        if (is.null(op2) || op2$value > op$value - 1e-4) { if (!is.null(op2) && op2$value < op$value) op <- op2; break }
        op <- op2
      }
      if (!is.null(op) && is.finite(op$value) && (is.null(best) || op$value < best$value)) best <- op
    }
    if (is.null(best)) stop("zi_cpb: all optimization starts failed.")
    ## EM cross-check: the EM climbs the same observed-data likelihood by a
    ## different route and cheaply escapes the degenerate pi -> 0 basin that can
    ## strand every Nelder-Mead start. Evaluate its solution under the ML
    ## objective, polish, and keep the better optimum.
    bestem <- tryCatch(run_em(), error = function(e) NULL)
    if (!is.null(bestem)) {
      st <- c(bestem$b, stats::qlogis(min(max(bestem$a, 1e-6), 1 - 1e-6)), bestem$g)
      v0 <- tryCatch(nll(st), error = function(e) Inf)
      ope <- tryCatch(stats::optim(st, nll, method = "Nelder-Mead",
                                   control = list(maxit = 20L * maxit, reltol = tol)),
                      error = function(e) NULL)
      cand <- if (!is.null(ope) && ope$value <= v0) ope
              else if (is.finite(v0)) list(par = st, value = v0, counts = c(0L, 0L))
              else NULL
      if (!is.null(cand) && cand$value < best$value) best <- cand
    }
    best <- .ud_bfgs_pass(best, nll)
    b <- best$par[1:pb]; a <- stats::plogis(best$par[pb + 1]); g <- best$par[(pb + 2):(pb + 1 + pg)]
    ll <- -best$value; it <- as.integer(best$counts[1]); conv <- isTRUE(best$convergence == 0)
  }
  ## degeneracy guard: the untruncated CPB cannot always satisfy its feasibility
  ## constraint under heavy zero-inflation with a wide count range
  if (mean(cpb_lp(b, a)$lp[!is0] <= -699) > 0.01)
    warning("zi_cpb: the count component cannot satisfy the CPB feasibility constraint for a ",
            "nontrivial share of observations; the mixture fit may be unreliable and a hurdle ",
            "model is likely more appropriate.")

  p0_full <- cpb_lp(b, a)$p0                                            # the count component's P(0) at the fit
  b <- b / sx; g <- g / sz                                              # back to the covariates' units
  names(b) <- colnames(X); names(g) <- colnames(Z)
  lambda_full <- as.numeric(exp(off + X %*% b))
  structure(list(coefficients = b, alpha = a, zero_coef = g, zero.coefficients = g, loglik = ll,
                 iterations = it, offset = off, df = length(b) + 1L + length(g), n = n, y = y, converged = conv,
                 weights = if (is.null(weights)) NULL else w, nobs_weighted = if (is.null(weights)) n else sum(w),
                 pi_full = as.numeric(stats::plogis(Z %*% g)), lambda_full = lambda_full, p0_full = p0_full,
                 zero_group = .ud_fe_group(data, zero_fe),
                 ceiling = lambda_full / (1 - a), X = X, Z = Z, fe = NULL, support = support,
                 int_beta = b, int_xref = colMeans(X), int_fe_ref = mean(off), max.support = ms,
                 zero_beta = g, zero_xref = .fe_dummy_ref(colMeans(Z), Z, y, zero_fe, .ident_binary),
                 int_terms = int_terms0, zero_terms = Zt, int_xlev = int_xlev, zero_xlev = zero_xlev,
                 int_contrasts = int_contrasts, zero_contrasts = zero_contrasts,
                 formula = formula, zero_formula = zero,
                 call = match.call()),
            class = "zi_cpb")
}

#' @method print zi_cpb
#' @export
print.zi_cpb <- function(x, ...) {
  cat("Zero-Inflated Continuous Parameter Binomial")
  if (!is.null(x$fe)) cat(" (intensity fixed effects on '", x$fe, "')", sep = "")
  cat("\nCall:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("N: %d   %s: %d   logLik: %.1f\n", x$n,
              if (identical(x$refit_args$method, "em")) "EM iterations" else "optimizer iterations",
              x$iterations, x$loglik))
  cat("\nIntensity (CPB) coefficients:\n"); print(round(x$coefficients, 4))
  cat("Intensity alpha (shape):", round(x$alpha, 4), "\n")
  .cpb_lower_note(x$alpha)
  cat("\nZero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n")
  print(round(x$zero_coef, 4))
  cat("Mean structural-zero probability:", round(mean(x$pi_full, na.rm = TRUE), 3), "\n")
  if (isTRUE(x$n_zero_units > 0L)) cat(.zi_zero_note(x), "\n")
  .ud_zero_df_note(x$zero_fe_zero_units, "zero-equation dummy", "zero-equation dummies")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  .zi_separation_note(x)
  invisible(x)
}
.zi_zero_note <- function(x)
  if (isTRUE(x$n_zero_units > 0L))
    paste0("Note: ", .zi_zero_units_text(x), " left out of the estimation (their count effect at its limit, ",
           "their zeros certain).")

#' Predict from a zi_cpb fit
#'
#' @param object A `"zi_cpb"` object.
#' @param newdata Data frame of covariate profiles. A `zero_fe` unit the fit has
#'   not seen, or every row when that column is absent, enters at the zero
#'   equation's reference unit, and a unit at a limit at that limit (see
#'   `zero_fe` in [zi_cpb()]); for the count equation's `fe`, see `offset`.
#' @param type `"response"` for the marginal expected count E(Y), `"zero"` for
#'   the structural-zero probability, `"intensity"` for the mean of the CPB
#'   component (on the floor support, its rate lambda), or `"rate"` for the
#'   component's rate lambda.
#' @param offset Optional offset (log scale) of the count component for the
#'   `newdata` rows; see [predict.cpb()]. With `fe`, a row whose `fe` column
#'   names a unit of the fit enters at that unit's count effect, a row of an
#'   unseen unit (or every row, when `newdata` has no `fe` column) at the
#'   reference unit, the mean effect of the units the data identify. A unit
#'   left out for having no positive count (see [zi_cpb()], `fe`) predicts a
#'   mean of zero, and a structural-zero probability of `NA` where its own
#'   `zero_fe` dummy would be needed.
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict zi_cpb
#' @export
predict.zi_cpb <- function(object, newdata = NULL, type = c("response", "zero", "intensity", "rate"), offset = NULL,
                           se.fit = FALSE, interval = c("none", "confidence"),
                           level = 0.95, ...) {
  type <- match.arg(type)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  bad <- FALSE
  if (is.null(newdata)) {
    pit <- object$pi_full; lam <- object$lambda_full
  } else {
    miss <- setdiff(unique(c(all.vars(object$int_terms), all.vars(object$zero_terms))),
                    c(names(newdata), .ud_zero_fe(object)))
    if (length(miss))
      stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    Xn <- .ud_newdata_matrix(object$int_terms, newdata, object$int_xlev, object$int_contrasts, names(object$coefficients))
    zn <- .zi_cpb_newdata_zero(object, newdata)
    off <- if (is.null(object$fe)) .ud_newdata_offset(object, newdata, offset, nrow(Xn))
           else .fe_newdata_effect(newdata, object$fe, object$fe_hat, object$int_fe_ref)
    bad <- .ud_na_rows(Xn, off) | is.na(zn$eta)                # a missing covariate or offset predicts NA
    if (any(bad)) { Xn[bad, ] <- 0; off[bad] <- 0 }
    lam <- as.numeric(exp(off + Xn[, names(object$coefficients), drop = FALSE] %*% object$coefficients))
    pit <- .zero_at_limit(as.numeric(stats::plogis(zn$eta)), object, newdata); pit[zn$lack] <- NA
  }
  at0 <- .count_at_limit(object, newdata)                        # a unit left out for having no positive count
  out <- switch(type, zero = pit, rate = replace(lam, at0, 0), {
    m <- replace(.zi_cpb_cmean(object, lam), at0, 0)
    if (type == "intensity") m else (1 - .zi_pi0(pit)) * m })
  .ud_mask(out, bad)
}

## The zero equation's linear predictor at `newdata` rows. A row of a unit the fit left out for having no positive
## count (see zi_cpb(), `fe`) whose own zero_fe dummy the fit lacks is built at a level the fit has and flagged in
## `lack`, so that its structural-zero probability is reported as NA. Any other zero_fe unit the fit has not seen
## (or every row, when `newdata` lacks the zero_fe column) enters at the zero equation's reference unit
## (`zero_xref`), and with `all_factors` (held-out scoring on King's support) an unseen level of any factor.
.zi_cpb_newdata_zero <- function(object, newdata, all_factors = FALSE) {
  zfe <- .ud_zero_fe(object); lack <- rep(FALSE, nrow(newdata))
  if (length(object$zero_units) && length(zfe) == 1L && zfe %in% names(newdata) && object$fe %in% names(newdata)) {
    have <- unique(as.character(object$zero_group[object$fit_rows]))
    lack <- as.character(newdata[[object$fe]]) %in% object$zero_units & !(as.character(newdata[[zfe]]) %in% have)
    if (any(lack)) { v <- as.character(newdata[[zfe]]); v[lack] <- have[1L]; newdata[[zfe]] <- v }
  }
  eta <- if (!.ud_has_unseen(object$zero_xlev, newdata, zfe, all_factors))
           as.numeric(.ud_newdata_matrix(object$zero_terms, newdata, object$zero_xlev, object$zero_contrasts,
                                         names(object$zero_coef)) %*% object$zero_coef)
         else .ud_dummy_eta(object$zero_terms, object$zero_xlev, object$zero_coef, newdata, object$zero_contrasts,
                            object$zero_xref, zfe, all_factors)
  list(eta = eta, lack = lack)
}

## the count component's mean at rate lam: the exact mean of King's CPB; on the
## floor support the rate itself, as in versions before 0.2.0 (the two differ
## where lambda/(1-alpha) is not an integer)
.zi_cpb_cmean <- function(object, lam) {
  sup <- .cpb_support(object)
  if (.cpb_fl(sup)) lam else .cpb_mean(lam, object$alpha, FALSE, sup)
}
