## ---------------------------------------------------------------------------
## King's generalized event count (Katz-family) regression: one dispersion
## parameter delta (Var/Mean) spanning under-, equi-, and over-dispersion.
## ---------------------------------------------------------------------------

## The GEC support guard. The recursion refuses a rate whose terms have not
## fallen below relative precision by `max.support` (the pmf is never
## renormalized over a cut support): a fit reports whether its largest fitted
## support reaches the guard, and a prediction beyond it is NA with a warning.
.gec_guard_binding <- function(mu, delta, max.support) {
  kl <- gec_klen_cpp(mu, delta, as.integer(max.support))
  b <- anyNA(kl) || max(kl) >= max.support
  if (b) warning("the fitted support reaches max.support (", max.support, "): the dispersion is bounded by the guard, ",
                 "not the data. Raise max.support.", call. = FALSE)
  b
}
.gec_mean_guarded <- function(mu, delta, max.support, truncated, p0fun = NULL) {
  m <- gec_mean_cpp(mu, delta, as.integer(max.support))
  if (anyNA(m)) warning(sum(is.na(m)), " prediction(s) lie beyond max.support (", max.support,
                        ") and are NA; raise max.support.", call. = FALSE)
  if (truncated) { p0 <- p0fun(); m <- m / pmax(1 - p0, 1e-8) }
  m
}

#' Generalized event count (Katz-family) regression
#'
#' Fits King's generalized event count model, a Katz-family count regression
#' whose single dispersion parameter `delta` (the variance-to-mean ratio) is
#' estimated freely and spans underdispersion (`delta < 1`, a finite-support
#' member; the continuous parameter binomial is this cell), equidispersion
#' (`delta = 1`, Poisson), and overdispersion (`delta > 1`, the negative
#' binomial). Unlike [cpb()], which fixes the direction of dispersion to
#' under, `gec()` lets the data choose. The Katz recursion delivers exact first
#' and second moments---the Winkelmann--Signorino--King correction realized
#' directly on an unbounded support---and the likelihood is evaluated in C++.
#' For `delta < 1` the support is finite and the renormalized distribution's
#' mean and variance equal `exp(x'b)` and `delta * exp(x'b)` only when
#' `exp(x'b)/(1 - delta)` is an integer, so `exp(x'b)` is the rate parameter of
#' the recursion, the fitted mean is computed exactly from the pmf, and `delta`
#' is the Katz dispersion parameter (the variance-to-mean ratio on an unbounded
#' support). The optimizer works on covariates scaled to unit standard
#' deviation and maps the coefficients back, so the fit does not depend on the
#' covariates' units.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param truncated Logical; if `TRUE`, fit the zero-truncated GEC (all `Y >= 1`),
#'   the intensity model of [hurdle_gec()].
#' @param se Coefficient inference: `"none"` (default; fast) or `"bootstrap"`.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector) for a cluster/block bootstrap; see [cpb()].
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [cpb()].
#' @param cores Worker processes for the bootstrap; see [cpb()].
#' @param max.support Guard on the maximum evaluated support.
#' @param maxit,reltol Optimizer controls.
#' @return An object of class `"gec"` with `coefficients`, `delta` (the estimated
#'   dispersion), `loglik`, bootstrap standard errors, and bookkeeping.
#' @examples
#' set.seed(1); x <- rnorm(400)
#' y <- rpois(400, exp(1 + 0.5 * x))
#' gec(y ~ x, data = data.frame(y = y, x = x), se = "none")
#' @seealso [cpb()], [compare_dispersion()], [dispersion_test()]
#' @export
gec <- function(formula, data, truncated = FALSE, se = c("none", "bootstrap"), B = 500, cluster = NULL,
                offset = NULL, weights = NULL, cores = 1L, max.support = 500, maxit = 20000, reltol = 1e-8) {
  .ud_no_formula_offset(formula)
  se <- match.arg(se); cl <- match.call()
  .ud_check_whole(max.support, "max.support")
  if (se == "bootstrap") .ud_check_whole(B, "B", lower = 2)
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit, drop.unused.levels = TRUE)
  Y  <- stats::model.response(mf); X <- stats::model.matrix(formula, mf)
  n  <- length(Y); p <- ncol(X); rows <- .ud_kept_rows(mf, data)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (n == 0L) stop("No complete cases on the model variables.")
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  .ud_warn_all_zero(Y)
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")
  Y <- as.integer(Y)
  .ud_rank_check(X)
  off <- rep_len(0, n)                                          # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    if (is.null(ov)) stop("'offset' must name a column of 'data' or be a numeric vector.")
    if (length(ov) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    off <- as.numeric(ov[rows])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  w <- .ud_weights(weights, data, rows); wv <- .ud_w1(w, n)
  clid <- NULL
  if (!is.null(cluster)) {
    if (se != "bootstrap") {
      warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
    } else {
      cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
      if (is.null(cv)) stop("'cluster' must be a column name in 'data' or a vector.")
      if (length(cv) != nrow(data)) stop("'cluster' must have one value per row of 'data'.")
      clid <- cv[rows]
      if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
    }
  }
  ms <- as.integer(max.support)
  objf <- if (!truncated) function(par, Xm, Ym, om, wm) gec_wnll_cpp(par, Xm, Ym, wm, om, ms)
          else function(par, Xm, Ym, om, wm) {       # zero-truncated Katz likelihood (P(Y>=1))
            m <- gec_lp0_cpp(par, Xm, Ym, om, ms)
            if (any(m[, 2] >= 1 - 1e-12) || any(m[, 1] <= -1e299)) return(1e10)   # infeasible point
            v <- -sum(wm * (m[, 1] - log1p(-m[, 2]))); if (is.finite(v)) v else 1e10
          }
  fitone <- function(Xm, Ym, om, wm) {
    s <- .ud_colscale(Xm, wm); Xs <- sweep(Xm, 2, s, "/")
    b0 <- tryCatch({ v <- stats::glm.fit(Xs, Ym, weights = wm, offset = om, family = stats::poisson())$coefficients
                     v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
    best <- NULL
    for (ld in c(-0.5, 0, 0.5)) {                     # under / Poisson / over starts
      op <- tryCatch(stats::optim(c(b0, ld), function(par) objf(par, Xs, Ym, om, wm),
                                  method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)),
                     error = function(e) NULL)
      if (!is.null(op) && op$value < 1e9 && (is.null(best) || op$value < best$value)) best <- op
    }
    if (!is.null(best)) {
      best <- .ud_nm_polish(best, function(par) objf(par, Xs, Ym, om, wm), maxit, reltol)
      best$par[seq_len(p)] <- best$par[seq_len(p)] / s
    }
    best
  }
  fit <- fitone(X, Y, off, wv)
  if (is.null(fit) || fit$value >= 1e9)
    stop("no feasible fit: at every start some count lay outside the model's support or the recursion ",
         "reached max.support (currently ", ms, "); raise max.support or check the data.")
  beta   <- setNames(fit$par[1:p], colnames(X)); delta <- exp(fit$par[p + 1]); loglik <- -fit$value
  mu     <- as.numeric(exp(off + X %*% beta))
  fitted <- if (!truncated) gec_mean_cpp(mu, delta, ms)
            else { p0 <- gec_lp0_cpp(c(beta, log(delta)), X, Y, off, ms)[, 2]; gec_mean_cpp(mu, delta, ms) / pmax(1 - p0, 1e-8) }
  support_binding <- .gec_guard_binding(mu, delta, ms)

  se.beta <- setNames(rep(NA_real_, p), colnames(X)); se.delta <- NA_real_; boot <- NULL; ci.beta <- NULL
  if (se == "bootstrap") {
    if (!is.null(clid)) .ud_cluster_guard(clid)
    grp <- if (!is.null(clid)) split(seq_len(n), clid) else NULL
    one <- function(b) {
      idx <- if (is.null(grp)) sample.int(n, n, replace = TRUE)
             else unlist(grp[sample.int(length(grp), length(grp), replace = TRUE)], use.names = FALSE)
      fb <- fitone(X[idx, , drop = FALSE], Y[idx], off[idx], wv[idx])
      if (!is.null(fb)) c(fb$par[1:p], exp(fb$par[p + 1])) else rep(NA_real_, p + 1L)
    }
    boot <- do.call(rbind, .ud_lapply(seq_len(B), one, cores))
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta  <- setNames(apply(ok[, 1:p, drop = FALSE], 2, stats::sd), colnames(X))
      se.delta <- stats::sd(ok[, p + 1])
      ci.beta  <- t(apply(ok[, 1:p, drop = FALSE], 2, stats::quantile, c(.025, .975)))
      rownames(ci.beta) <- colnames(X)
    } else warning("Too few bootstrap resamples converged for stable inference.")
    attr(boot, "nboot_ok") <- nrow(ok)
    if (nrow(ok) < B / 2) warning(sprintf("only %d of %d bootstrap replicates converged; the standard errors rest on the survivors.", nrow(ok), B), call. = FALSE)
  }

  structure(list(coefficients = beta, delta = delta, dispersion = delta, loglik = loglik,
                 se.beta = se.beta, se.delta = se.delta, boot = boot, ci.beta = ci.beta,
                 fitted.values = fitted, linear.predictors = as.numeric(off + X %*% beta), offset = off,
                 weights = w, nobs_weighted = if (is.null(w)) n else sum(w),
                 n = n, df = p + 1, p = p, se.type = se, clustered = !is.null(clid),
                 converged = fit$convergence == 0,
                 truncated = truncated, max.support = ms, support_binding = support_binding,
                 formula = formula, terms = attr(mf, "terms"),
                 levels = .cpb_xlevels(mf), contrasts = attr(X, "contrasts"),
                 X = X, Y = Y, call = cl), class = "gec")
}

#' @method print gec
## direction label for the printed dispersion: the pooled fit's label reports the
## likelihood-ratio test when the calibration rule admits its asymptotic
## distribution (a print method runs no bootstrap); otherwise, and for the
## concentrated and two-part fits, it is the point estimate with its band stated
.gec_direction <- function(x) {
  d <- x$delta
  p <- if (inherits(x, "gec") && !inherits(x, "gec_fe"))
         tryCatch({ dt <- suppressWarnings(dispersion_test(x, B = 0))
                    if (isTRUE(dt$asymptotic_ok)) dt$p.value else NA_real_ },
                  error = function(e) NA_real_) else NA_real_
  if (is.finite(p)) {
    if (p >= 0.05) sprintf("equidispersion not rejected (LR p = %.2f)", p)
    else if (d < 1) sprintf("underdispersed (LR p = %.2g)", p) else sprintf("overdispersed (LR p = %.2g)", p)
  } else if (d < 0.97) "underdispersed (point estimate)" else if (d > 1.03) "overdispersed (point estimate)"
  else "within 3% of equidispersion (point estimate)"
}

#' @export
print.gec <- function(x, ...) {
  disp <- .gec_direction(x)
  cat("Generalized event count (Katz family) regression\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta)))
    print(round(cbind(Estimate = x$coefficients, `Std. Error` = x$se.beta), 4))
  else print(round(x$coefficients, 4))
  cat(sprintf("\ndispersion delta (Katz; Var/Mean on an unbounded support) = %.3f  [%s]", x$delta, disp))
  if (is.finite(x$se.delta)) cat(sprintf("  (SE %.3f)", x$se.delta))
  cat(sprintf("\nlogLik = %.2f,  n = %d\n", x$loglik, x$n))
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' @method coef gec
#' @export
coef.gec <- function(object, ...) object$coefficients

#' @method logLik gec
#' @export
logLik.gec <- function(object, ...) {
  val <- object$loglik; attr(val, "df") <- object$df; attr(val, "nobs") <- nobs(object)
  class(val) <- "logLik"; val
}

#' @method nobs gec
#' @export
nobs.gec <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted

#' @method fitted gec
#' @export
fitted.gec <- function(object, ...) object$fitted.values

#' @method vcov gec
#' @export
vcov.gec <- function(object, ...) {
  if (is.null(object$boot)) return(NULL)                        # no bootstrap requested
  p  <- length(object$coefficients)                            # covariate block (works for gec AND gec_fe)
  ok <- object$boot[stats::complete.cases(object$boot), 1:p, drop = FALSE]
  if (nrow(ok) < 2L) return(NULL)
  v <- stats::cov(ok); dimnames(v) <- list(names(object$coefficients), names(object$coefficients)); v
}

#' Predictions from a generalized event count fit
#'
#' @param object A `"gec"` object.
#' @param newdata Optional covariate profiles.
#' @param type `"response"` (the mean), `"link"` (the log mean), or `"prob"`
#'   (the probability of the count `at`).
#' @param at Count value(s) for `type = "prob"`.
#' @param offset Optional offset (log scale) for the count component when
#'   predicting on `newdata`; a numeric vector or a column name in `newdata`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict gec
#' @export
predict.gec <- function(object, newdata = NULL, type = c("response", "link", "prob"), at = NULL,
                        offset = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    X <- object$X; off <- if (!is.null(object$offset)) object$offset else rep_len(0, nrow(X))
  } else {
    Terms <- stats::delete.response(object$terms)
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    off <- if (is.null(offset)) rep_len(0, nrow(X))
           else as.numeric(if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset)
  }
  bad <- if (is.null(newdata)) FALSE else .ud_na_rows(X, off)     # a missing covariate or offset predicts NA
  if (any(bad)) { X[bad, ] <- 0; off[bad] <- 0 }
  mu <- as.numeric(exp(off + X %*% object$coefficients))
  out <- switch(type,
    link = log(mu),
    response = .gec_mean_guarded(mu, object$delta, object$max.support, isTRUE(object$truncated),
                                 function() gec_lp0_cpp(c(object$coefficients, log(object$delta)), X, rep(0L, length(mu)), off, object$max.support)[, 2]),
    prob = {
      if (is.null(at)) stop("For type = \"prob\", supply 'at' (the count value).")
      yv <- if (length(at) == 1) rep(as.integer(at), length(mu)) else as.integer(at)
      if (length(yv) != length(mu)) stop("'at' must be length 1 or nrow(newdata).")
      m <- gec_lp0_cpp(c(object$coefficients, log(object$delta)), X, yv, off, object$max.support)
      pk <- exp(pmax(m[, 1], -700))
      if (isTRUE(object$truncated)) pk <- ifelse(yv == 0L, 0, pk / (1 - m[, 2]))  # condition on Y > 0
      pk
    })
  .ud_mask(out, bad)
}

## ---------------------------------------------------------------------------
## GEC with high-dimensional unit fixed effects (concentrated likelihood)
## ---------------------------------------------------------------------------

#' GEC (Katz-family) regression with high-dimensional unit fixed effects
#'
#' Fits [gec()] with a full set of unit fixed effects, concentrating (profiling)
#' out the unit intercepts by a one-dimensional inner maximization per unit, so
#' the outer optimizer handles only the covariate coefficients and the free
#' dispersion parameter `delta`. This is the [cpb_fe()] concentration generalized
#' to the whole Katz family: the panel can be under-, equi-, or overdispersed and
#' the direction is estimated, not presumed.
#' For `delta < 1` the Katz support is finite, but the likelihood is continuous across an
#' integer ceiling (the entering support point's mass grows from zero) and only kinks
#' there, so each unit's objective is unimodal in its intercept and is maximized by golden
#' section on a bracket around the unit's Poisson intercept, expanded while the maximum
#' sits at an edge; a maximum at a kink is tracked in the analytic gradient of the
#' concentrated likelihood, which the outer BFGS uses. The feasibility floor is exact:
#' every count must lie in `0..ceiling(mu/(1-delta))`.
#' Runtime: about twice that of [cpb_fe()] on the same panel (the Katz recursion
#' runs the whole support), so a 2,600-row panel in 146 units takes about two minutes.
#'
#' Limitation: unlike [cpb_fe()], `gec_fe()` has no `truncated` argument -- a
#' concentrated zero-truncated GEC is not currently implemented. For a
#' zero-truncated GEC with fixed effects, enter the unit factor as dummies in
#' [gec()]'s formula (feasible for moderate unit counts).
#'
#' Short panels: `delta` and the fixed effects carry the incidental-parameters
#' bias of nonlinear fixed-effects estimation, of order 1/T; treat `delta`
#' cautiously below about T = 30. The covariate coefficients are not materially
#' affected. `bias_correct = "jackknife"` removes the leading 1/T term by the
#' split-panel jackknife exactly as in [cpb_fe()]: refit on each unit's temporal
#' halves and report `2 * full - mean(halves)`, with the unit effects and fitted
#' values re-concentrated at the corrected parameters and `logLik`/`AIC` kept at
#' the maximum-likelihood fit (uncorrected estimates in `$uncorrected`). The
#' same time-homogeneity validity gate applies: on panels where the two halves
#' do not estimate a common parameter the correction is REFUSED with a warning
#' and the maximum-likelihood fit is returned (see [cpb_fe()], Details).
#'
#' @param formula A model formula for the covariates only (no unit factor, no
#'   intercept; the fixed effects absorb it).
#' @param data A data frame.
#' @param fe Name of the column holding the unit identifier.
#' @param se Inference for the covariate coefficients: `"none"` (default) or
#'   `"bootstrap"`, a pairs/cluster bootstrap over units.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Cluster for the bootstrap; `NULL` (default) resamples the
#'   fixed-effects units. See [cpb_fe()].
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [cpb()].
#' @param cores Worker processes for the bootstrap; see [cpb()].
#' @param max.support Guard on the maximum evaluated support.
#' @param inner_it Golden-section iterations for each unit's inner maximization.
#' @param maxit,reltol Outer optimizer controls.
#' @param bias_correct `"none"` (default) or `"jackknife"`, the split-panel
#'   jackknife correction for the 1/T incidental-parameters bias (see Details).
#' @return An object of class `c("gec_fe", "gec")`.
#' @examples
#' \donttest{
#' set.seed(5)
#' d <- do.call(rbind, lapply(1:30, function(i) {
#'   x <- rnorm(10); N <- pmax(round(exp(1 + rnorm(1, 0, 0.4) + 0.3 * x) / 0.5), 1)
#'   data.frame(unit = i, x = x, y = rbinom(10, N, 0.5))
#' }))
#' gec_fe(y ~ x, data = d, fe = "unit")   # delta ~ 0.5: underdispersed
#' }
#' @seealso [gec()], [cpb_fe()]
#' @export
gec_fe <- function(formula, data, fe, se = c("none", "bootstrap"), B = 500, cluster = NULL,
                   offset = NULL, weights = NULL, cores = 1L, max.support = NULL, inner_it = 30L,
                   maxit = 3000L, reltol = 1e-7, bias_correct = c("none", "jackknife")) {
  .ud_no_formula_offset(formula)
  data <- .ud_drop_na_fe(data, fe)
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  se <- match.arg(se); bias_correct <- match.arg(bias_correct)
  .ud_check_whole(max.support, "max.support", null_ok = TRUE); .ud_check_whole(inner_it, "inner_it")
  if (se == "bootstrap") .ud_check_whole(B, "B", lower = 2)
  if (!is.character(fe) || length(fe) != 1L || !fe %in% names(data)) stop("'fe' must name a column of 'data'.")
  if (!is.null(cluster) && se != "bootstrap")
    warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
  if (!is.null(offset) && !(is.character(offset) && length(offset) == 1L)) {
    if (length(offset) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    data[[".gec_off"]] <- as.numeric(offset); offset <- ".gec_off"   # vector -> column, sorts with data
  }
  if (!is.null(weights) && !(is.character(weights) && length(weights) == 1L)) {
    if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
    data[[".gec_w"]] <- as.numeric(weights); weights <- ".gec_w"
  }
  data <- data[order(data[[fe]]), , drop = FALSE]
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit, drop.unused.levels = TRUE)
  rows <- .ud_kept_rows(mf, data)
  Y  <- stats::model.response(mf)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  .ud_warn_all_zero(Y)
  Y <- as.integer(Y)
  off <- if (is.null(offset)) rep_len(0, length(Y)) else as.numeric(data[[offset]][rows])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  w <- .ud_weights(weights, data, rows); wv <- .ud_w1(w, length(Y))
  Xf <- stats::model.matrix(formula, mf)
  keep <- setdiff(colnames(Xf), "(Intercept)")
  if (!length(keep))
    stop("Provide at least one covariate in 'formula'; the fixed effects are the intercepts.")
  X  <- Xf[, keep, drop = FALSE]
  uf <- factor(data[[fe]][rows]); nu <- nlevels(uf)
  .ud_rank_check_fe(X, as.integer(uf))
  if (is.null(max.support)) max.support <- max(500L, 10L * max(Y))
  ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
  p  <- ncol(X); ms <- as.integer(max.support); it <- as.integer(inner_it)
  s  <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/")
  bs <- tryCatch({ v <- stats::glm.fit(cbind(1, Xs), Y, weights = wv, offset = off, family = stats::poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  ## concentrated objective, every unit searched from cold (src/gec_fe.cpp);
  ## .fe_outer() runs the starts, the polish, and the cold re-evaluation
  nll_raw <- function(par, awarm) gec_fe_nll_cpp(par, Xs, Y, off, wv, ustart, nu, ms, it, awarm)
  grad_raw <- function(par, a, edge) gec_fe_grad_cpp(par, Xs, Y, off, wv, ustart, nu, ms, a, edge)
  best <- .fe_outer(nll_raw, grad_raw, bs, c(-0.5, 0.2), maxit, reltol)
  if (is.null(best)) stop("All optimization starts failed.")
  if (best$value >= 1e9)
    stop("no feasible fit: at every start some count lay outside the model's support or the recursion ",
         "reached max.support (currently ", ms, "); raise max.support or check the data.")
  beta <- setNames(best$par[1:p] / s, keep); delta <- exp(best$par[p + 1])
  fe_hat <- best$a
  rate <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
  fitted <- gec_mean_cpp(rate, delta, ms)
  support_binding <- .gec_guard_binding(rate, delta, ms)

  se.beta <- setNames(rep(NA_real_, p), keep); ci.beta <- NULL; boot <- NULL; clab <- NULL
  if (se == "bootstrap") {
    fev  <- as.character(data[[fe]][rows]); dest <- data[rows, , drop = FALSE]
    cval <- if (is.null(cluster)) fev else .ud_cluster_values(cluster, data, rows)
    clab <- if (is.null(cluster)) fe else if (is.character(cluster) && length(cluster) == 1L) cluster else "custom"
    .ud_cluster_guard(cval)
    grp  <- split(seq_along(Y), cval)
    one <- function(b) {
      gs   <- sample.int(length(grp), length(grp), replace = TRUE)
      rws  <- unlist(grp[gs], use.names = FALSE)
      bunit <- unlist(lapply(seq_along(gs), function(k) paste0(k, "_", fev[grp[[gs[k]]]])), use.names = FALSE)
      bd <- dest[rws, , drop = FALSE]; bd[[".bootunit"]] <- bunit
      fb <- tryCatch(.ud_quiet_guard(gec_fe(formula, data = bd, fe = ".bootunit", se = "none", offset = offset,
                                            weights = weights, max.support = max.support, inner_it = inner_it,
                                            maxit = maxit, reltol = reltol)),
                     error = function(e) NULL)
      if (!is.null(fb)) fb$coefficients[keep] else rep(NA_real_, p)
    }
    boot <- do.call(rbind, .ud_lapply(seq_len(B), one, cores))
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta <- setNames(apply(ok, 2, stats::sd), keep)
      z <- stats::qnorm(0.975)
      ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta); rownames(ci.beta) <- keep
    } else warning("Too few bootstrap resamples converged for stable inference.")
    attr(boot, "nboot_ok") <- nrow(ok)
    if (nrow(ok) < B / 2) warning(sprintf("only %d of %d bootstrap replicates converged; the standard errors rest on the survivors.", nrow(ok), B), call. = FALSE)
  }

  ## split-panel jackknife: refit on each unit's temporal halves and remove the
  ## leading 1/T incidental-parameters bias term (see .fe_jackknife)
  uncorrected <- NULL
  if (bias_correct == "jackknife") {
    dest <- data[rows, , drop = FALSE]
    refit <- function(dd) tryCatch({
      f <- .ud_quiet_guard(gec_fe(formula, data = dd, fe = fe, se = "none", offset = offset, weights = weights,
                                  max.support = max.support, inner_it = inner_it, maxit = maxit, reltol = reltol))
      c(f$coefficients[keep], f$delta)
    }, error = function(e) NULL)
    jk <- .fe_jackknife(dest, uf, refit, c(beta, delta),
                        resid  = (Y - fitted) / sqrt(pmax(delta * fitted, 1e-12)),
                        torder = stats::ave(seq_along(Y), as.integer(uf), FUN = seq_along),
                        disp_lower = 0.01, disp_upper = Inf)
    if (!is.null(jk$refused)) {
      warning("bias_correct = \"jackknife\" REFUSED: ", jk$refused,
              ". Returning the uncorrected maximum-likelihood estimates.")
      bias_correct <- "none"
    } else {
      for (wmsg in jk$warn) warning("bias_correct = \"jackknife\": ", wmsg, ".")
      uncorrected <- list(coefficients = beta, delta = delta)
      beta  <- setNames(jk$par[seq_len(p)], keep)
      delta <- unname(jk$par[p + 1L])
      fe_hat <- gec_fe_intercepts_cpp(c(beta, log(delta)), X, Y, off, wv, ustart, nu, ms, it)
      rate   <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
      fitted <- gec_mean_cpp(rate, delta, ms)
      support_binding <- .gec_guard_binding(rate, delta, ms)
      if (!is.null(ci.beta)) {                        # recenter the normal-approx interval
        z <- stats::qnorm(0.975)
        ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
        rownames(ci.beta) <- keep
      }
    }
  }

  structure(list(coefficients = beta, delta = delta, dispersion = delta, loglik = -best$value,
                 se.beta = se.beta, ci.beta = ci.beta, boot = boot, se.type = se, cluster = clab,
                 fe = setNames(fe_hat, levels(uf)), fitted.values = fitted, rate = rate, offset = off,
                 weights = w, nobs_weighted = if (is.null(w)) length(Y) else sum(w),
                 xref = colMeans(X), fe_ref = mean(fe_hat),
                 linear.predictors = log(rate), n = length(Y), n_units = nu, Y = Y,
                 df = p + nu + 1L, max.support = ms, support_binding = support_binding, formula = formula, fe_var = fe,
                 X = X, unit = as.integer(uf), levels = .cpb_xlevels(mf), contrasts = attr(Xf, "contrasts"),
                 converged = best$convergence == 0,
                 bias_correct = bias_correct, uncorrected = uncorrected,
                 call = match.call()), class = c("gec_fe", "gec"))
}

#' @method print gec_fe
#' @export
print.gec_fe <- function(x, ...) {
  disp <- .gec_direction(x)
  cat("GEC (Katz family) regression with", x$n_units, "unit fixed effects (concentrated likelihood)\n")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta))) {
    cat("Coefficients (", x$cluster, "-clustered bootstrap SEs):\n", sep = "")
    print(round(cbind(Estimate = x$coefficients, `Std. Error` = x$se.beta[names(x$coefficients)]), 4))
  } else { cat("Coefficients:\n"); print(round(x$coefficients, 4)) }
  cat(sprintf("\ndispersion delta (Katz; Var/Mean on an unbounded support) = %.3f  [%s],  n = %d\n", x$delta, disp, x$n))
  if (identical(x$bias_correct, "jackknife"))
    cat("Estimates are split-panel jackknife bias-corrected; logLik/AIC refer to the\nuncorrected maximum-likelihood fit. See ?gec_fe.\n")
  else
    cat("Note: delta is subject to incidental-parameters bias for short panels; see ?gec_fe.\n")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' @method predict gec_fe
#' @export
predict.gec_fe <- function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) { rate <- object$rate } else {
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    bad <- .ud_na_rows(X); if (any(bad)) X[bad, ] <- 0         # a missing covariate predicts NA
    rate <- as.numeric(exp(mean(object$fe) + X %*% object$coefficients))
  }
  if (is.null(newdata)) bad <- FALSE
  .ud_mask(switch(type, link = log(rate), response = .gec_mean_guarded(rate, object$delta, object$max.support, FALSE)), bad)
}
