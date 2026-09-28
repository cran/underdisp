## ---------------------------------------------------------------------------
## CPB regression with high-dimensional unit fixed effects (concentrated likelihood)
## ---------------------------------------------------------------------------

#' CPB regression with high-dimensional unit fixed effects
#'
#' Fits the CPB with a full set of unit fixed effects by concentrating (profiling)
#' out the unit intercepts. Each unit's intercept is solved by a one-dimensional
#' inner maximization, so the outer optimizer handles only the covariate coefficients
#' and the dispersion parameter. This avoids an explicit unit-dummy design matrix and
#' scales to thousands of units, where dummy-based fitting exhausts memory. The
#' concentrated log-likelihood equals the full-dummy log-likelihood to the
#' tolerance of the inner search, which is solved exactly. Each unit's objective
#' is a saw-tooth in its intercept: whenever the intercept crosses
#' \eqn{\log(k(1-\alpha)) - o_{it}} for an integer \eqn{k} the support of
#' observation \eqn{t} gains the point \eqn{k}, its normalizing constant grows,
#' and the objective drops by about \eqn{(1-\alpha)^k}; between such
#' breakpoints it is smooth and unimodal. Its supremum is therefore attained at
#' an interior critical point of one tooth or as the left limit at a breakpoint,
#' and the inner search evaluates every breakpoint within a window of the
#' envelope's peak (located on a coarse grid around the unit's Poisson intercept)
#' and refines inside whole teeth by golden section; every evaluation of the
#' concentrated likelihood searches each unit afresh, so the objective is a
#' function of the parameters alone. Solving the inner problem to its supremum makes the concentrated likelihood
#' continuous in `alpha` (with kinks where a unit's maximizing tooth switches) but
#' not in the slopes: where a unit's supremum sits on its feasibility floor (an
#' observation whose count equals its ceiling) the intercept cannot retreat, so
#' when another observation's breakpoint crosses that floor the concentrated
#' likelihood drops by about \eqn{(1-\alpha)^k}, a cliff in the slopes on whose
#' edge the maximizer can rest with a non-zero one-sided gradient. The outer
#' optimizer is therefore a deterministic multistart: BFGS with the analytic
#' envelope gradient from the Poisson slopes and two dispersion starts, a
#' Nelder-Mead polish, and two perturbed restarts; alternative local optima
#' differ by the order of the jumps (a few hundredths of a log-likelihood unit
#' on a 12-unit panel), far inside the bootstrap variability of the estimates.
#'
#' Short panels: the dispersion parameter `alpha` and the fixed effects are subject to
#' the incidental-parameters bias of nonlinear fixed-effects estimation. The bias in
#' `alpha` is downward and of order 1/T, where T is the per-unit number of
#' observations. It is small for T at least about 30 and should be treated cautiously
#' for short panels. The covariate coefficients are not materially affected.
#' `bias_correct = "jackknife"` removes the leading 1/T term by the split-panel
#' (half-panel) jackknife: the model is refit on the first and second temporal
#' halves of every unit's series (rows are split in the order supplied, which
#' should be temporal order) and the corrected estimate is
#' `2 * full - mean(halves)`. In the package's fixed-effects bias Monte Carlo
#' (150 units, alpha = 0.5) the alpha bias at T = 6/10/20/40 falls from
#' -0.153/-0.100/-0.057/-0.033 to -0.046/-0.024/-0.013/-0.009. With correction on,
#' `coefficients` and `alpha`
#' are the corrected estimates (the maximum-likelihood values are kept in
#' `$uncorrected`), the unit effects and fitted values are re-concentrated at the
#' corrected parameters, and `logLik`/`AIC` continue to refer to the
#' maximum-likelihood fit.
#'
#' Validity gate: the split-panel identity requires the two half-panels to
#' estimate the same pseudo-true parameter (time-homogeneity; Dhaene & Jochmans
#' 2015, Review of Economic Studies). On trending or time-heterogeneous panels
#' the correction is invalid, so the function REFUSES it -- returning the
#' uncorrected maximum-likelihood fit with a warning naming the failed check --
#' whenever the corrected dispersion leaves its feasible space or the two
#' half-panel dispersion estimates disagree beyond what the 1/T bias can
#' explain. A refusal is diagnostic information about the panel, not an error.
#' The gate is deliberately powered over sized: in the package's calibration
#' Monte Carlo it refuses about 7 percent of genuinely time-homogeneous panels
#' (a conservative nuisance; the returned fit is exactly the ordinary
#' maximum-likelihood estimate) while catching 98 percent of dispersion regime
#' changes and 100 percent of smooth unmodeled trends -- the cases where an
#' uncaught correction would be silently wrong.
#'
#' One set of fixed effects is concentrated out. For two-way (e.g. unit and time)
#' fixed effects, add `factor(time)` to `formula` (entered as dummies), or use
#' [count_reg()], which supports two-way fixed effects natively.
#'
#' @param formula A model formula for the covariates only. Do not include the unit
#'   factor; supply it via `fe`. The intercept is absorbed by the fixed effects.
#' @param data A data frame.
#' @param fe Name of the column holding the unit identifier.
#' @param truncated Logical; fit the zero-truncated CPB if `TRUE`. **Note the
#'   sibling default differs:** `cpb_fe()` defaults to `FALSE` (whole-panel
#'   fits including zeros), while [cpb()] defaults to `TRUE` (positives-only
#'   fits). State `truncated` explicitly when moving a specification between
#'   the two.
#' @param se Inference for the covariate coefficients: `"none"` (default) or
#'   `"bootstrap"`, a pairs/cluster bootstrap.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Cluster for the bootstrap. `NULL` (default) resamples the
#'   fixed-effects units themselves (unit-clustered inference, the natural panel
#'   default); a column name or vector resamples those clusters while preserving
#'   the unit fixed effects. Duplicated blocks are relabeled so the concentrated
#'   likelihood treats them as distinct units.
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [cpb()].
#' @param cores Worker processes for the bootstrap; see [cpb()].
#' @param max.support Guard on the maximum evaluated support; `NULL` (default)
#'   sets `max(500, 10 * max(y))`. See [cpb()].
#' @param inner_it Golden-section iterations refining each unit's inner
#'   maximization after the grid scans (see Details).
#' @param maxit,reltol Outer optimizer controls.
#' @param bias_correct `"none"` (default) or `"jackknife"`, the split-panel
#'   jackknife correction for the 1/T incidental-parameters bias (see Details).
#' @return An object of class `"cpb_fe"` with `coefficients`, `alpha`, `loglik`, `fe`
#'   (the estimated unit intercepts), `fitted.values`, and the implied per-observation
#'   `ceiling`.
#' @seealso [cpb()], [ud_screen()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- do.call(rbind, lapply(1:25, function(i) {
#'   x <- rnorm(12); lam <- exp(rnorm(1, 0, 0.5) + 0.5 * x)
#'   N <- pmax(round(lam / 0.5), 1)
#'   data.frame(unit = i, x = x, y = rbinom(12, N, 0.5))
#' }))
#' fit <- cpb_fe(y ~ x, data = d, fe = "unit")
#' fit
#' }
#' @export
cpb_fe <- function(formula, data, fe, truncated = FALSE, se = c("none", "bootstrap"),
                   B = 500, cluster = NULL, offset = NULL, weights = NULL, cores = 1L,
                   max.support = NULL, inner_it = 30L, maxit = 3000L, reltol = 1e-7,
                   bias_correct = c("none", "jackknife")) {
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
    data[[".cpb_off"]] <- as.numeric(offset); offset <- ".cpb_off"   # vector -> column, so it sorts with data
  }
  if (!is.null(weights) && !(is.character(weights) && length(weights) == 1L)) {
    if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
    data[[".cpb_w"]] <- as.numeric(weights); weights <- ".cpb_w"
  }
  data <- data[order(data[[fe]]), , drop = FALSE]
  mf <- model.frame(formula, data, na.action = na.omit, drop.unused.levels = TRUE)
  rows <- .ud_kept_rows(mf, data)
  Y  <- model.response(mf)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  .ud_warn_all_zero(Y)
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")
  Y <- as.integer(Y)
  off <- if (is.null(offset)) rep_len(0, length(Y)) else as.numeric(data[[offset]][rows])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  w <- .ud_weights(weights, data, rows); wv <- .ud_w1(w, length(Y))
  Xf <- model.matrix(formula, mf)
  keep <- setdiff(colnames(Xf), "(Intercept)")
  if (!length(keep))
    stop("Provide at least one covariate in 'formula'; the fixed effects are the intercepts.")
  X  <- Xf[, keep, drop = FALSE]
  uf <- factor(data[[fe]][rows]); nu <- nlevels(uf)
  .ud_rank_check_fe(X, as.integer(uf))
  if (is.null(max.support)) max.support <- max(500L, 10L * max(Y))
  ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
  p  <- ncol(X)
  s  <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/")
  bs <- tryCatch({ v <- glm.fit(cbind(1, Xs), Y, weights = wv, offset = off, family = poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  ## concentrated objective, every unit searched from cold (src/cpb_fe.cpp);
  ## .fe_outer() runs the starts, the polish, and the cold re-evaluation
  nll_raw <- function(par, awarm) cpb_fe_nll_cpp(par, Xs, Y, off, wv, ustart, nu, as.integer(max.support),
                                                 truncated, as.integer(inner_it), awarm)
  grad_raw <- function(par, a, edge) cpb_fe_grad_cpp(par, Xs, Y, off, wv, ustart, nu, as.integer(max.support),
                                                     truncated, a, edge)
  best <- .fe_outer(nll_raw, grad_raw, bs, qlogis(c(0.4, 0.6)), maxit, reltol)
  if (is.null(best)) stop("All optimization starts failed.")
  if (best$value >= 1e9)
    stop("no feasible fit: every start left some count above its implied ceiling (raise max.support, ",
         "currently ", max.support, ", or check the data).")
  beta  <- setNames(best$par[1:p] / s, keep); alpha <- plogis(best$par[p + 1])
  fe_hat <- best$a
  lam <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
  fitted <- .cpb_mean(lam, alpha, truncated)
  binding <- max(lam / (1 - alpha)) >= 0.95 * max.support
  if (binding)
    warning("the fitted ceiling reaches max.support (", max.support, "): alpha is bounded by the guard, not the data. ",
            "The panel may not be underdispersed (compare gec_fe()), or raise max.support.")

  ## pairs/cluster bootstrap for the covariate coefficients: resample whole
  ## clusters (default = the fixed-effects unit) and relabel each sampled block
  ## so the concentration treats duplicated draws as distinct units.
  se.beta <- setNames(rep(NA_real_, p), keep); ci.beta <- NULL; boot <- NULL
  clab <- NULL
  if (se == "bootstrap") {
    fev  <- as.character(data[[fe]][rows]); dest <- data[rows, , drop = FALSE]
    cval <- if (is.null(cluster)) fev else .ud_cluster_values(cluster, data, rows)
    clab <- if (is.null(cluster)) fe else if (is.character(cluster) && length(cluster) == 1L) cluster else "custom"
    .ud_cluster_guard(cval)
    grp  <- split(seq_along(Y), cval)
    one <- function(b) {
      gs   <- sample.int(length(grp), length(grp), replace = TRUE)
      rws  <- unlist(grp[gs], use.names = FALSE)
      bunit <- unlist(lapply(seq_along(gs), function(k) paste0(k, "_", fev[grp[[gs[k]]]])),
                      use.names = FALSE)
      bd <- dest[rws, , drop = FALSE]; bd[[".bootunit"]] <- bunit
      fb <- tryCatch(.ud_quiet_guard(cpb_fe(formula, data = bd, fe = ".bootunit", truncated = truncated,
                                            se = "none", offset = offset, weights = weights, max.support = max.support,
                                            inner_it = inner_it, maxit = maxit, reltol = reltol)),
                     error = function(e) NULL)
      if (!is.null(fb)) fb$coefficients[keep] else rep(NA_real_, p)
    }
    boot <- do.call(rbind, .ud_lapply(seq_len(B), one, cores))
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta <- setNames(apply(ok, 2, stats::sd), keep)
      ## normal-approximation interval centered on the point estimate: FE bootstrap
      ## replicates share the incidental-parameters bias, so a percentile interval
      ## is off-center; the bootstrap SD is the robust quantity to keep.
      z <- stats::qnorm(0.975)
      ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
      rownames(ci.beta) <- keep
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
      f <- .ud_quiet_guard(cpb_fe(formula, data = dd, fe = fe, truncated = truncated, se = "none",
                                  offset = offset, weights = weights, max.support = max.support, inner_it = inner_it,
                                  maxit = maxit, reltol = reltol))
      c(f$coefficients[keep], f$alpha)
    }, error = function(e) NULL)
    jk <- .fe_jackknife(dest, uf, refit, c(beta, alpha),
                        resid  = (Y - fitted) / sqrt(pmax(alpha * lam, 1e-12)),
                        torder = stats::ave(seq_along(Y), as.integer(uf), FUN = seq_along),
                        disp_lower = 1e-3, disp_upper = 1 - 1e-3)
    if (!is.null(jk$refused)) {
      warning("bias_correct = \"jackknife\" REFUSED: ", jk$refused,
              ". Returning the uncorrected maximum-likelihood estimates.")
      bias_correct <- "none"
    } else {
      for (wmsg in jk$warn) warning("bias_correct = \"jackknife\": ", wmsg, ".")
      uncorrected <- list(coefficients = beta, alpha = alpha)
      beta  <- setNames(jk$par[seq_len(p)], keep)
      alpha <- unname(jk$par[p + 1L])
      fe_hat <- cpb_fe_intercepts_cpp(c(beta, qlogis(alpha)), X, Y, off, wv, ustart, nu,
                                      as.integer(max.support), truncated, as.integer(inner_it))
      lam <- as.numeric(exp(off + fe_hat[as.integer(uf)] + X %*% beta))
      fitted <- .cpb_mean(lam, alpha, truncated)
      if (!is.null(ci.beta)) {                        # recenter the normal-approx interval
        z <- stats::qnorm(0.975)
        ci.beta <- cbind(lower = beta - z * se.beta, upper = beta + z * se.beta)
        rownames(ci.beta) <- keep
      }
    }
  }

  structure(list(coefficients = beta, alpha = alpha, loglik = -best$value,
                 se.beta = se.beta, ci.beta = ci.beta, boot = boot, se.type = se, cluster = clab,
                 fe = setNames(fe_hat, levels(uf)), fitted.values = fitted, rate = lam, offset = off,
                 weights = w, nobs_weighted = if (is.null(w)) length(Y) else sum(w),
                 xref = colMeans(X), fe_ref = mean(fe_hat),
                 ceiling = lam / (1 - alpha), n = length(Y), n_units = nu, Y = Y,
                 df = p + nu + 1L, truncated = truncated, max.support = as.integer(max.support),
                 support_binding = binding,
                 X = X, unit = as.integer(uf), levels = .cpb_xlevels(mf), contrasts = attr(Xf, "contrasts"),
                 formula = formula, fe_var = fe, converged = best$convergence == 0,
                 bias_correct = bias_correct, uncorrected = uncorrected,
                 call = match.call()),
            class = "cpb_fe")
}

#' @export
print.cpb_fe <- function(x, ...) {
  cat("CPB regression with", x$n_units, "unit fixed effects (concentrated likelihood)\n")
  if (!is.null(x$se.beta) && any(is.finite(x$se.beta))) {
    cat("Coefficients (", x$cluster, "-clustered bootstrap SEs):\n", sep = "")
    print(round(cbind(Estimate = x$coefficients,
                      `Std. Error` = x$se.beta[names(x$coefficients)]), 4))
  } else { cat("Coefficients:\n"); print(round(x$coefficients, 4)) }
  cat("\nalpha (shape parameter):", round(x$alpha, 4),
      "  median implied bound:", round(stats::median(x$ceiling), 2), "\n")
  if (identical(x$bias_correct, "jackknife"))
    cat("Estimates are split-panel jackknife bias-corrected; logLik/AIC refer to the\nuncorrected maximum-likelihood fit. See ?cpb_fe.\n")
  else
    cat("Note: alpha is subject to incidental-parameters bias for short panels; see ?cpb_fe.\n")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}

#' Predictions from a fixed-effects CPB fit
#'
#' @param object A `"cpb_fe"` object.
#' @param newdata Optional covariate profiles. With `newdata`, predictions use the
#'   average unit (the mean fixed effect), since the panel's own unit effects do
#'   not apply to new rows.
#' @param type `"response"` (the mean; for a zero-truncated fit the conditional
#'   mean E(Y | Y >= 1)), `"rate"` (the CPB rate lambda), `"link"` (the log
#'   rate), or `"ceiling"`.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict cpb_fe
#' @export
predict.cpb_fe <- function(object, newdata = NULL, type = c("response", "rate", "link", "ceiling"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    lam <- if (is.null(object$rate)) object$fitted.values else object$rate
  } else {
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss))
      stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X  <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    bad <- .ud_na_rows(X); if (any(bad)) X[bad, ] <- 0         # a missing covariate predicts NA
    lam <- as.numeric(exp(mean(object$fe) + X %*% object$coefficients))
  }
  if (is.null(newdata)) bad <- FALSE
  .ud_mask(switch(type, link = log(lam), rate = lam,
         response = .cpb_mean(lam, object$alpha, isTRUE(object$truncated)),
         ceiling = lam / (1 - object$alpha)), bad)
}
