## ---------------------------------------------------------------------------
## Hurdle and zero-inflated variants of the generalized event count (GEC / Katz)
## model. The count component's dispersion delta is estimated freely, so these
## two-part models span under-, equi-, and overdispersed intensities. The zero-
## inflated fit uses the same robust direct-ML estimation as zi_cpb (several starts that span the structural share
## of zeros, then joint maximization), not coordinate-wise EM.
## ---------------------------------------------------------------------------

#' Hurdle GEC (Katz-family) regression
#'
#' Joins a participation (hurdle) logit to a zero-truncated GEC intensity on the
#' positive counts. Unlike [hurdle_cpb()], whose intensity is fixed to the
#' underdispersed CPB, the GEC intensity estimates its own dispersion direction.
#'
#' @param formula Intensity formula, `y ~ x`.
#' @param data A data frame.
#' @param participation One-sided formula for the participation logit; defaults
#'   to the intensity's right-hand side.
#' @param part_fe Optional column name for fixed effects in the participation
#'   model (added as factor dummies), whose slopes carry the dummy logit's
#'   incidental-parameters bias in short panels (see [hurdle_cpb()]). The dummy
#'   of a unit that never participates is not counted in the degrees of freedom
#'   (see [compare_models()], Degrees of freedom). In new data such a unit, and
#'   one whose every count is positive, get the participation probability at
#'   their limits (zero and one), and an unseen unit enters at the reference
#'   unit, as in [hurdle_cpb()]; the participation equation is checked for
#'   separation (`part_separation`).
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`.
#' @param offset Optional offset (log scale) for the intensity: a numeric vector
#'   or a column name in `data`.
#' @param cluster Optional cluster (column name or vector): cluster-robust
#'   standard errors for the participation equation, and the cluster bootstrap
#'   of the intensity with `se = "bootstrap"`.
#' @param se Intensity inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples.
#' @param weights Optional frequency weights (a numeric vector or a column name),
#'   applied to both margins; see [cpb()].
#' @param cores Worker processes for the intensity bootstrap; see [cpb()].
#' @param max.support Guard on the evaluated support of the intensity; `NULL`
#'   (default) sets `max(500, 10 * max(y))`. See [gec()].
#' @return An object of class `"hurdle_gec"`.
#' @section Limitation:
#' Unlike [hurdle_cpb()], `hurdle_gec()` has no `fe` argument for a
#' concentrated fixed-effects intensity (no concentrated zero-truncated GEC is
#' implemented). For within-unit intensities, enter the unit factor as dummies
#' in `formula` (feasible for moderate unit counts), or use [hurdle_cpb()].
#' @examples
#' \donttest{
#' set.seed(3); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
#' hurdle_gec(y ~ x, data.frame(y = y, x = x, z = z), participation = ~ z)
#' }
#' @seealso [hurdle_cpb()], [gec()], [zi_gec()]
#' @export
hurdle_gec <- function(formula, data, participation = NULL, part_fe = NULL,
                       link = c("logit", "probit", "cloglog"), offset = NULL, cluster = NULL,
                       se = c("none", "bootstrap"), B = 500, weights = NULL, cores = 1L, max.support = NULL) {
  .ud_no_formula_offset(formula, participation)
  data <- .ud_drop_na_fe(data, list(part_fe))
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  se <- match.arg(se); link <- match.arg(link)
  mv <- unique(c(all.vars(formula), all.vars(if (is.null(participation)) formula[-2L] else participation), part_fe))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  cl_vals <- .ud_cluster_values(cluster, data, keep)           # checked before the rows are reduced
  if (!all(keep)) {
    if (!is.null(offset) && !is.character(offset)) offset <- offset[keep]
    if (!is.null(weights) && !is.character(weights)) weights <- weights[keep]
    data <- data[keep, , drop = FALSE]
  }
  y <- stats::model.response(stats::model.frame(formula, data))
  if (!is.numeric(y)) stop("Response must be a numeric count; got ", class(y)[1L], ".")
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  .ud_warn_all_zero(y)
  d <- as.integer(y > 0)
  int_offset <- if (is.null(offset)) NULL
                else if (is.character(offset) && length(offset) == 1L) offset else offset[y > 0]
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) {
    data[[".cluster"]] <- cl_vals; cluster <- ".cluster"
  }
  if (!is.null(weights) && !(is.character(weights) && length(weights) == 1L)) {
    if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
    data[[".hw"]] <- as.numeric(weights); weights <- ".hw"
  }
  w_all <- .ud_weights(weights, data, seq_len(nrow(data)))
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d; pdata[[".w"]] <- if (is.null(w_all)) rep(1, nrow(data)) else w_all
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, weights = .w,
                     family = stats::binomial(link = link))
  pf <- .part_fields(pfit, part_fe, data, if (is.null(cluster)) "model" else "cluster",
                     if (!is.null(cluster)) data[[cluster]], "hurdle_gec")

  ## For within-unit (fixed-effects) intensity, add factor(unit) to `formula`
  ## (dummies); a concentrated zero-truncated GEC is not currently available.
  pos  <- data[y > 0, , drop = FALSE]
  ifit <- gec(formula, data = pos, truncated = TRUE, se = se, B = B, cluster = if (se == "bootstrap") cluster,
              offset = int_offset, weights = weights, cores = cores, max.support = max.support)

  Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
  beta <- ifit$coefficients
  off_all <- if (is.null(offset)) 0 else if (is.character(offset)) as.numeric(data[[offset]]) else as.numeric(offset)
  lambda_full <- as.numeric(exp(off_all + Xall %*% beta))
  zl <- .ud_fe_zero_levels(y, data, part_fe)   # never participate: not counted (.ud_zero_levels)
  ol <- .ud_fe_one_levels(y, data, part_fe)    # always participate: predicted at the limit
  structure(c(list(participation = pfit, intensity = ifit,
                 converged = isTRUE(pfit$converged) && isTRUE(ifit$converged), Zpart = stats::model.matrix(pfit),
                 y = as.integer(y), p_full = as.numeric(stats::fitted(pfit)),
                 lambda_full = lambda_full, delta = ifit$delta, n = length(y), n_participate = sum(d),
                 int_beta = beta, int_xref = colMeans(Xall), int_fe_ref = mean(off_all),
                 part_beta = stats::coef(pfit),
                 part_xref = .fe_dummy_ref(colMeans(stats::model.matrix(pfit)), stats::model.matrix(pfit), d, part_fe,
                                           .ident_binary),
                 loglik = as.numeric(stats::logLik(pfit)) + ifit$loglik,
                 df = length(stats::coef(pfit)) - length(zl) + ifit$df, part_zero_units = zl, part_one_units = ol,
                 se.type = se, cluster = if (is.character(cluster)) cluster else NULL, max.support = ifit$max.support,
                 link = link, part_fe = part_fe,
                 weights = w_all, nobs_weighted = if (is.null(w_all)) length(y) else sum(w_all),
                 formula = formula, part_formula = participation, call = match.call()), pf),
            class = "hurdle_gec")
}

#' @method print hurdle_gec
#' @export
print.hurdle_gec <- function(x, ...) {
  disp <- .gec_direction(x)
  cat("Hurdle Generalized Event Count (Katz family)\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("Units: %d (%d participate, %.0f%%)\n", x$n, x$n_participate, 100 * x$n_participate / x$n))
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n",
              if (is.null(x$link)) "logit" else x$link))
  pc <- stats::coef(x$participation); print(round(pc[!.is_fe_col(names(pc), x$part_fe)], 4))
  cat("\nIntensity (zero-truncated GEC) coefficients:\n"); print(round(x$intensity$coefficients, 4))
  cat(sprintf("Intensity dispersion delta: %.3f  [%s]\n", x$delta, disp))
  for (nt in .sep_note(x$part_separation, "participation")) cat(nt, "\n")
  .ud_zero_df_note(x$part_zero_units, "participation dummy", "participation dummies")
  invisible(x)
}

#' Predict from a hurdle-GEC fit
#'
#' @param object A `"hurdle_gec"` object.
#' @param newdata Optional covariate profiles (intensity and participation). A
#'   `part_fe` unit the fit has not seen, or every row when that column is
#'   absent, enters at the reference unit, and a unit at a limit at that limit
#'   (see `part_fe` in [hurdle_gec()]).
#' @param type `"response"` (marginal E(Y)), `"participation"` (P(Y>0)), or
#'   `"intensity"` (E(Y | Y>0)).
#' @param offset Optional offset (log scale) of the intensity for the `newdata`
#'   rows; see [predict.cpb()].
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict hurdle_gec
#' @export
predict.hurdle_gec <- function(object, newdata = NULL,
                               type = c("response", "participation", "intensity"), offset = NULL,
                               se.fit = FALSE, interval = c("none", "confidence"),
                               level = 0.95, ...) {
  type <- match.arg(type)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  if (is.null(newdata)) { p <- object$p_full; lam <- object$lambda_full } else {
    p <- .part_at_limit(.ud_part_prob(object, newdata), object, newdata)
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    off <- .ud_newdata_offset(object$intensity, newdata, offset)
    lam <- as.numeric(exp(off + .ud_newdata_matrix(Terms, newdata, object$intensity$levels, object$intensity$contrasts,
                                                   names(object$int_beta)) %*% object$int_beta))
  }
  bad <- is.na(lam) | is.na(p); lam[bad] <- 1                 # a missing covariate predicts NA
  p0    <- gec_pmf_cpp(lam, object$delta, 0L, object$max.support)[, 1]
  tmean <- gec_mean_cpp(lam, object$delta, object$max.support) / pmax(1 - p0, 1e-12)   # exact E(Y | Y > 0)
  .ud_mask(as.numeric(switch(type, participation = p, intensity = tmean, response = p * tmean)), bad)
}

#' @method logLik hurdle_gec
#' @export
logLik.hurdle_gec <- function(object, ...) {
  ll <- as.numeric(stats::logLik(object$participation)) + object$intensity$loglik
  attr(ll, "df") <- .ud_twopart_df(object); attr(ll, "nobs") <- nobs(object); class(ll) <- "logLik"; ll
}

#' @method nobs hurdle_gec
#' @export
nobs.hurdle_gec <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted


#' Zero-inflated GEC (Katz-family) regression
#'
#' A structural-zero mixture with a GEC count component whose dispersion `delta`
#' is estimated freely. Estimated by direct maximum likelihood from several
#' starts, as in [zi_cpb()]: the single-equation GEC with a logit of the zero
#' indicator for the inflation, at four dispersions; the zero-inflated Poisson's
#' coefficients; and a low-inflation start. The best is polished, which avoids
#' the degenerate `pi -> 0` basin that coordinate-wise EM falls into.
#'
#' @param formula Count formula, `y ~ x`.
#' @param data A data frame.
#' @param zero One-sided formula for the structural-zero logit; defaults to the
#'   count formula's right-hand side.
#' @param zero_fe Optional column name for fixed effects in the zero-inflation
#'   equation, entered as factor dummies (opt-in). The dummy of a unit with no
#'   positive count is not counted in the degrees of freedom (see
#'   [compare_models()], Degrees of freedom); in new data such a unit gets the
#'   structural-zero probability one, and one whose every count is positive
#'   zero.
#' @param se Coefficient inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples.
#' @param cluster Optional cluster (column name or vector) for the bootstrap.
#' @param max.support Guard on the evaluated support of the count component;
#'   `NULL` (default) sets `max(500, 10 * max(y))`. See [gec()].
#' @param offset Optional exposure offset (log scale) for the count component:
#'   a numeric vector or the name of a column in `data`; the bootstrap carries
#'   it through every replicate.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [cpb()].
#' @param cores Worker processes for the bootstrap; see [cpb()].
#' @param maxit,tol Optimizer controls.
#' @section Separation:
#' The fit checks its zero-inflation equation for separation as [zi_cpb()] does
#' (see its Separation section). It warns when the data put no useful bound on a
#' covariate coefficient, and stores the check in the component `separation`.
#' @return An object of class `"zi_gec"`.
#' @examples
#' \donttest{
#' set.seed(4); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
#' zi_gec(y ~ x, data.frame(y = y, x = x, z = z), zero = ~ z)
#' }
#' @seealso [zi_cpb()], [gec()], [hurdle_gec()]
#' @export
zi_gec <- function(formula, data, zero = NULL, zero_fe = NULL, se = c("none", "bootstrap"), B = 500,
                   cluster = NULL, max.support = NULL, maxit = 200, tol = 1e-6, offset = NULL,
                   weights = NULL, cores = 1L) {
  .ud_no_formula_offset(formula, zero)
  data <- .ud_drop_na_fe(data, list(zero_fe))
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  se <- match.arg(se); cl <- match.call()
  .ud_check_whole(max.support, "max.support", null_ok = TRUE)
  if (se == "bootstrap") .ud_check_whole(B, "B", lower = 2)
  offset_name <- .ud_offset_name(offset)                        # before a vector becomes a column
  if (!is.null(offset) && !(is.character(offset) && length(offset) == 1L)) {
    if (length(offset) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    data[[".zi_off"]] <- as.numeric(offset); offset <- ".zi_off"    # vector -> column, subsets with data
  }
  if (!is.null(weights) && !(is.character(weights) && length(weights) == 1L)) {
    if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
    data[[".zi_w"]] <- as.numeric(weights); weights <- ".zi_w"
  }
  if (!is.null(cluster) && se != "bootstrap")
    warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
  zrhs0 <- if (is.null(zero)) formula[-2L] else zero
  mv <- unique(c(all.vars(formula), all.vars(zrhs0), zero_fe, offset, weights))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  cl_vals <- .ud_cluster_values(cluster, data, keep)           # checked before the rows are reduced
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) cluster <- cl_vals
  data <- data[keep, , drop = FALSE]
  if (is.null(max.support))                                     # the replicates keep the fit's guard
    max.support <- max(500L, 10L * max(stats::model.response(stats::model.frame(formula, data)), 0L))
  ## `start`: c(b, log(delta), g) on the original scale, from which only the polish runs (the restart from the
  ## separation check, .zi_restart)
  fitfun <- function(dat, start = NULL) {
    mf <- stats::model.frame(formula, dat, drop.unused.levels = TRUE)
    y  <- stats::model.response(mf)
    if (!is.numeric(y)) stop("Response must be a numeric count; got ", class(y)[1L], ".")
    if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
    .ud_warn_all_zero(y)
    y  <- as.integer(y)
    X  <- stats::model.matrix(formula, dat); .ud_rank_check(X, "count design")
    off <- if (is.null(offset)) rep(0, nrow(X)) else as.numeric(dat[[offset]])
    w   <- .ud_w1(.ud_weights(weights, dat, seq_len(nrow(dat))), nrow(X))
    zrhs <- zrhs0
    if (!is.null(zero_fe))                                        # fixed effects in the zero equation
      zrhs <- stats::reformulate(c(labels(stats::terms(zrhs)), paste0("factor(", zero_fe, ")")))
    Zt <- stats::terms(stats::update(zrhs, ~ .)); Z <- stats::model.matrix(Zt, dat); .ud_rank_check(Z, "zero design")
    n  <- length(y); is0 <- y == 0; ms <- as.integer(max.support)
    pb <- ncol(X); pg <- ncol(Z)
    pf <- if (is.null(start))
            tryCatch(.ud_quiet_guard(gec(formula, dat, se = "none", max.support = max.support,
                                         offset = if (is.null(offset)) NULL else offset, weights = weights)),
                     error = function(e) NULL)
    b0 <- if (!is.null(pf)) as.numeric(pf$coefficients)
          else { v <- stats::glm.fit(X, y, offset = off, family = stats::poisson())$coefficients
                 v[!is.finite(v)] <- 0; v }
    d0 <- if (!is.null(pf)) pf$delta else 0.8
    g0 <- suppressWarnings(stats::glm.fit(Z, as.numeric(is0), weights = w, family = stats::binomial())$coefficients)
    g0[!is.finite(g0)] <- 0
    nll <- function(par) {
      m   <- gec_lp0_cpp(par[1:(pb + 1)], X, y, off, ms)
      pit <- as.numeric(stats::plogis(Z %*% par[(pb + 2):(pb + 1 + pg)]))
      -sum(w * ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
    }
    ## the optimizer works on the covariates of both equations scaled to unit
    ## (weighted) standard deviation, and the best start is polished by chained
    ## Nelder-Mead restarts; the coefficients are mapped back at the end
    sc <- c(.ud_colscale(X, w), 1, .ud_colscale(Z, w))
    nll_s <- function(v) nll(v / sc)
    ## starts: the separate-stage seeds at four dispersions, and the zero-inflated
    ## Poisson (delta = 1, inside the GEC family) with its neighbours (.zi_zip_seeds)
    zs <- if (is.null(start))
            tryCatch(.zi_zip_seeds(X, Z, y, off, w, support = "ceiling"), error = function(e) NULL)   # GEC: King's support
    ## and a low-inflation start: the single-equation fit with a constant
    ## structural-zero probability of 0.02 (see .zi_low_start in zi.R)
    icz <- match("(Intercept)", colnames(Z))
    glow <- if (is.na(icz)) NULL else replace(rep(0, pg), icz, stats::qlogis(0.02))
    starts <- if (!is.null(start)) list(start * sc)
              else c(lapply(unique(c(log(d0), log(0.5), 0, log(1.5))), function(ld0) c(b0, ld0, g0) * sc),
                     if (!is.null(zs)) lapply(c(log(0.8), 0, log(1.25)), function(ld0) c(zs$b, ld0, zs$g) * sc),
                     if (!is.null(glow)) lapply(unique(c(log(d0), 0)), function(ld0) c(b0, ld0, glow) * sc))
    best <- NULL
    for (st in starts) {
      op <- tryCatch(stats::optim(st, nll_s, method = "Nelder-Mead",
                                  control = list(maxit = 20L * maxit, reltol = tol)), error = function(e) NULL)
      if (!is.null(op) && is.finite(op$value) && (is.null(best) || op$value < best$value)) best <- op
    }
    if (is.null(best)) stop("zi_gec: all optimization starts failed.")
    best <- .ud_nm_polish(best, nll_s, 20L * maxit, tol)
    best <- .ud_bfgs_pass(best, nll_s)
    best$par <- best$par / sc
    b <- best$par[1:pb]; delta <- exp(best$par[pb + 1]); g <- best$par[(pb + 2):(pb + 1 + pg)]
    names(b) <- colnames(X); names(g) <- colnames(Z)
    list(b = b, delta = delta, g = g, ll = -best$value, it = as.integer(best$counts[1]),
         n = n, y = y, X = X, Z = Z, Zt = Zt, pb = pb, pg = pg, off = off, w = w,
         converged = best$convergence == 0)
  }
  ## the fit, its separation check, and a restart of the polish from any better point it finds (see ?zi_cpb,
  ## Separation); the bootstrap replicates run the same
  p0_of <- function(r) gec_lp0_cpp(c(r$b, log(r$delta)), r$X, r$y, r$off, as.integer(max.support))[, 2]
  fit_full <- function(dat) {
    r <- fitfun(dat)
    chk <- .zi_restart(r,
      check = function(r) .zi_separation(r$Z, r$g, r$y, p0_of(r), if (is.null(weights)) NULL else r$w,
                                         .ud_fe_group(dat, zero_fe), .zi_fe_cols(r$Z, zero_fe)),
      refit = function(r, g) fitfun(dat, start = c(r$b, log(r$delta), g)),
      loglik = function(r) r$ll)
    r <- chk$fit; r$separation <- chk$separation
    r
  }
  r <- fit_full(data); sep <- r$separation
  .zi_warn_separation(sep, "zi_gec")
  lambda_full <- as.numeric(exp(r$off + r$X %*% r$b)); p0_full <- p0_of(r)
  se.beta <- setNames(rep(NA_real_, r$pb), names(r$b))
  se.zero <- setNames(rep(NA_real_, r$pg), names(r$g)); boot <- NULL
  if (se == "bootstrap") {                        # rows or clusters (.ud_boot); see .zi_cpb_boot for the layout
    clid <- if (is.null(cluster)) NULL
            else if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    zcov <- setdiff(names(r$g), .zi_fe_cols(r$Z, zero_fe))
    boot <- .ud_boot(nrow(data), clid, function(idx, blk) {
      fb <- .ud_quiet_guard(fit_full(data[idx, , drop = FALSE]))
      if (isTRUE(fb$converged)) c(fb$b, fb$g[zcov], fb$delta) else NULL
    }, c(paste0("count:", names(r$b)), paste0("zero:", zcov), "delta"), B, cores,
    label = if (is.null(clid)) "obs" else .ud_cluster_label(cluster))
    se.beta <- .ud_boot_stats(boot, paste0("count:", names(r$b)), r$b, normal = TRUE)$se
    se.zero[zcov] <- .ud_boot_stats(boot, paste0("zero:", zcov), r$g[zcov], normal = TRUE)$se
  }
  Ti0 <- stats::delete.response(stats::terms(formula))
  zl <- .ud_fe_zero_levels(r$y, data, zero_fe)   # their zero dummies at the limit: not counted
  ol <- .ud_fe_one_levels(r$y, data, zero_fe)    # every count positive: the dummy at the other limit, in new data
  structure(list(coefficients = r$b, delta = r$delta, dispersion = r$delta, zero_coef = r$g, zero.coefficients = r$g,
                 loglik = r$ll, iterations = r$it, df = r$pb + 1L + r$pg - length(zl), zero_fe_zero_units = zl,
                 zero_fe_one_units = ol,
                 n = r$n, y = r$y,
                 offset = r$off, offset_name = offset_name,
                 pi_full = as.numeric(stats::plogis(r$Z %*% r$g)), lambda_full = lambda_full, p0_full = p0_full,
                 separation = sep,
                 se.beta = se.beta, se.zero = se.zero, boot = boot, X = r$X, Z = r$Z,
                 weights = if (is.null(weights)) NULL else r$w,
                 nobs_weighted = if (is.null(weights)) r$n else sum(r$w), converged = r$converged,
                 int_terms = Ti0, zero_terms = r$Zt,
                 int_xlev = stats::.getXlevels(Ti0, stats::model.frame(Ti0, data, drop.unused.levels = TRUE)),
                 zero_xlev = stats::.getXlevels(r$Zt, stats::model.frame(r$Zt, data, drop.unused.levels = TRUE)),
                 int_contrasts = attr(stats::model.matrix(Ti0, data), "contrasts"),
                 zero_contrasts = attr(stats::model.matrix(r$Zt, data), "contrasts"),
                 xref = colMeans(r$X), zref = .fe_dummy_ref(colMeans(r$Z), r$Z, r$y, zero_fe, .ident_binary),
                 fe = NULL, zero_fe = zero_fe, se.type = se,
                 max.support = as.integer(max.support), formula = formula, zero_formula = zero,
                 call = cl), class = "zi_gec")
}

#' @method print zi_gec
#' @export
print.zi_gec <- function(x, ...) {
  disp <- .gec_direction(x)
  cat("Zero-Inflated Generalized Event Count (Katz family)\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("N: %d   logLik: %.1f\n", x$n, x$loglik))
  cat("\nCount (GEC) coefficients:\n"); print(round(x$coefficients, 4))
  cat(sprintf("Count dispersion delta: %.3f  [%s]\n", x$delta, disp))
  cat("\nZero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n")
  print(round(x$zero_coef, 4))
  cat("Mean structural-zero probability:", round(mean(x$pi_full), 3), "\n")
  .ud_zero_df_note(x$zero_fe_zero_units, "zero-equation dummy", "zero-equation dummies")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  .zi_separation_note(x)
  invisible(x)
}

#' Predict from a zi_gec fit
#'
#' @param object A `"zi_gec"` object.
#' @param newdata Optional covariate profiles. A `zero_fe` unit the fit has not
#'   seen, or every row when that column is absent, enters at the zero
#'   equation's reference unit; a unit with no positive count gets the
#'   structural-zero probability one, and one whose every count is positive
#'   zero (see [compare_models()], Degrees of freedom).
#' @param type `"response"` (marginal E(Y)), `"zero"` (structural-zero
#'   probability), or `"intensity"` (the count mean E of the GEC component).
#' @param offset Optional offset (log scale) of the count component for the
#'   `newdata` rows; see [predict.cpb()].
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict zi_gec
#' @export
predict.zi_gec <- function(object, newdata = NULL, type = c("response", "zero", "intensity"), offset = NULL,
                           se.fit = FALSE, interval = c("none", "confidence"),
                           level = 0.95, ...) {
  type <- match.arg(type)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  if (is.null(newdata)) { pit <- object$pi_full; lam <- object$lambda_full } else {
    Ti <- stats::delete.response(object$int_terms); Tz <- object$zero_terms
    miss <- setdiff(unique(c(all.vars(Ti), all.vars(Tz))), c(names(newdata), object$zero_fe))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    off <- .ud_newdata_offset(object, newdata, offset)
    lam <- as.numeric(exp(off + .ud_newdata_matrix(Ti, newdata, object$int_xlev, object$int_contrasts,
                                                   names(object$coefficients)) %*% object$coefficients))
    pit <- as.numeric(stats::plogis(.zi_gec_zero_eta(object, newdata)))
    pit <- .zero_at_limit(pit, object, newdata)
  }
  bad <- is.na(lam) | is.na(pit); lam[bad] <- 1               # a missing covariate predicts NA
  cmean <- gec_mean_cpp(lam, object$delta, object$max.support)
  .ud_mask(switch(type, zero = pit, intensity = cmean, response = (1 - pit) * cmean), bad)
}

## the zero equation's linear predictor in new data; a zero_fe unit the fit has not seen (or no zero_fe column)
## enters at the reference unit, and with `all_factors` (held-out scoring) an unseen level of any factor
.zi_gec_zero_eta <- function(object, newdata, all_factors = FALSE) {
  if (!.ud_has_unseen(object$zero_xlev, newdata, object$zero_fe, all_factors))
    return(as.numeric(.ud_newdata_matrix(object$zero_terms, newdata, object$zero_xlev, object$zero_contrasts,
                                         names(object$zero_coef)) %*% object$zero_coef))
  .ud_dummy_eta(object$zero_terms, object$zero_xlev, object$zero_coef, newdata, object$zero_contrasts,
                object$zref, object$zero_fe, all_factors)
}

#' @method logLik zi_gec
#' @export
logLik.zi_gec <- function(object, ...) {
  ll <- object$loglik; attr(ll, "df") <- object$df; attr(ll, "nobs") <- nobs(object)
  class(ll) <- "logLik"; ll
}

#' @method nobs zi_gec
#' @export
nobs.zi_gec <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted

#' @method coef zi_gec
#' @export
coef.zi_gec <- function(object, ...) .zi_coef(object)
