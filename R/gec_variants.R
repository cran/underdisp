## ---------------------------------------------------------------------------
## Hurdle and zero-inflated variants of the generalized event count (GEC / Katz)
## model. The count component's dispersion delta is estimated freely, so these
## two-part models span under-, equi-, and overdispersed intensities. The zero-
## inflated fit uses the same robust direct-ML estimation as zi_cpb (seeded from
## separate fits, then joint maximization), not coordinate-wise EM.
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
#'   model (added as factor dummies).
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`.
#' @param offset Optional offset (log scale) for the intensity: a numeric vector
#'   or a column name in `data`.
#' @param cluster Optional cluster (column name or vector) for the intensity's
#'   cluster bootstrap.
#' @param se Intensity inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples.
#' @param max.support Guard on the maximum evaluated support.
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
                       se = c("none", "bootstrap"), B = 500, max.support = 500) {
  se <- match.arg(se); link <- match.arg(link)
  y <- stats::model.response(stats::model.frame(formula, data))
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  d <- as.integer(y > 0)
  int_offset <- if (is.null(offset)) NULL
                else if (is.character(offset) && length(offset) == 1L) offset else offset[y > 0]
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) {
    data[[".cluster"]] <- cluster; cluster <- ".cluster"
  }
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, family = stats::binomial(link = link))

  ## For within-unit (fixed-effects) intensity, add factor(unit) to `formula`
  ## (dummies); a concentrated zero-truncated GEC is not currently available.
  pos  <- data[y > 0, , drop = FALSE]
  ifit <- gec(formula, data = pos, truncated = TRUE, se = se, B = B, cluster = cluster,
              offset = int_offset, max.support = max.support)

  Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
  beta <- ifit$coefficients
  lambda_full <- as.numeric(exp(Xall %*% beta))
  structure(list(participation = pfit, intensity = ifit, Zpart = stats::model.matrix(pfit),
                 y = as.integer(y), p_full = as.numeric(stats::fitted(pfit)),
                 lambda_full = lambda_full, delta = ifit$delta, n = length(y), n_participate = sum(d),
                 int_beta = beta, int_xref = colMeans(Xall), int_fe_ref = 0,
                 part_beta = stats::coef(pfit), part_xref = colMeans(stats::model.matrix(pfit)),
                 cluster = if (is.character(cluster)) cluster else NULL, max.support = as.integer(max.support),
                 link = link, part_fe = part_fe,
                 formula = formula, part_formula = participation, call = match.call()),
            class = "hurdle_gec")
}

#' @method print hurdle_gec
#' @export
print.hurdle_gec <- function(x, ...) {
  disp <- if (x$delta < 0.97) "underdispersed" else if (x$delta > 1.03) "overdispersed" else "~ equidispersed"
  cat("Hurdle Generalized Event Count (Katz family)\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("Units: %d (%d participate, %.0f%%)\n", x$n, x$n_participate, 100 * x$n_participate / x$n))
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n",
              if (is.null(x$link)) "logit" else x$link))
  print(round(stats::coef(x$participation), 4))
  cat("\nIntensity (zero-truncated GEC) coefficients:\n"); print(round(x$intensity$coefficients, 4))
  cat(sprintf("Intensity dispersion delta: %.3f  [%s]\n", x$delta, disp))
  invisible(x)
}

#' Predict from a hurdle-GEC fit
#'
#' @param object A `"hurdle_gec"` object.
#' @param newdata Optional covariate profiles (intensity and participation).
#' @param type `"response"` (marginal E(Y)), `"participation"` (P(Y>0)), or
#'   `"intensity"` (E(Y | Y>0)).
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict hurdle_gec
#' @export
predict.hurdle_gec <- function(object, newdata = NULL,
                               type = c("response", "participation", "intensity"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) { p <- object$p_full; lam <- object$lambda_full } else {
    p <- as.numeric(stats::predict(object$participation, newdata = newdata, type = "response"))
    Terms <- stats::delete.response(stats::terms(object$formula))
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    lam <- as.numeric(exp(stats::model.matrix(Terms, newdata) %*% object$int_beta))
  }
  p0    <- gec_pmf_cpp(lam, object$delta, 0L, object$max.support)[, 1]
  tmean <- lam / pmax(1 - p0, 1e-8)                     # E(Y | Y > 0) for the Katz intensity
  switch(type, participation = p, intensity = tmean, response = p * tmean)
}

#' @method logLik hurdle_gec
#' @export
logLik.hurdle_gec <- function(object, ...) {
  ll <- as.numeric(stats::logLik(object$participation)) + object$intensity$loglik
  df <- length(stats::coef(object$participation)) + object$intensity$df
  attr(ll, "df") <- df; attr(ll, "nobs") <- object$n; class(ll) <- "logLik"; ll
}

#' @method nobs hurdle_gec
#' @export
nobs.hurdle_gec <- function(object, ...) object$n


#' Zero-inflated GEC (Katz-family) regression
#'
#' A structural-zero mixture with a GEC count component whose dispersion `delta`
#' is estimated freely. Estimated by robust direct maximum likelihood, seeded
#' from separate fits (a GEC on the positive counts for the count parameters and
#' a logit of the zero indicator for the inflation), exactly as in [zi_cpb()];
#' this avoids the degenerate `pi -> 0` basin that coordinate-wise EM falls into.
#'
#' @param formula Count formula, `y ~ x`.
#' @param data A data frame.
#' @param zero One-sided formula for the structural-zero logit; defaults to the
#'   count formula's right-hand side.
#' @param zero_fe Optional column name for fixed effects in the zero-inflation
#'   equation, entered as factor dummies (opt-in).
#' @param se Coefficient inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples.
#' @param cluster Optional cluster (column name or vector) for the bootstrap.
#' @param max.support Guard on the maximum evaluated support.
#' @param offset Optional exposure offset (log scale) for the count component:
#'   a numeric vector or the name of a column in `data`; the bootstrap carries
#'   it through every replicate.
#' @param maxit,tol Optimizer controls.
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
                   cluster = NULL, max.support = 500, maxit = 200, tol = 1e-6, offset = NULL) {
  se <- match.arg(se); cl <- match.call()
  if (!is.null(offset) && !(is.character(offset) && length(offset) == 1L)) {
    data[[".zi_off"]] <- as.numeric(offset); offset <- ".zi_off"    # vector -> column, subsets with data
  }
  fitfun <- function(dat) {
    mf <- stats::model.frame(formula, dat)
    y  <- as.integer(stats::model.response(mf))
    if (any(y < 0)) stop("The response must be nonnegative integer counts.")
    X  <- stats::model.matrix(formula, dat)
    off <- if (is.null(offset)) rep(0, nrow(X))
           else as.numeric(dat[[offset]][match(rownames(mf), rownames(dat))])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
    zrhs <- if (is.null(zero)) formula[-2L] else zero
    if (!is.null(zero_fe))                                        # fixed effects in the zero equation
      zrhs <- stats::reformulate(c(labels(stats::terms(zrhs)), paste0("factor(", zero_fe, ")")))
    Zt <- stats::terms(stats::update(zrhs, ~ .)); Z <- stats::model.matrix(Zt, dat)
    n  <- length(y); is0 <- y == 0; ms <- as.integer(max.support)
    pb <- ncol(X); pg <- ncol(Z)
    pf <- tryCatch(gec(formula, dat, se = "none", max.support = max.support,
                       offset = if (is.null(offset)) NULL else offset),
                   error = function(e) NULL)
    b0 <- if (!is.null(pf)) as.numeric(pf$coefficients)
          else { v <- stats::glm.fit(X, y, offset = off, family = stats::poisson())$coefficients
                 v[!is.finite(v)] <- 0; v }
    d0 <- if (!is.null(pf)) pf$delta else 0.8
    g0 <- suppressWarnings(stats::glm.fit(Z, as.numeric(is0), family = stats::binomial())$coefficients)
    g0[!is.finite(g0)] <- 0
    nll <- function(par) {
      m   <- gec_lp0_cpp(par[1:(pb + 1)], X, y, off, ms)
      pit <- as.numeric(stats::plogis(Z %*% par[(pb + 2):(pb + 1 + pg)]))
      -sum(ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
    }
    best <- NULL
    for (ld0 in unique(c(log(d0), log(0.5), 0, log(1.5)))) {
      op <- tryCatch(stats::optim(c(b0, ld0, g0), nll, method = "Nelder-Mead",
                                  control = list(maxit = 20L * maxit, reltol = tol)), error = function(e) NULL)
      if (!is.null(op) && (is.null(best) || op$value < best$value)) best <- op
    }
    if (is.null(best)) stop("zi_gec: all optimization starts failed.")
    b <- best$par[1:pb]; delta <- exp(best$par[pb + 1]); g <- best$par[(pb + 2):(pb + 1 + pg)]
    names(b) <- colnames(X); names(g) <- colnames(Z)
    list(b = b, delta = delta, g = g, ll = -best$value, it = as.integer(best$counts[1]),
         n = n, y = y, X = X, Z = Z, Zt = Zt, pb = pb, pg = pg, off = off)
  }
  r <- fitfun(data)
  lambda_full <- as.numeric(exp(r$off + r$X %*% r$b))
  se.beta <- setNames(rep(NA_real_, r$pb), names(r$b))
  se.zero <- setNames(rep(NA_real_, r$pg), names(r$g)); boot <- NULL
  if (se == "bootstrap") {
    clid <- if (is.null(cluster)) NULL
            else if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    grp  <- if (!is.null(clid)) split(seq_len(nrow(data)), clid) else NULL
    boot <- matrix(NA_real_, B, r$pb + 1 + r$pg)
    for (b in seq_len(B)) {
      idx <- if (is.null(grp)) sample.int(nrow(data), nrow(data), replace = TRUE)
             else unlist(grp[sample.int(length(grp), length(grp), replace = TRUE)], use.names = FALSE)
      fb <- tryCatch(fitfun(data[idx, , drop = FALSE]), error = function(e) NULL)
      if (!is.null(fb)) boot[b, ] <- c(fb$b, log(fb$delta), fb$g)
    }
    ok <- boot[stats::complete.cases(boot), , drop = FALSE]
    if (nrow(ok) >= 2) {
      se.beta <- setNames(apply(ok[, 1:r$pb, drop = FALSE], 2, stats::sd), names(r$b))
      se.zero <- setNames(apply(ok[, (r$pb + 2):(r$pb + 1 + r$pg), drop = FALSE], 2, stats::sd), names(r$g))
    } else warning("Too few bootstrap resamples converged for stable inference.")
  }
  Ti0 <- stats::delete.response(stats::terms(formula))
  structure(list(coefficients = r$b, delta = r$delta, dispersion = r$delta, zero_coef = r$g, zero.coefficients = r$g,
                 loglik = r$ll, iterations = r$it, df = r$pb + 1L + r$pg, n = r$n, y = r$y,
                 pi_full = as.numeric(stats::plogis(r$Z %*% r$g)), lambda_full = lambda_full,
                 se.beta = se.beta, se.zero = se.zero, boot = boot, X = r$X, Z = r$Z,
                 int_terms = Ti0, zero_terms = r$Zt,
                 int_xlev = stats::.getXlevels(Ti0, stats::model.frame(Ti0, data)),
                 zero_xlev = stats::.getXlevels(r$Zt, stats::model.frame(r$Zt, data)),
                 max.support = as.integer(max.support), formula = formula, zero_formula = zero,
                 call = cl), class = "zi_gec")
}

#' @method print zi_gec
#' @export
print.zi_gec <- function(x, ...) {
  disp <- if (x$delta < 0.97) "underdispersed" else if (x$delta > 1.03) "overdispersed" else "~ equidispersed"
  cat("Zero-Inflated Generalized Event Count (Katz family)\n")
  cat("Call:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("N: %d   logLik: %.1f\n", x$n, x$loglik))
  cat("\nCount (GEC) coefficients:\n"); print(round(x$coefficients, 4))
  cat(sprintf("Count dispersion delta: %.3f  [%s]\n", x$delta, disp))
  cat("\nZero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n")
  print(round(x$zero_coef, 4))
  cat("Mean structural-zero probability:", round(mean(x$pi_full), 3), "\n")
  invisible(x)
}

#' Predict from a zi_gec fit
#'
#' @param object A `"zi_gec"` object.
#' @param newdata Optional covariate profiles.
#' @param type `"response"` (marginal E(Y)), `"zero"` (structural-zero
#'   probability), or `"intensity"` (the count mean E of the GEC component).
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict zi_gec
#' @export
predict.zi_gec <- function(object, newdata = NULL, type = c("response", "zero", "intensity"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) { pit <- object$pi_full; lam <- object$lambda_full } else {
    Ti <- stats::delete.response(object$int_terms); Tz <- object$zero_terms
    miss <- setdiff(unique(c(all.vars(Ti), all.vars(Tz))), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    lam <- as.numeric(exp(stats::model.matrix(Ti, newdata, xlev = object$int_xlev) %*% object$coefficients))
    pit <- as.numeric(stats::plogis(stats::model.matrix(Tz, newdata, xlev = object$zero_xlev) %*% object$zero_coef))
  }
  cmean <- gec_mean_cpp(lam, object$delta, object$max.support)
  switch(type, zero = pit, intensity = cmean, response = (1 - pit) * cmean)
}

#' @method logLik zi_gec
#' @export
logLik.zi_gec <- function(object, ...) {
  ll <- object$loglik; attr(ll, "df") <- object$df; attr(ll, "nobs") <- object$n
  class(ll) <- "logLik"; ll
}

#' @method nobs zi_gec
#' @export
nobs.zi_gec <- function(object, ...) object$n

#' @method coef zi_gec
#' @export
coef.zi_gec <- function(object, ...) .zi_coef(object)
