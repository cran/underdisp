## ---------------------------------------------------------------------------
## Family comparison: adjudicate hard-ceiling (CPB) vs soft-tail (COM-Poisson)
## against the Poisson / negative binomial defaults, on fit and on calibration.
## ---------------------------------------------------------------------------

#' Compare count-model dispersion across the model family
#'
#' Fits the Poisson and negative binomial defaults alongside the underdispersed
#' models---the soft-tail Conway--Maxwell--Poisson (fit natively, no external
#' dependency), the hard-ceiling CPB, and the GEC, which spans both
#' directions---and returns a comparison on both fit
#' (log-likelihood, AIC, BIC) and calibration (mean logarithmic score and ranked
#' probability score, lower is better, plus the fitted share of zeros against the
#' observed share). Whether the tail is better described by an accelerated decay
#' (COM-Poisson) or a hard ceiling (CPB) is a question for the data's fit and
#' calibration, not an assumption; the scoring rules are the calibration
#' counterpart to the information criteria. Either family's better fit
#' describes the distribution, not the process behind it. The CPB's bound can
#' be checked only on held-out counts, with [ceiling_exceedance()], and only
#' against a ceiling that is too low: on the data it was fit to every feasible
#' CPB respects its bound, since its likelihood forbids exceedance.
#'
#' With `hurdle = TRUE` and zeros present, a hurdle-CPB is added, so a
#' zero-inflated underdispersed process can be compared to the single-equation
#' models on the same footing.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param max.support Passed to the CPB and GEC fits; the default `NULL` lets
#'   each set its own (see [cpb()]); on the floor support it is 500.
#' @param hurdle If `TRUE` and the data contain zeros, also fit and score a
#'   [hurdle_cpb()].
#' @param zi If `TRUE` and the data contain zeros, also fit and score a
#'   [zi_cpb()]. The mixture EM is slow on large panels, so it is a separate
#'   opt-in from `hurdle`.
#' @param support The CPB's support for the CPB-family fits, `"ceiling"` (King's
#'   CPB, the default) or `"floor"`; see [cpb()].
#' @return A list with `table` (a data frame of `df`, `logLik`, `AIC`, `BIC`,
#'   `logscore`, `rps`, `zero_fit`, and `converged` per model), `obs_zero` (the
#'   observed zero share), `nu` (the COM-Poisson dispersion, `>1` =
#'   underdispersion, `NA` if the COM-Poisson fit failed), `alpha` (the CPB shape
#'   parameter), `gec_delta` (the GEC's dispersion, `NA` if its fit failed), and
#'   `cpb_ok` (`FALSE` if the single-equation CPB could not satisfy its feasibility
#'   constraint on the data, e.g. under heavy zero-inflation with a wide count
#'   range: the CPB's support cannot hold the counts, which another family, or
#'   a hurdle or zero-inflated model, may fit).
#' @examples
#' set.seed(1); x <- rnorm(120)
#' N <- pmax(round(exp(1.6 + 0.4 * x) / 0.5), 1); y <- rbinom(120, N, 0.5)
#' compare_dispersion(y ~ x, data = data.frame(y = y, x = x))$table
#' @seealso [ud_screen()], [cpb()], [score()]
#' @export
compare_dispersion <- function(formula, data, max.support = NULL, hurdle = FALSE, zi = FALSE,
                               support = c("ceiling", "floor")) {
  .ud_no_formula_offset(formula)
  support <- match.arg(support)
  .ud_check_whole(max.support, "max.support", null_ok = TRUE)
  if (is.null(max.support) && .cpb_fl(support)) max.support <- 500L
  y <- stats::model.response(stats::model.frame(formula, data))
  n <- length(y); kmax <- max(y); obs_zero <- mean(y == 0)
  aic <- function(ll, df) -2 * ll + 2 * df
  bic <- function(ll, df) -2 * ll + df * log(n)

  ## log score, RPS, and fitted zero share from an n x (kmax+1) pmf matrix
  scr <- function(P) {
    yk <- pmin(y, kmax)
    ls <- -mean(log(pmax(P[cbind(seq_len(n), yk + 1L)], 1e-12)))
    cdf <- t(apply(P, 1, cumsum))
    ind <- outer(yk, 0:kmax, function(a, b) as.numeric(a <= b))
    c(logscore = ls, rps = mean(rowSums((cdf - ind)^2)), zero_fit = mean(P[, 1]))
  }

  rows <- list(); nu <- NA_real_; alpha <- NA_real_

  ## Poisson
  pois <- stats::glm(formula, data = data, family = stats::poisson())
  p <- length(stats::coef(pois)); llp <- as.numeric(stats::logLik(pois))
  Pp <- t(vapply(stats::fitted(pois), function(m) stats::dpois(0:kmax, m), numeric(kmax + 1)))
  rows$Poisson <- c(df = p, logLik = llp, AIC = aic(llp, p), BIC = bic(llp, p), scr(Pp), converged = isTRUE(pois$converged))

  ## Negative binomial
  nb <- tryCatch(suppressWarnings(MASS::glm.nb(formula, data = data)), error = function(e) NULL)
  if (!is.null(nb)) {
    lln <- as.numeric(stats::logLik(nb))
    Pn <- t(vapply(stats::fitted(nb), function(m) stats::dnbinom(0:kmax, mu = m, size = nb$theta),
                   numeric(kmax + 1)))
    rows$NegBinomial <- c(df = p + 1, logLik = lln, AIC = aic(lln, p + 1), BIC = bic(lln, p + 1), scr(Pn),
                          converged = isTRUE(nb$converged))
  }

  ## COM-Poisson (native, via count_reg -- no external dependency, and with the
  ## per-observation pmf so the proper scores are reported like every other row)
  cc <- tryCatch(suppressWarnings(count_reg(formula, data = data, family = "compois", se = "none")),
                 error = function(e) NULL)
  if (!is.null(cc) && is.finite(cc$loglik) && cc$loglik > -1e9) {
    famc <- .count_fam("compois"); rate <- exp(cc$linear.predictors)
    Pcmp <- t(vapply(rate, function(l) famc$pvec(l, cc$theta, kmax), numeric(kmax + 1)))
    rows$`COM-Poisson` <- c(df = cc$df, logLik = cc$loglik, AIC = aic(cc$loglik, cc$df),
                            BIC = bic(cc$loglik, cc$df), scr(Pcmp), converged = isTRUE(cc$converged))
    nu <- cc$theta
  }

  ## CPB (untruncated)
  cf <- tryCatch(suppressWarnings(cpb(formula, data = data, truncated = FALSE, se = "none",
                                      max.support = max.support, support = support)), error = function(e) NULL)
  cpb_ok <- !is.null(cf) && is.finite(cf$loglik) && cf$loglik > -1e9
  if (cpb_ok) {
    lambda <- as.numeric(exp(cf$X %*% cf$coefficients))
    Pc <- t(vapply(lambda, .cpb_pmf1, numeric(kmax + 1), alpha = cf$alpha, kmax = kmax, truncated = FALSE,
                   support = support))
    rows$CPB <- c(df = cf$df, logLik = cf$loglik, AIC = aic(cf$loglik, cf$df),
                  BIC = bic(cf$loglik, cf$df), scr(Pc), converged = isTRUE(cf$converged))
    alpha <- cf$alpha
  }

  ## GEC / Katz family: dispersion estimated freely (spans under/equi/over)
  gf <- tryCatch(suppressWarnings(gec(formula, data = data, se = "none", max.support = max.support)),
                 error = function(e) NULL)
  gec_delta <- NA_real_
  if (!is.null(gf) && is.finite(gf$loglik) && gf$loglik > -1e9) {
    Pg <- gec_pmf_cpp(as.numeric(exp(gf$X %*% gf$coefficients)), gf$delta, kmax, gf$max.support)
    rows$GEC <- c(df = gf$df, logLik = gf$loglik, AIC = aic(gf$loglik, gf$df),
                  BIC = bic(gf$loglik, gf$df), scr(Pg), converged = isTRUE(gf$converged))
    gec_delta <- gf$delta
  }

  ## optional hurdle-CPB (only meaningful with zeros)
  if (isTRUE(hurdle) && obs_zero > 0) {
    hc <- tryCatch(.ud_quiet_guard(hurdle_cpb(formula, data = data, se = "none", support = support)),
                   error = function(e) NULL)
    if (!is.null(hc)) {
      dfh <- .ud_twopart_df(hc)
      llh <- as.numeric(stats::logLik(hc$participation)) + hc$intensity$loglik
      Ph  <- .pmf_and_y(hc, kmax)$P
      rows$`Hurdle-CPB` <- c(df = dfh, logLik = llh, AIC = aic(llh, dfh), BIC = bic(llh, dfh), scr(Ph),
                             converged = isTRUE(hc$converged))
    }
  }
  ## optional ZI-CPB (mixture; EM is slow on large panels, so a separate opt-in)
  if (isTRUE(zi) && obs_zero > 0) {
    zc <- tryCatch(.ud_quiet_guard(zi_cpb(formula, data = data, support = support)), error = function(e) NULL)
    if (!is.null(zc)) {
      Pz <- .pmf_and_y(zc, kmax)$P
      rows$`ZI-CPB` <- c(df = zc$df, logLik = zc$loglik, AIC = aic(zc$loglik, zc$df),
                         BIC = bic(zc$loglik, zc$df), scr(Pz), converged = isTRUE(zc$converged))
    }
  }

  tab <- as.data.frame(do.call(rbind, rows)); tab$converged <- as.logical(tab$converged)
  list(table = tab, obs_zero = obs_zero,
       nu = nu, alpha = alpha, gec_delta = gec_delta, cpb_ok = cpb_ok)
}

#' Compare fitted underdispersed-count models
#'
#' Compares fitted models from this package --- `cpb`, `cpb_fe`, `hurdle_cpb`,
#' `zi_cpb`, `gec`, `gec_fe`, `hurdle_gec`, `zi_gec`, `count_reg`,
#' `hurdle_count`, and `zi_count`, in any combination of fixed effects and
#' robust/clustered standard errors --- on information criteria and, where a
#' predicted distribution is
#' available, on proper scores. Unlike [zi_test()] this is not a hypothesis test:
#' the models need not be nested, so it is the appropriate tool for the
#' non-nested hurdle-versus-mixture choice. Passing an object that is not a
#' fitted model from this package is an error, and models fit on different data
#' raise a warning.
#'
#' @param ... Two or more fitted models of the classes above, optionally named.
#'   With frequency weights, BIC uses the weight total ([nobs()]).
#' @return A data frame with `df`, `logLik`, `AIC`, `BIC`, and (where the
#'   predicted distribution is available) `logscore` and `rps`, one row per
#'   model, ordered by AIC, and the flags `truncated` (a zero-truncated
#'   likelihood, conditional on `Y >= 1`) and `separated` (a zero-inflation or
#'   participation equation the data do not bound, whose log-likelihood is a
#'   supremum along a divergent ridge; see [zi_cpb()], Separation). The
#'   function warns when the models mix truncated and untruncated likelihoods,
#'   whose log-likelihoods are of different events, or when a model is
#'   separated.
#' @section Degrees of freedom:
#' Each model's `df` is its own (`logLik()`'s). The hurdle and zero-inflated
#' fits do not count the parameters of a unit with no positive count: with unit
#' fixed effects in the participation or zero equation (`part_fe`, `zero_fe`)
#' or in the count equation (`fe`; dummies written into a formula are
#' covariates to the fit and are counted), such a unit's dummy and its count effect sit
#' at their limits, where every one of its rows has probability one, so it adds
#' nothing to the likelihood. As fixest leaves such units out of a
#' fixed-effects fit, these parameters are left out of `df`, so that a hurdle
#' and a zero-inflated model fit to the same panel are charged alike. The fits
#' record the units (`part_zero_units` in the hurdles; `zero_units` and
#' `zero_fe_zero_units` in [zi_cpb()]; `fe_zero_units` and `zero_fe_zero_units`
#' in [zi_count()]; `zero_fe_zero_units` in [zi_gec()]), and their print
#' methods give the number. In new data (`predict()`, [score()], [cv_score()])
#' the hurdles give such a unit's rows the participation probability at that
#' limit, zero, as the zero-inflated fits give them a count effect at its limit.
#' A unit whose every count is positive has its hurdle participation dummy at
#' the other limit, plus infinity; that dummy is counted, as the zero-inflated
#' fits count their zero-equation dummy for the unit, and in new data the
#' hurdles give the unit's rows the participation probability one (the fits
#' record these units in `part_one_units`). The zero-inflated fits give a
#' `zero_fe` unit's rows the structural-zero probability at its limit in the
#' same way: one for a unit with no positive count, zero for a unit whose every
#' count is positive (`zero_fe_one_units`). CPB fits on the floor support
#' (`support = "floor"`) count these parameters, and in new data evaluate these
#' units at the coefficients the fit stopped at.
#'
#' The rule is the two-part fits'. The single-equation fixed-effects fits
#' ([cpb_fe()], [gec_fe()], `count_reg(fe = )`) count every unit's effect, a
#' unit with no positive count included;
#' set beside a two-part fit to the same panel, such a fit is charged two AIC
#' points per unit with no positive count that the two-part fit is not, and
#' `compare_models()` says so in a message.
#' @seealso [zi_test()], [compare_dispersion()]
#' @examples
#' \donttest{
#' set.seed(1); n <- 400; x <- rnorm(n); z <- rnorm(n)
#' y <- rhurdle_cpb(n, exp(1.2 + 0.5 * x), 0.5, plogis(-0.2 + 0.8 * z))
#' d <- data.frame(y = y, x = x, z = z)
#' compare_models(hurdle = hurdle_cpb(y ~ x, data = d, participation = ~ z),
#'                zi = zi_cpb(y ~ x, data = d, zero = ~ z))
#' }
#' @export
compare_models <- function(...) {
  models <- list(...)
  if (length(models) < 2L) stop("Provide at least two fitted models to compare.")
  cls <- c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb", "gec", "gec_fe", "hurdle_gec", "zi_gec",
           "count_reg", "hurdle_count", "zi_count")
  ok  <- vapply(models, function(m) inherits(m, cls), logical(1))
  if (!all(ok))
    stop("compare_models() takes fitted models from this package (cpb, cpb_fe, hurdle_cpb, ",
         "zi_cpb, gec, gec_fe, hurdle_gec, zi_gec, count_reg, hurdle_count, zi_count); argument(s) ",
         paste(which(!ok), collapse = ", "), " are not.")
  ## the observations each model counts: the weight total for a weighted fit (nobs()), as BIC needs
  getn <- function(m) { v <- tryCatch(stats::nobs(m), error = function(e) NULL)
    as.numeric(if (!is.null(v)) v else if (!is.null(m$n)) m$n else length(m$y)) }
  ns <- vapply(models, getn, numeric(1))
  if (length(unique(ns)) > 1L)
    warning("Models were fit on different numbers of observations (", paste(unique(ns), collapse = ", "),
            "); AIC/BIC and scores are only comparable on the same data.")
  ## the same data: the responses each fit keeps, compared where they have the same length
  ys <- lapply(models, function(m) { v <- if (!is.null(m$y)) m$y else m$Y; if (is.null(v)) NULL else as.numeric(v) })
  if (!any(vapply(ys, is.null, logical(1))) && length(unique(lengths(ys))) == 1L &&
      !all(vapply(ys[-1L], function(v) identical(v, ys[[1L]]), logical(1))))
    warning("Models were fit to different responses; AIC/BIC and scores are only comparable on the same data.")
  info <- function(m) {
    if (inherits(m, c("hurdle_cpb", "hurdle_gec"))) {
      df  <- .ud_twopart_df(m)
      ll  <- as.numeric(stats::logLik(m$participation)) + m$intensity$loglik
    } else { df <- m$df; ll <- m$loglik }
    n  <- getn(m)
    sc <- tryCatch(score(m), error = function(e) c(logscore = NA_real_, rps = NA_real_))
    c(df = df, logLik = ll, AIC = -2 * ll + 2 * df, BIC = -2 * ll + df * log(n),
      logscore = unname(sc["logscore"]), rps = unname(sc["rps"]))
  }
  tab <- t(vapply(models, info, numeric(6)))
  nm  <- names(models); if (is.null(nm)) nm <- rep("", length(models))
  rn  <- ifelse(nzchar(nm), nm, vapply(models, function(m) class(m)[1L], character(1)))
  rownames(tab) <- make.unique(rn)
  out <- as.data.frame(tab)
  out$truncated <- vapply(models, function(m) isTRUE(m$truncated), logical(1))
  out$separated <- vapply(models, function(m) length(.ud_separated(m)) > 0L, logical(1))
  if (any(out$truncated) && !all(out$truncated))
    warning("Some models are zero-truncated (", paste(rownames(out)[out$truncated], collapse = ", "), "): their ",
            "likelihood is conditional on Y >= 1, so their log-likelihoods, AIC and BIC are not comparable with the ",
            "untruncated models'.", call. = FALSE)
  if (any(out$separated))
    warning("The zero-inflation or participation equation of ", paste(rownames(out)[out$separated], collapse = ", "),
            " is separated: its log-likelihood is a supremum along a divergent ridge, so its rank on AIC/BIC rests on ",
            "terms the data do not bound (see ?zi_cpb, Separation).", call. = FALSE)
  ## a single-equation fixed-effects fit counts every unit's effect, a two-part fit not those of a unit with no
  ## positive count (Degrees of freedom): say so where the two kinds meet on such units
  k1 <- vapply(models, .ud_single_zero_units, numeric(1))
  two <- vapply(models, function(m) length(c(m$part_zero_units, m$zero_units, m$zero_fe_zero_units,
                                              m$fe_zero_units)) > 0L, logical(1))
  if (any(k1 > 0) && any(two))
    message("compare_models(): ", paste(rownames(out)[k1 > 0], collapse = ", "), " count", if (sum(k1 > 0) == 1L) "s",
            " the effect of each fixed-effects unit with no positive count (", paste(unique(k1[k1 > 0]), collapse = ", "),
            " units), which the two-part fits leave out of their degrees of freedom, so ", if (sum(k1 > 0) == 1L) "its" else "their",
            " AIC carries two points per such unit that theirs does not (see ?compare_models, Degrees of freedom).")
  out[order(tab[, "AIC"]), , drop = FALSE]
}

## the number of fixed-effects units with no positive count whose effect a single-equation, untruncated
## fixed-effects fit counts in its degrees of freedom (cpb_fe(), gec_fe(), count_reg(fe = )); 0 for other fits
.ud_single_zero_units <- function(m) {
  if (isTRUE(m$truncated)) return(0)
  if (inherits(m, c("cpb_fe", "gec_fe")) && !is.null(m$fe_identified)) return(sum(!m$fe_identified))
  if (inherits(m, "count_reg") && !is.null(m$fe) && !is.null(m$X)) {
    k <- 0
    for (v in m$fe) {
      cols <- colnames(m$X)[startsWith(colnames(m$X), paste0("factor(", v, ")"))]
      if (!length(cols)) next
      lev <- as.integer(m$X[, cols, drop = FALSE] %*% seq_along(cols))
      k <- k + sum(!tapply(m$Y > 0, lev, any))
    }
    return(k)
  }
  0
}
