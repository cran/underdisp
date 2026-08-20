## ---------------------------------------------------------------------------
## Family comparison: adjudicate hard-ceiling (CPB) vs soft-tail (COM-Poisson)
## against the Poisson / negative binomial defaults, on fit and on calibration.
## ---------------------------------------------------------------------------

#' Compare count-model dispersion across the model family
#'
#' Fits the Poisson and negative binomial defaults alongside the two underdispersed
#' workhorses---the soft-tail Conway--Maxwell--Poisson (fit natively, no external
#' dependency) and the hard-ceiling CPB---and returns a comparison on both fit
#' (log-likelihood, AIC, BIC) and calibration (mean logarithmic score and ranked
#' probability score, lower is better, plus the fitted share of zeros against the
#' observed share). Whether the tail is better described by an accelerated decay
#' (COM-Poisson) or a binding ceiling (CPB) is a testable question, not an
#' assumption; the scoring rules are the calibration counterpart to the
#' information criteria. The CPB's ceiling-exceedance share (observations above
#' the implied ceiling) is reported as a falsification check on the hard bound.
#'
#' With `hurdle = TRUE` and zeros present, a hurdle-CPB is added, so a
#' zero-inflated underdispersed process can be compared to the single-equation
#' models on the same footing.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param max.support Passed to [cpb()].
#' @param hurdle If `TRUE` and the data contain zeros, also fit and score a
#'   [hurdle_cpb()].
#' @param zi If `TRUE` and the data contain zeros, also fit and score a
#'   [zi_cpb()]. The mixture EM is slow on large panels, so it is a separate
#'   opt-in from `hurdle`.
#' @return A list with `table` (a data frame of `df`, `logLik`, `AIC`, `BIC`,
#'   `logscore`, `rps`, and `zero_fit` per model), `obs_zero` (the observed zero
#'   share), `nu` (the COM-Poisson dispersion, `>1` = underdispersion, `NA` if the
#'   COM-Poisson fit failed), `alpha` (the CPB shape parameter),
#'   `ceiling_exceedance` (share of observations above the CPB ceiling), and
#'   `cpb_ok` (`FALSE` if the single-equation CPB could not satisfy its feasibility
#'   constraint on the data, e.g. under heavy zero-inflation with a wide count
#'   range --- itself a signal that a hurdle or zero-inflated model is needed).
#' @examples
#' set.seed(1); x <- rnorm(120)
#' N <- pmax(round(exp(1.6 + 0.4 * x) / 0.5), 1); y <- rbinom(120, N, 0.5)
#' compare_dispersion(y ~ x, data = data.frame(y = y, x = x))$table
#' @seealso [ud_screen()], [cpb()], [score()]
#' @export
compare_dispersion <- function(formula, data, max.support = 500, hurdle = FALSE, zi = FALSE) {
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

  rows <- list(); nu <- NA_real_; alpha <- NA_real_; exceed <- NA_real_

  ## Poisson
  pois <- stats::glm(formula, data = data, family = stats::poisson())
  p <- length(stats::coef(pois)); llp <- as.numeric(stats::logLik(pois))
  Pp <- t(vapply(stats::fitted(pois), function(m) stats::dpois(0:kmax, m), numeric(kmax + 1)))
  rows$Poisson <- c(df = p, logLik = llp, AIC = aic(llp, p), BIC = bic(llp, p), scr(Pp))

  ## Negative binomial
  nb <- tryCatch(suppressWarnings(MASS::glm.nb(formula, data = data)), error = function(e) NULL)
  if (!is.null(nb)) {
    lln <- as.numeric(stats::logLik(nb))
    Pn <- t(vapply(stats::fitted(nb), function(m) stats::dnbinom(0:kmax, mu = m, size = nb$theta),
                   numeric(kmax + 1)))
    rows$NegBinomial <- c(df = p + 1, logLik = lln, AIC = aic(lln, p + 1), BIC = bic(lln, p + 1), scr(Pn))
  }

  ## COM-Poisson (native, via count_reg -- no external dependency, and with the
  ## per-observation pmf so the proper scores are reported like every other row)
  cc <- tryCatch(suppressWarnings(count_reg(formula, data = data, family = "compois", se = "none")),
                 error = function(e) NULL)
  if (!is.null(cc) && is.finite(cc$loglik) && cc$loglik > -1e9) {
    famc <- .count_fam("compois"); rate <- exp(cc$linear.predictors)
    Pcmp <- t(vapply(rate, function(l) famc$pvec(l, cc$theta, kmax), numeric(kmax + 1)))
    rows$`COM-Poisson` <- c(df = cc$df, logLik = cc$loglik, AIC = aic(cc$loglik, cc$df),
                            BIC = bic(cc$loglik, cc$df), scr(Pcmp))
    nu <- cc$theta
  }

  ## CPB (untruncated)
  cf <- tryCatch(suppressWarnings(cpb(formula, data = data, truncated = FALSE, se = "none",
                                      max.support = max.support)), error = function(e) NULL)
  cpb_ok <- !is.null(cf) && is.finite(cf$loglik) && cf$loglik > -1e9
  if (cpb_ok) {
    lambda <- as.numeric(exp(cf$X %*% cf$coefficients))
    Pc <- t(vapply(lambda, .cpb_pmf1, numeric(kmax + 1), alpha = cf$alpha, kmax = kmax, truncated = FALSE))
    rows$CPB <- c(df = cf$df, logLik = cf$loglik, AIC = aic(cf$loglik, cf$df),
                  BIC = bic(cf$loglik, cf$df), scr(Pc))
    alpha <- cf$alpha; exceed <- mean(cf$Y > floor(cf$ceiling))
  }

  ## GEC / Katz family: dispersion estimated freely (spans under/equi/over)
  gf <- tryCatch(suppressWarnings(gec(formula, data = data, se = "none", max.support = max.support)),
                 error = function(e) NULL)
  gec_delta <- NA_real_
  if (!is.null(gf) && is.finite(gf$loglik) && gf$loglik > -1e9) {
    Pg <- gec_pmf_cpp(as.numeric(exp(gf$X %*% gf$coefficients)), gf$delta, kmax, gf$max.support)
    rows$GEC <- c(df = gf$df, logLik = gf$loglik, AIC = aic(gf$loglik, gf$df),
                  BIC = bic(gf$loglik, gf$df), scr(Pg))
    gec_delta <- gf$delta
  }

  ## optional hurdle-CPB (only meaningful with zeros)
  if (isTRUE(hurdle) && obs_zero > 0) {
    hc <- tryCatch(hurdle_cpb(formula, data = data, se = "none"), error = function(e) NULL)
    if (!is.null(hc)) {
      dfh <- length(stats::coef(hc$participation)) + hc$intensity$df
      llh <- as.numeric(stats::logLik(hc$participation)) + hc$intensity$loglik
      Ph  <- .pmf_and_y(hc, kmax)$P
      rows$`Hurdle-CPB` <- c(df = dfh, logLik = llh, AIC = aic(llh, dfh), BIC = bic(llh, dfh), scr(Ph))
    }
  }
  ## optional ZI-CPB (mixture; EM is slow on large panels, so a separate opt-in)
  if (isTRUE(zi) && obs_zero > 0) {
    zc <- tryCatch(zi_cpb(formula, data = data), error = function(e) NULL)
    if (!is.null(zc)) {
      Pz <- .pmf_and_y(zc, kmax)$P
      rows$`ZI-CPB` <- c(df = zc$df, logLik = zc$loglik, AIC = aic(zc$loglik, zc$df),
                         BIC = bic(zc$loglik, zc$df), scr(Pz))
    }
  }

  list(table = as.data.frame(do.call(rbind, rows)), obs_zero = obs_zero,
       nu = nu, alpha = alpha, gec_delta = gec_delta, ceiling_exceedance = exceed, cpb_ok = cpb_ok)
}

#' Compare fitted underdispersed-count models
#'
#' Compares fitted models from this package --- `cpb`, `cpb_fe`, `hurdle_cpb`, or
#' `zi_cpb`, in any combination of fixed effects and robust/clustered standard
#' errors --- on information criteria and, where a predicted distribution is
#' available, on proper scores. Unlike [zi_test()] this is not a hypothesis test:
#' the models need not be nested, so it is the appropriate tool for the
#' non-nested hurdle-versus-mixture choice. Passing an object that is not a
#' fitted model from this package is an error, and models fit on different data
#' raise a warning.
#'
#' @param ... Two or more fitted models (`cpb`, `cpb_fe`, `hurdle_cpb`,
#'   `zi_cpb`), optionally named.
#' @return A data frame with `df`, `logLik`, `AIC`, `BIC`, and (where the
#'   predicted distribution is available) `logscore` and `rps`, one row per
#'   model, ordered by AIC.
#' @seealso [zi_test()], [compare_dispersion()]
#' @examples
#' \donttest{
#' set.seed(1); n <- 800; x <- rnorm(n); z <- rnorm(n)
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
  getn <- function(m) as.integer(if (!is.null(m$n)) m$n else length(m$y))
  ns <- vapply(models, getn, integer(1))
  if (length(unique(ns)) > 1L)
    warning("Models were fit on different numbers of observations (", paste(unique(ns), collapse = ", "),
            "); AIC/BIC and scores are only comparable on the same data.")
  info <- function(m) {
    if (inherits(m, c("hurdle_cpb", "hurdle_gec"))) {
      idf <- if (!is.null(m$intensity$df)) m$intensity$df else length(m$intensity$coefficients) + 1L
      df  <- length(stats::coef(m$participation)) + idf
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
  as.data.frame(tab)[order(tab[, "AIC"]), , drop = FALSE]
}
