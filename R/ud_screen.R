## ---------------------------------------------------------------------------
## Underdispersion screening diagnostic
## ---------------------------------------------------------------------------

## Zero-truncated Poisson MLE by Fisher scoring with a log link on the rate:
## the score contribution is (y - m(lambda)) x and the expected information
## weight is m'(lambda) * lambda. Started from the positives-only Poisson GLM
## and step-halved on the log-likelihood, this converges where VGAM's IRLS can
## stop internally on dummy-heavy designs; a converged fit is the same MLE.
.ztp_fisher <- function(formula, data_pos, maxit = 200L) {
  mf <- tryCatch(stats::model.frame(formula, data_pos, drop.unused.levels = TRUE), error = function(e) NULL)
  if (is.null(mf)) return(NULL)
  X <- stats::model.matrix(attr(mf, "terms"), mf)
  y <- as.numeric(stats::model.response(mf))
  .ztp_fisher_fit(X, y, maxit)
}

## Core Fisher-scoring iteration on an explicit design; reused by the
## parametric-bootstrap threshold, which refits many simulated responses on
## one fixed design matrix.
.ztp_fisher_fit <- function(X, y, maxit = 200L) {
  gs <- tryCatch(suppressWarnings(stats::glm.fit(X, y, family = stats::poisson())),
                 error = function(e) NULL)
  if (is.null(gs)) return(NULL)
  b <- gs$coefficients; b[!is.finite(b)] <- 0
  llf <- function(lam) sum(y * log(lam) - lam - log1p(-exp(-lam)) - lgamma(y + 1))
  ## eta floored at -8: a unit whose replicate counts are all ones drives its
  ## ZTP effect to -Inf (quasi-separation; P(Y=1|Y>0) -> 1 as lambda -> 0), so
  ## the likelihood climbs forever without the floor. At the floor the unit's
  ## Pearson contribution is ~0, its correct limit.
  lam <- exp(pmin(pmax(as.numeric(X %*% b), -8), 30)); ll <- llf(lam)
  if (!is.finite(ll)) return(NULL)
  for (it in seq_len(maxit)) {
    em <- exp(-lam); m <- lam / (1 - em)
    w  <- pmax(((1 - em * (1 + lam)) / (1 - em)^2) * lam, 1e-10)
    A <- crossprod(X, X * w); g <- crossprod(X, y - m)
    step <- tryCatch(solve(A, g), error = function(e)
      tryCatch(solve(A + diag(1e-8 * max(diag(A)), ncol(A)), g),  # ridge fallback
               error = function(e2) NULL))
    if (is.null(step)) return(NULL)
    sz <- 1; ll2 <- -Inf; b2 <- b
    for (h in 1:30) {
      b2 <- b + sz * as.numeric(step)
      lam2 <- exp(pmin(pmax(as.numeric(X %*% b2), -8), 30)); ll2 <- llf(lam2)
      if (is.finite(ll2) && ll2 >= ll - 1e-8) break
      sz <- sz / 2
    }
    if (!is.finite(ll2) || ll2 < ll - 1e-8) break
    conv <- abs(ll2 - ll) < 1e-9 * (abs(ll) + 1)
    b <- b2; lam <- lam2; ll <- ll2
    if (conv) return(list(lam = lam, ll = ll, converged = TRUE, rank = qr(X)$rank))
  }
  list(lam = lam, ll = ll, converged = FALSE, rank = qr(X)$rank)
}

## Conditional dispersion + direction, given response y, fitted mu, #params k.
.ud_disp <- function(y, mu, k) {
  n <- length(y)
  pearson <- sum((y - mu)^2 / mu) / (n - k)
  w  <- (y - mu)^2 - y
  cf <- summary(stats::lm(w ~ mu + 0))$coefficients
  list(n = n, pearson = pearson, slope = cf[1, 1], t = cf[1, 3], p = cf[1, 4])
}

## Calibrated one-sided 5% threshold for the at-risk ZTP-Pearson under an
## equidispersed (zero-truncated Poisson) null: 1 - 2.27/sqrt(n_pos), fit to the
## size/power simulation (0.84 at n=200, 0.90 at n=500, 0.94 at n=1500).
.ud_ztp_threshold <- function(n_pos) 1 - 2.27 / sqrt(n_pos)

.ud_verdict_marg <- function(dd, p_nb = NA) {
  over  <- (dd$slope > 0 && dd$p < 0.10) &&
           (dd$pearson > 1.10 || (!is.na(p_nb) && p_nb < 0.05))
  under <- (dd$pearson < 0.85) || (dd$slope < 0 && dd$p < 0.10 && dd$pearson < 1.0)
  if (over) "OVERDISPERSED" else if (under) "UNDERDISPERSED" else "~equi / inconclusive"
}

#' Screen a count outcome for underdispersion
#'
#' Applies the diagnostic sequence developed in Bagozzi (2026): a *marginal* verdict
#' from the conditional Pearson statistic and a through-origin score regression (which
#' gives the direction of dispersion), a negative-binomial-versus-Poisson test for the
#' overdispersion call, and---most importantly---an *at-risk* verdict on the positive
#' counts benchmarked against a **zero-truncated Poisson**. The last step is what
#' separates genuine underdispersion from the artifact of conditioning on \eqn{Y>0}.
#' When the data are not zero-dominated it also fits the CPB and reports a
#' ceiling-exceedance diagnostic, and it fits a generalized-Poisson soft-tail
#' comparator.
#'
#' Runtime: the marginal and at-risk arms take seconds. `ztp_threshold = "bootstrap"`
#' refits the zero-truncated Poisson `ztp_boot_B` times (use `cores`); the CPB
#' comparator is fit only up to `cpb_max_n` rows, and the generalized-Poisson and
#' COM-Poisson comparators only when the mean model carries at most
#' `comp_max_par` parameters and at most `comp_max_n` rows, beyond which their
#' rows print `NA`; each is a full model fit (the two soft-tail families have no
#' concentrated fixed-effects path, so a dummy-heavy screen would take minutes).
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param run_cpb Logical; fit the CPB when the data are not zero-dominated
#'   (default `TRUE`).
#' @param cpb_max_n Skip the CPB fit above this sample size (default 3000).
#' @param run_gp Logical; fit the generalized-Poisson comparator (default
#'   `TRUE`); gated by `comp_max_par` and `comp_max_n` like the COM-Poisson.
#' @param run_comp Logical; fit the native COM-Poisson comparator (default
#'   `TRUE`). Skipped when the mean model carries more than `comp_max_par`
#'   parameters (the COM-Poisson has no concentrated fixed-effects path, so
#'   dummy-heavy screens would be slow) or when `n` exceeds `comp_max_n`.
#' @param comp_max_par,comp_max_n Parameter and sample-size gates for the
#'   generalized-Poisson and COM-Poisson comparators (defaults 30 and 5000).
#' @param ztp_threshold How to set the at-risk test's underdispersion cutoff.
#'   `"calibrated"` (default) uses the simulation-calibrated rule
#'   \eqn{1 - 2.27/\sqrt{n_+}}, whose constant is an estimated standard
#'   deviation fitted to one calibration grid (no fixed effects, and that
#'   grid's rate profile); `"bootstrap"` calibrates the cutoff on the data at
#'   hand by a parametric bootstrap of the fitted zero-truncated Poisson null
#'   (simulate `ztp_boot_B` at-risk panels at the fitted rates, refit, and take
#'   the empirical 5\% quantile of the statistic). The calibrated rule is
#'   anti-conservative in two separable regimes: when the mean model's
#'   parameter count is a nontrivial share of the positive observations
#'   (dummy-heavy fixed-effects screens, where the null's center drifts with
#'   the parameter share), and when the fitted rates concentrate at small
#'   values (roughly \eqn{\lambda} below 2, and the more severely the smaller
#'   the rates), where the null's spread exceeds the fitted constant even
#'   without fixed effects (simulated size roughly 0.07--0.10 against the
#'   nominal 0.05). The bootstrap absorbs both
#'   departures by construction and is the recommended choice in either
#'   regime.
#' @param ztp_boot_B Number of parametric-bootstrap replicates (default 199).
#' @param cores Worker processes for the parametric-bootstrap replicates
#'   (default 1); see [cpb()].
#' @param digits Printing precision.
#' @param nb_boot Whether the NB-vs-Poisson p-value is a parametric bootstrap
#'   under the fitted Poisson (`ztp_boot_B` replicates over `cores`), the
#'   package default at a boundary null (see [dispersion_test()]). `NULL`
#'   (default) bootstraps it when `ztp_threshold = "bootstrap"` and the mean
#'   model is within `comp_max_par` and `comp_max_n`; otherwise the p-value is
#'   the asymptotic boundary mixture, which is conservative for this test. The
#'   bootstrap restores the random-number state, so every other result is the
#'   same with or without it.
#' @return An object of class `"ud_screen"` with `verdict_marginal`, `verdict_atrisk`,
#'   the conditional and at-risk (ZTP-benchmarked) Pearson statistics, the NB-vs-Poisson
#'   LR test (`p_nb_method` records how its p-value was computed), a log-likelihood
#'   comparison, (when fit) the CPB alpha and
#'   ceiling-exceedance share, and the over-conditioning guard state:
#'   `atrisk_skipped` (`TRUE` when the mean model nearly saturates the positive
#'   counts, so the at-risk statistic is not computed and the printout says why),
#'   `sat_ratio` (the fitted parameter share of the positives), and
#'   `overconditioned` (`TRUE` when that share reaches 0.10, the region where the
#'   calibrated threshold is anti-conservative; the printout then flags the
#'   verdict as diagnostic rather than probative and recommends
#'   `ztp_threshold = "bootstrap"`).
#' @references King, G. (1989). Variance specification in event count models.
#'   \emph{AJPS} 33(3):762-784.
#'
#'   Bagozzi, B. E. (2026). Revisiting underdispersion in political science.
#'   Companion manuscript.
#' @examples
#' set.seed(1); x <- rnorm(250)
#' N <- pmax(round(exp(1.4 + 0.4 * x) / 0.6), 1); y <- rbinom(250, N, 0.6)
#' ud_screen(y ~ x, data = data.frame(y = y, x = x))
#' @seealso [cpb()], [compare_dispersion()]
#' @export
ud_screen <- function(formula, data, run_cpb = TRUE, cpb_max_n = 3000,
                      run_gp = TRUE, run_comp = TRUE, comp_max_par = 30,
                      comp_max_n = 5000,
                      ztp_threshold = c("calibrated", "bootstrap"),
                      ztp_boot_B = 199L, cores = 1L, digits = 3, nb_boot = NULL) {
  .ud_no_formula_offset(formula)
  ztp_threshold <- match.arg(ztp_threshold)
  mf <- model.frame(formula, data, na.action = na.omit, drop.unused.levels = TRUE)
  y  <- model.response(mf)
  if (any(y < 0) || any(y != floor(y)))
    stop("Response must be a non-negative integer count.")
  n <- length(y)

  pois <- glm(formula, data = data, family = poisson())
  X <- model.matrix(pois); k <- ncol(X); yv <- pois$y
  dd <- .ud_disp(yv, fitted(pois), k)
  nb <- tryCatch(suppressWarnings(MASS::glm.nb(formula, data = data)), error = function(e) NULL)
  ## the NB nests the Poisson at theta -> Inf, so the LR is at least zero: a
  ## negative difference is glm.nb() stopping short of that boundary
  lr_nb <- if (!is.null(nb)) max(as.numeric(2 * (logLik(nb) - logLik(pois))), 0) else NA_real_
  ## boundary-corrected: the NB nests the Poisson at theta -> Inf (a boundary), so
  ## the LR null is the 1/2 chi^2_0 + 1/2 chi^2_1 mixture -- same correction the
  ## alpha-existence test uses (methods.R). Halving keeps the two coherent.
  p_nb  <- if (!is.na(lr_nb)) 0.5 * pchisq(max(lr_nb, 0), df = 1, lower.tail = FALSE) else NA_real_
  ## In bootstrap mode the NB p-value is a parametric bootstrap under the fitted
  ## Poisson, the package default at a boundary null (dispersion_test()), when the
  ## mean model is within the comparator gates; otherwise it is the asymptotic
  ## mixture above, which is conservative for this test. The bootstrap restores the
  ## random-number state, so every other result is the same with or without it.
  p_nb_method <- "asymptotic"; nb_B_ok <- NA_integer_
  if (is.null(nb_boot)) nb_boot <- identical(ztp_threshold, "bootstrap")
  if (isTRUE(nb_boot) && is.finite(lr_nb) && k <= comp_max_par && n <= comp_max_n) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
    mu0 <- fitted(pois)
    LRb <- tryCatch(unlist(.ud_lapply(seq_len(ztp_boot_B), function(b) {
      yb <- stats::rpois(length(mu0), mu0)
      pf <- tryCatch(stats::glm.fit(X, yb, family = stats::poisson()), error = function(e) NULL)
      nf <- tryCatch(suppressWarnings(MASS::glm.nb(yb ~ 0 + X)), error = function(e) NULL)
      if (is.null(pf) || is.null(nf)) return(NA_real_)
      max(2 * (as.numeric(stats::logLik(nf)) - sum(stats::dpois(yb, pf$fitted.values, log = TRUE))), 0)
    }, cores), use.names = FALSE),
    finally = {
      if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
      else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
    })
    nb_B_ok <- sum(is.finite(LRb))
    if (nb_B_ok >= 0.8 * ztp_boot_B) {
      p_nb <- (1 + sum(LRb[is.finite(LRb)] >= lr_nb - 1e-10)) / (1 + nb_B_ok)
      p_nb_method <- "parametric bootstrap"
    } else {
      warning("ud_screen: the NB-vs-Poisson bootstrap failed on ", ztp_boot_B - nb_B_ok, " of ", ztp_boot_B,
              " replicates; the p-value stays asymptotic.")
    }
  }
  verdict_marg <- .ud_verdict_marg(dd, p_nb)
  idx_uncond <- var(yv) / mean(yv); pct_zero <- mean(yv == 0)

  ## at-risk (positive counts) vs zero-truncated Poisson, calibrated threshold
  pos <- yv > 0; n_pos <- sum(pos)
  pearson_ztp <- NA_real_; verdict_pos <- NA_character_; ztp_ll <- NA_real_
  thr <- .ud_ztp_threshold(n_pos)
  thr_hi <- 1.15   # conventional over cutoff; bootstrap mode replaces it with the null's 95th percentile
  ## over-conditioning guard state: `atrisk_skipped` = the mean model (nearly)
  ## saturates the positives, so the at-risk statistic is not computed at all;
  ## `overconditioned` = the fit ran but the parameter share of the positives,
  ## k_fit / n_pos, is at or above 0.10 -- the region where the calibrated
  ## threshold's size inflates with p/n_+ (the companion paper's App. B drift
  ## result), so the verdict is flagged and the bootstrap threshold recommended.
  atrisk_skipped <- FALSE; overconditioned <- NA; sat_ratio <- NA_real_
  if (n_pos <= (k + 2) || n_pos < 10L) {           # too few positives for an at-risk verdict
    atrisk_skipped <- TRUE; overconditioned <- TRUE
    sat_ratio <- k / max(n_pos, 1L)
  }
  if (!atrisk_skipped) {
    used <- rownames(mf); data_pos <- data[used[pos], , drop = FALSE]
    ztp <- tryCatch(suppressWarnings(VGAM::vglm(formula, family = VGAM::pospoisson, data = data_pos)),
                    error = function(e) NULL)
    if (is.null(ztp)) {
      ## VGAM's IRLS can fail internally on dummy-heavy fits from its default
      ## start (an if(NA) stop in its step-halving). A positives-only Poisson
      ## GLM always converges and sits close to the ZTP optimum, so retry once
      ## from its coefficients; converged first-try fits never reach this path.
      gs <- tryCatch(stats::glm(formula, family = stats::poisson(), data = data_pos),
                     error = function(e) NULL)
      if (!is.null(gs)) {
        cs <- gs$coefficients; cs[!is.finite(cs)] <- 0
        ztp <- tryCatch(suppressWarnings(VGAM::vglm(formula, family = VGAM::pospoisson,
                                                    data = data_pos, coefstart = cs)),
                        error = function(e) NULL)
      }
    }
    ## A fit is usable only if its rates are sane: a diverged IRLS (VGAM's
    ## half-step collapse, or a silent divergence that raises no error) can
    ## leave a unit's fitted rate orders of magnitude above any observed count
    ## and corrupt the Pearson statistic. Fitted rates are the decisive
    ## criterion -- an iteration-cap or warning-based kill would also discard
    ## slow-converging but numerically sound fixed-effects screens. When the
    ## vglm fit is missing or unusable, fall through to our own ZTP Fisher
    ## scoring (score (y - m)x, weight m'(lambda)*lambda, step-halved on the
    ## log-likelihood): a converged Fisher fit is the same MLE, and it survives
    ## the dummy-heavy designs where VGAM's IRLS does not.
    ypos <- yv[pos]; rate_cap <- 20 * max(ypos)
    lam <- NULL; ll_fit <- NA_real_; k_ztp <- NA_integer_
    if (!is.null(ztp)) {
      lv <- exp(VGAM::predict(ztp)[, 1])
      if (max(lv) <= rate_cap) {
        lam <- lv; ll_fit <- as.numeric(VGAM::logLik(ztp))
        k_ztp <- length(stats::coef(ztp))
      }
    }
    if (is.null(lam)) {
      fz <- .ztp_fisher(formula, data_pos)
      if (!is.null(fz) && isTRUE(fz$converged) && max(fz$lam) <= rate_cap) {
        lam <- fz$lam; ll_fit <- fz$ll; k_ztp <- fz$rank
      }
    }
    if (!is.null(lam)) {
      ## The divisor's parameter count is the rank of the fit that produced the
      ## residuals: on zero-heavy panels the positives-only design drops the
      ## all-zero units' levels, so the full-data mean model's column count
      ## overstates it (and inflates the statistic).
      m    <- lam / (1 - exp(-lam)); v <- m * (1 + lam - m)
      pearson_ztp <- sum((ypos - m)^2 / v) / (n_pos - k_ztp)
      ztp_ll <- ll_fit
      sat_ratio <- k_ztp / n_pos
      overconditioned <- sat_ratio >= 0.10
      if (ztp_threshold == "bootstrap") {
        ## Parametric bootstrap of the ZTP null at the fitted rates: the exact
        ## finite-sample null of THIS design, so the cutoff carries the
        ## fixed-effects estimation drift that the calibrated rule does not.
        mfz <- stats::model.frame(formula, data_pos, drop.unused.levels = TRUE)
        Xz  <- stats::model.matrix(attr(mfz, "terms"), mfz)
        p0  <- exp(-lam)
        Tb <- unlist(.ud_lapply(seq_len(ztp_boot_B), function(b) {
          yb <- stats::qpois(p0 + stats::runif(n_pos) * (1 - p0), lam)
          fb <- .ztp_fisher_fit(Xz, yb, maxit = 500L)
          if (is.null(fb) || !isTRUE(fb$converged)) return(NA_real_)
          if (max(fb$lam) > 20 * max(yb)) return(NA_real_)
          mb <- fb$lam / (1 - exp(-fb$lam)); vb <- mb * (1 + fb$lam - mb)
          sum((yb - mb)^2 / vb) / (n_pos - fb$rank)
        }, cores), use.names = FALSE)
        n_ok <- sum(is.finite(Tb))
        if (n_ok >= 0.8 * ztp_boot_B) {
          thr    <- as.numeric(stats::quantile(Tb, 0.05, na.rm = TRUE))
          thr_hi <- as.numeric(stats::quantile(Tb, 0.95, na.rm = TRUE))
        } else {
          warning("ud_screen: parametric-bootstrap threshold failed on ",
                  ztp_boot_B - n_ok, " of ", ztp_boot_B,
                  " replicates; falling back to the calibrated threshold.")
          ztp_threshold <- "calibrated"
        }
      }
      verdict_pos <- if (pearson_ztp < thr) "UNDERDISPERSED" else
                     if (pearson_ztp > thr_hi) "OVERDISPERSED" else
                     "~equi (only as tight as truncation forces)"
    }
  }

  ## CPB (untruncated) only when the data are not zero-dominated
  cpb_ll <- NA_real_; ceil_exceed <- NA_real_; cpb_alpha <- NA_real_
  if (run_cpb && pct_zero < 0.20 && n <= cpb_max_n) {
    cf <- tryCatch(suppressWarnings(cpb(formula, data = data, truncated = FALSE, se = "none")),
                   error = function(e) NULL)
    if (!is.null(cf) && is.finite(cf$loglik) && !is.na(cf$alpha) && cf$alpha < 0.999) {
      cpb_ll <- cf$loglik; cpb_alpha <- cf$alpha
      ceil_exceed <- mean(yv > floor(cf$ceiling))
    }
  }
  gp_ll <- NA_real_
  if (run_gp && k <= comp_max_par && n <= comp_max_n) {
    ## the package's own generalized Poisson, whose dispersion may be negative
    ## (an overdispersion-only GP would report the Poisson fit on underdispersed
    ## data); gated like the COM-Poisson, since it is a full numerical-gradient fit
    gp <- tryCatch(suppressWarnings(count_reg(formula, data = data, family = "genpois", se = "none")),
                   error = function(e) NULL)
    if (!is.null(gp) && is.finite(gp$loglik)) gp_ll <- gp$loglik
  }
  ## native COM-Poisson comparator (soft tail); gated because it has no
  ## concentrated fixed-effects path and the normalizing constant is per-obs work
  comp_ll <- NA_real_; comp_nu <- NA_real_
  if (run_comp && k <= comp_max_par && n <= comp_max_n) {
    cmp <- tryCatch(suppressWarnings(count_reg(formula, data = data, family = "compois", se = "none")),
                    error = function(e) NULL)
    if (!is.null(cmp) && is.finite(cmp$loglik)) { comp_ll <- cmp$loglik; comp_nu <- cmp$theta }
  }

  structure(list(
    formula = formula, n = n, n_pos = n_pos,
    verdict_marginal = verdict_marg, verdict_atrisk = verdict_pos,
    idx_uncond = idx_uncond, pct_zero = pct_zero, mean = mean(yv), max = max(yv),
    dd = dd, pearson_ztp = pearson_ztp, ztp_threshold = thr, ztp_threshold_hi = thr_hi,
    ztp_threshold_method = ztp_threshold, ztp_ll = ztp_ll,
    lr_nb = lr_nb, p_nb = p_nb, p_nb_method = p_nb_method, nb_boot_B_ok = nb_B_ok,
    overconditioned = overconditioned, sat_ratio = sat_ratio,
    atrisk_skipped = atrisk_skipped,
    ll = c(Poisson = as.numeric(logLik(pois)),
           NB = if (!is.null(nb)) as.numeric(logLik(nb)) else NA_real_,
           GP = gp_ll, COMP = comp_ll, CPB = cpb_ll),
    comp_nu = comp_nu,
    cpb_alpha = cpb_alpha, ceil_exceed = ceil_exceed,
    run_cpb = run_cpb, cpb_eligible = (pct_zero < 0.20), digits = digits),
    class = "ud_screen")
}

#' @export
print.ud_screen <- function(x, ...) {
  d <- x$digits
  cat("\n=== Underdispersion screen ===\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat(sprintf("N=%d  mean=%.3f  max=%d  %%zero=%.1f  (unconditional var/mean=%.2f)\n",
              x$n, x$mean, x$max, 100 * x$pct_zero, x$idx_uncond))
  cat("\nMARGINAL verdict:", x$verdict_marginal, "\n")
  cat(sprintf("   Pearson=%.3f  prop.slope=%.3f (p=%s)\n", x$dd$pearson, x$dd$slope,
              format.pval(x$dd$p, digits = 2)))
  if (!is.na(x$p_nb)) {
    how <- if (identical(x$p_nb_method, "parametric bootstrap")) "by parametric bootstrap"
           else "asymptotic, conservative at this boundary"
    if (round(x$lr_nb, 2) > 0)
      cat("   NB vs Poisson LR =", round(x$lr_nb, 2), "(p=", format.pval(x$p_nb, digits = 2), how,
          "; sig => overdispersion)\n")
    else cat("   NB vs Poisson LR = 0, at the Poisson boundary (p=", format.pval(x$p_nb, digits = 2), how, ")\n")
  }
  if (isTRUE(x$atrisk_skipped))
    cat(sprintf("\nAT-RISK screen SKIPPED -- over-conditioning: the mean model (nearly)\nsaturates the positive counts (%d parameters vs n_pos = %d); any within-unit\ntightness at this saturation would be manufactured by the specification,\nnot measured.\n",
                as.integer(round(x$sat_ratio * max(x$n_pos, 1L))), x$n_pos))
  if (!is.na(x$pearson_ztp)) {
    cat("\nAT-RISK (y>0) verdict:", x$verdict_atrisk, "  [n_pos=", x$n_pos, "]\n", sep = "")
    cat(sprintf("   ZTP-Pearson = %.3f  (underdispersed if < %.3f, the %s 5%% threshold)\n",
                x$pearson_ztp, x$ztp_threshold,
                if (identical(x$ztp_threshold_method, "bootstrap")) "parametric-bootstrap"
                else "calibrated"))
    if (isTRUE(x$overconditioned) && identical(x$ztp_threshold_method, "calibrated"))
      cat(sprintf("   over-conditioning caution: the mean model absorbs %.0f%% of the positives'\n   df (share %.2f >= 0.10); the calibrated threshold is anti-conservative at this\n   share -- use ztp_threshold = \"bootstrap\" and read the verdict as diagnostic,\n   not probative.\n",
                  100 * x$sat_ratio, x$sat_ratio))
    else if (isTRUE(x$overconditioned))
      cat(sprintf("   over-conditioning note: parameter share of positives = %.2f (>= 0.10);\n   bootstrap threshold in use, verdict remains diagnostic rather than probative.\n",
                  x$sat_ratio))
  }
  cat("\nModel comparison (log-lik):\n"); print(round(x$ll, d))
  if (!is.null(x$comp_nu) && !is.na(x$comp_nu))
    cat(sprintf("\nCOM-Poisson (full data): nu = %.2f  (> 1 = underdispersed, soft tail)\n", x$comp_nu))
  if (!is.na(x$cpb_alpha))
    cat(sprintf("\nCPB (full data): alpha = %.3f; ceiling-exceedance = %.1f%% of obs\n",
                x$cpb_alpha, 100 * x$ceil_exceed))
  cat("\n")
  invisible(x)
}
