## ---------------------------------------------------------------------------
## Split-panel (half-panel) jackknife for the concentrated fixed-effects
## estimators (Dhaene & Jochmans 2015), WITH the assumption checks the method
## requires.
##
## The concentrated FE-MLE's dispersion parameter carries an incidental-
## parameters bias of order 1/T. The correction refits the model on the first
## and second temporal halves of every unit's series and reports
##   theta_jack = 2 * theta_full - (theta_half1 + theta_half2) / 2,
## which removes the leading 1/T bias term. In the package's fixed-effects bias
## Monte Carlo (150 units, alpha 0.5):
## on King's support the alpha bias -0.122/-0.071/-0.036/-0.018 at
## T = 6/10/20/40 falls to +0.008/+0.004/+0.000/-0.001; on the floor support
## -0.153/-0.100/-0.057/-0.033 falls to -0.046/-0.024/-0.013/-0.009. Within each
## unit, rows are split in time order (the fits' `time` argument, which the
## correction requires).
##
## VALIDITY GATE. The jackknife identity requires the two half-panels to
## estimate the SAME pseudo-true parameter, differing only through the 1/T
## bias term (time-homogeneity; Dhaene & Jochmans 2015). On trending or
## time-heterogeneous panels the halves converge to different values and
## 2*full - mean(halves) removes the trend, not the bias.
##
## Two failure modes need two different detectors:
##   * ASYMMETRIC heterogeneity (regime change): the two temporal halves
##     estimate different dispersions. Detected by R2, PLACEBO-CALIBRATED:
##     five deterministic UNIT splits -- exchangeable under any time pattern --
##     measure the pure sampling noise of a half-sample fit difference (the
##     MEAN absolute placebo disagreement), and the time-split disagreement
##     must not exceed KAPPA times that noise.
##   * SYMMETRIC heterogeneity (a smooth trend): the halves agree with each
##     other but both differ systematically from the full fit, biasing the
##     full-vs-halves contrast the correction is built on, while R2 stays
##     blind. Detected by R2b, a residual time-trend check on the FULL fit:
##     a significant correlation between the Pearson residuals and the
##     within-unit time order means the mean model omits time structure; the
##     refusal message says to put time controls in the formula (which absorbs
##     the trend and makes the jackknife valid) rather than to trust the
##     correction.
## The gate REFUSES the correction (returning the ML fit with a warning) on
## R1 (corrected dispersion outside its feasible space), R2, or R2b. It WARNS
## (but proceeds) when the implied correction exceeds 0.15 * scale and when the
## median panel length is below 6. KAPPA = 2.5 and the R2b threshold (p < .01)
## are calibrated by Monte Carlo: stationary
## panels should rarely refuse (on King's support 9 to 13 percent do, on the
## floor support 5.5 to 9); regime-change panels must refuse (98.5 percent); and
## on trending panels the corrections that survive must not be worse than the
## uncorrected fit (every trending panel refuses). The placebo splits cost ten
## extra half-sample fits.
## ---------------------------------------------------------------------------

## A bootstrap replicate's refit runs inside .ud_replicate(). The validity gate is a check on the data, made once
## on the original fit, so a replicate applies the correction without it (its standard errors then carry the
## correction's variance). A replicate whose half fits fail fails (a missing replicate) instead of falling back to
## the uncorrected fit; one whose corrected dispersion leaves its space is held at the space's edge, since dropping
## it would censor one tail of the replicates. The floor support keeps its earlier replicates (cpb_fe() passes
## `replicate = FALSE` there).
.ud_replicate <- function(expr) {
  op <- options(underdisp.jackknife_replicate = TRUE); on.exit(options(op))
  expr
}
.fe_jack_replicate <- function() isTRUE(getOption("underdisp.jackknife_replicate"))

## refit(half_data) returns the c(coefficients, dispersion) vector of a half
## fit (dispersion LAST), or NULL on failure. `resid` and `torder` are the full
## fit's Pearson residuals and the (per-unit centered) time order, for R2b.
## Returns list(par = corrected vector, warn = character()) on success, or
## list(refused = "<reason>").
.fe_jackknife <- function(data_est, unit, refit, full_par, resid, torder,
                          disp_lower = 1e-3, disp_upper = Inf, kappa = 2.5, replicate = .fe_jack_replicate()) {
  iu  <- as.integer(factor(unit))
  ord <- stats::ave(seq_along(iu), iu, FUN = seq_along)
  len <- stats::ave(seq_along(iu), iu, FUN = length)
  if (replicate) {                            # a bootstrap replicate: the correction without the gate
    h <- ord <= len / 2; k <- length(full_par)
    q1 <- refit(data_est[h, , drop = FALSE]); q2 <- refit(data_est[!h, , drop = FALSE])
    if (is.null(q1) || is.null(q2)) return(list(refused = "a half-panel refit failed"))
    jk <- 2 * full_par - (q1 + q2) / 2
    eps <- 1e-6 * max(1, abs(disp_lower))
    jk[k] <- min(max(jk[k], disp_lower + eps), if (is.finite(disp_upper)) disp_upper - eps else Inf)
    return(list(par = jk, warn = character()))
  }
  warn <- character()
  if (stats::median(tabulate(iu)) < 6)
    warn <- c(warn, "median panel length is below 6; half-panel fits are short and the correction may be unreliable")
  ## R2b: residual time trend on the full fit (symmetric-heterogeneity check)
  tt <- torder - stats::ave(torder, iu)
  ct <- tryCatch(stats::cor.test(resid, tt), error = function(e) NULL)
  if (!is.null(ct) && is.finite(ct$p.value) && ct$p.value < 0.01)
    return(list(refused = sprintf(
      "the full fit's residuals trend with time (r = %+.3f, p = %.2g), so the mean model omits time structure and the split-panel jackknife's time-homogeneity requirement fails; include time controls (a trend or period effects) in the formula and re-fit",
      unname(ct$estimate), ct$p.value)))

  h1 <- ord <= len / 2
  p1 <- refit(data_est[h1, , drop = FALSE])
  p2 <- refit(data_est[!h1, , drop = FALSE])
  if (is.null(p1) || is.null(p2))
    return(list(refused = "a half-panel refit failed"))
  k  <- length(full_par)                      # dispersion is the last element
  d0 <- full_par[k]; d1 <- p1[k]; d2 <- p2[k]
  scale <- max(1, abs(d0))

  ## placebo: five deterministic unit splits (exchangeable under any time
  ## pattern) measure the sampling noise of a half-sample dispersion difference;
  ## the mean absolute disagreement is the noise estimate (five draws make it
  ## stable enough for a sharp threshold; a median of few draws is too crude)
  uid <- sort(unique(iu)); m <- seq_along(uid)
  splits <- list(uid[m %% 2L == 1L],
                 uid[m <= length(uid) / 2],
                 uid[m %% 4L %in% c(1L, 2L)],
                 uid[m %% 8L %in% c(1L, 2L, 3L, 4L)],
                 uid[m %% 8L %in% c(1L, 2L, 5L, 6L)])
  placebo <- vapply(splits, function(s) {
    q1 <- refit(data_est[iu %in% s, , drop = FALSE])
    q2 <- refit(data_est[!(iu %in% s), , drop = FALSE])
    if (is.null(q1) || is.null(q2)) NA_real_ else abs(q1[k] - q2[k])
  }, numeric(1))
  noise <- mean(placebo, na.rm = TRUE)
  if (!is.finite(noise) || sum(is.finite(placebo)) < 3L)
    return(list(refused = "the placebo unit-split fits failed, so the time-homogeneity check cannot be performed"))
  if (abs(d1 - d2) > max(0.05 * scale, kappa * noise))
    return(list(refused = sprintf(
      "the half-panel dispersion estimates disagree (%.3f vs %.3f) far beyond the exchangeable unit-split noise (mean %.3f); the panel appears time-heterogeneous and the split-panel jackknife is not valid here (Dhaene & Jochmans 2015)",
      d1, d2, noise)))

  jk <- 2 * full_par - (p1 + p2) / 2
  corr <- d0 - (d1 + d2) / 2                  # implied dispersion correction
  if (jk[k] <= disp_lower || jk[k] >= disp_upper)
    return(list(refused = sprintf(
      "the corrected dispersion (%.3f) leaves its feasible space; the two half-panels do not estimate a common parameter (the split-panel jackknife requires time-homogeneity; Dhaene & Jochmans 2015)",
      jk[k])))
  if (abs(corr) > 0.15 * scale)
    warn <- c(warn, sprintf(
      "the implied dispersion correction (%+.3f) is large relative to the validated 1/T magnitudes; verify the panel is time-homogeneous before relying on it", corr))
  list(par = jk, warn = warn)
}
