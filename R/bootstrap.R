## ---------------------------------------------------------------------------
## The package's one nonparametric bootstrap. Every estimator with
## se = "bootstrap" (cpb, cpb_fe, gec, gec_fe, zi_cpb, zi_gec, and the hurdles
## through their intensity fits) resamples rows, or whole clusters, B times
## through .ud_boot(), refits, and stores the draws in one layout:
##   - a matrix with one row per replicate and one named column per quantity:
##     the coefficients (prefixed "count:" and "zero:" for the zero-inflated
##     classes), the dispersion on its natural scale ("alpha" or "delta"), and,
##     for the fixed-effects fits, "fe_ref", the replicate's reference unit
##     effect (the mean of the unit effects the data identify);
##   - NA rows for replicates whose refit failed or did not report convergence;
##   - attributes `nboot_ok` (the converged replicates), `B`, and `cluster` (the
##     resampling unit: "obs", a column name, or "custom").
## The quantities of interest (first_difference(), irr()) and vcov(), confint()
## and summary() read the draws by column name.
## ---------------------------------------------------------------------------

## n rows, resampled as rows (cl = NULL) or as whole clusters (cl a vector of
## cluster ids on the rows). fit_one(idx, block) refits on the rows `idx` and
## returns the replicate's vector of length(cols), or NULL when the refit failed
## or did not converge; `block` gives, for each resampled row, the cluster draw
## it came from (so a fixed-effects refit can relabel a unit drawn twice as two
## units), and is NULL for a row bootstrap. A serial run times its first
## replicate and says how long the whole run will take when that is over a
## minute. A factor's unused levels are not clusters (split() would make them
## empty ones, drawn like the others).
.ud_boot <- function(n, cl, fit_one, cols, B, cores = 1L, label = "obs") {
  if (is.factor(cl)) cl <- droplevels(cl)
  grp <- if (is.null(cl)) NULL else split(seq_len(n), cl)
  if (!is.null(grp)) .ud_cluster_guard(cl)
  .ud_boot_size_note(B)
  one <- function(b) {
    if (is.null(grp)) { idx <- sample.int(n, n, replace = TRUE); blk <- NULL }
    else {
      gs <- sample.int(length(grp), length(grp), replace = TRUE)
      idx <- unlist(grp[gs], use.names = FALSE); blk <- rep(seq_along(gs), lengths(grp[gs]))
    }
    v <- tryCatch(fit_one(idx, blk), error = function(e) NULL)
    if (is.null(v) || length(v) != length(cols) || anyNA(v)) rep(NA_real_, length(cols)) else unname(as.numeric(v))
  }
  D <- matrix(unlist(.ud_run_reps(B, one, cores, "bootstrap")), nrow = B, byrow = TRUE, dimnames = list(NULL, cols))
  nok <- sum(stats::complete.cases(D))
  if (nok < 2L) warning("Too few bootstrap resamples converged for stable inference.", call. = FALSE)
  else if (nok < B / 2)
    warning(sprintf("only %d of %d bootstrap replicates converged; the standard errors rest on the survivors.", nok, B),
            call. = FALSE)
  attr(D, "nboot_ok") <- nok; attr(D, "B") <- as.integer(B); attr(D, "cluster") <- label
  D
}

## The replicates of a bootstrap, `one(b)` for b = 1..B, through .ud_lapply(). A serial run times its first
## replicate and says how long the whole run will take when that is over a minute (the same draws, in the same
## order, as running them all at once).
.ud_run_reps <- function(B, one, cores = 1L, what = "bootstrap") {
  serial <- is.na(as.integer(cores)[1L]) || as.integer(cores)[1L] <= 1L
  if (!serial || B < 2L) return(.ud_lapply(seq_len(B), one, cores))
  t0 <- proc.time()[["elapsed"]]; first <- one(1L); dt <- proc.time()[["elapsed"]] - t0
  if (dt * B > 60)
    message(sprintf("%s: %d replicates at about %s each, so about %s in all (the `cores` argument runs them in ",
                    what, B, .ud_duration(dt), .ud_duration(dt * B)), "parallel).")
  c(list(first), lapply(seq_len(B)[-1L], one))
}

## a bootstrap standard error has a Monte Carlo error of about 1/sqrt(2B) of
## itself, 10 percent at B = 50; below B = 100 the run says so
.ud_boot_size_note <- function(B) {
  if (B < 100)
    message(sprintf(paste0("B = %d bootstrap replicates: a bootstrap standard error carries a Monte Carlo error of ",
                           "about 1/sqrt(2B) of itself (%.0f%% here), and percentile intervals need more; use B of ",
                           "at least 200 for reported results."), B, 100 / sqrt(2 * B)))
  invisible(NULL)
}

.ud_duration <- function(s) {
  if (s < 90) sprintf("%.0f seconds", s) else if (s < 5400) sprintf("%.0f minutes", s / 60)
  else sprintf("%.1f hours", s / 3600)
}

## standard errors, covariance and 95% intervals of the columns `cols` of the
## stored draws: percentile intervals, or (normal = TRUE, the fixed-effects and
## mixture rule of confint.underdisp) normal ones centred on the estimate `est`
.ud_boot_stats <- function(D, cols, est, normal = FALSE) {
  ok <- D[stats::complete.cases(D), cols, drop = FALSE]
  if (nrow(ok) < 2L) return(list(se = stats::setNames(rep(NA_real_, length(cols)), names(est)), vcov = NULL, ci = NULL))
  se <- stats::setNames(apply(ok, 2, stats::sd), names(est))
  ci <- if (normal) cbind(lower = est - stats::qnorm(0.975) * se, upper = est + stats::qnorm(0.975) * se)
        else { q <- t(apply(ok, 2, stats::quantile, c(0.025, 0.975))); colnames(q) <- c("lower", "upper"); q }
  rownames(ci) <- names(est)
  v <- stats::cov(ok); dimnames(v) <- list(names(est), names(est))
  list(se = se, vcov = v, ci = ci)
}

## The converged draws of the columns `cols` (all columns when NULL), or NULL when the fit has fewer than two or
## does not store those columns. An object saved before version 0.2.0 stored its draws in the layout
## .ud_old_draw_names() gives (no dispersion column for the fixed-effects and zero-inflated CPB fits, and the
## zero-inflated GEC's on the log scale), with the column names rbind() took from the replicates' vectors (the
## coefficients' own, an empty one for the dispersion) and without the attributes .ud_boot() sets.
.ud_draws <- function(object, cols = NULL) {
  D <- object$boot
  if (is.null(D)) return(NULL)
  if (is.null(attr(D, "B"))) {
    nm <- .ud_old_draw_names(object)
    if (length(nm) != ncol(D)) return(NULL)
    colnames(D) <- nm
  }
  if (!is.null(cols)) {
    if (!all(cols %in% colnames(D))) return(NULL)
    D <- D[, cols, drop = FALSE]
  }
  D <- D[stats::complete.cases(D), , drop = FALSE]
  if (nrow(D) < 2L) NULL else D
}
## a bootstrapped fit whose draws have fewer than two converged replicates: its intervals are NA, and it says so
.ud_boot_too_few <- function()
  warning("fewer than two bootstrap replicates converged, so the intervals are NA.", call. = FALSE)
.ud_old_draw_names <- function(object) {
  b <- names(object$coefficients); z <- names(object$zero_coef)
  if (inherits(object, "zi_gec")) c(paste0("count:", b), ".log_delta", paste0("zero:", z))
  else if (inherits(object, "zi_cpb")) c(paste0("count:", b), paste0("zero:", z))
  else if (inherits(object, c("cpb_fe", "gec_fe"))) b
  else if (inherits(object, "gec")) c(b, "delta")
  else c(b, "alpha")
}

## the cluster label a fit records: a column name, or "custom" for a vector
.ud_cluster_label <- function(cluster)
  if (is.null(cluster)) NULL else if (is.character(cluster) && length(cluster) == 1L) cluster else "custom"
