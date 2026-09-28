## ---------------------------------------------------------------------------
## The dispersion profile: conditional variance-to-mean ratio against the
## fitted mean, empirically and as each fitted family implies it.
##
## Underdispersion has more than one mechanism, and the mechanisms leave
## different signatures in how the conditional dispersion moves with the mean:
## a hard ceiling (CPB) and the Katz family hold Var/Mean constant (alpha,
## delta); the COM-Poisson's ratio drifts toward 1/nu only as the mean grows;
## the gamma-count and double Poisson sit near 1/alpha and 1/theta but bend
## at small means; the negative binomial's ratio rises linearly in the mean.
## Binning the data by the fitted mean and plotting the empirical Pearson
## ratio against each family's implied curve shows which mechanism the data
## follow and where a family's variance function fails.
## ---------------------------------------------------------------------------

## per-observation untruncated mean, variance, and P(0) for a single-equation fit
.disp_moments <- function(fit) {
  ## the CPB and the Katz family: exact moments of the fitted pmf (the rate
  ## and alpha * rate, delta * rate are the moments only on an unbounded or
  ## integer-lattice support)
  if (inherits(fit, c("cpb", "cpb_fe"))) {
    lam <- .cpb_rate(fit); a <- fit$alpha
    p0 <- vapply(lam, function(l) .cpb_pmf1(l, a, 0L, FALSE)[1L], numeric(1))
    m <- .cpb_moments(lam, a, FALSE)
    list(mean = m$mean, var = m$var, p0 = p0, truncated = isTRUE(fit$truncated))
  } else if (inherits(fit, c("gec", "gec_fe"))) {
    lam <- if (!is.null(fit$rate)) fit$rate else exp(fit$linear.predictors)
    ms <- if (!is.null(fit$max.support)) fit$max.support else 500L
    km <- .gec_kmax_qr(lam, fit$delta, ms)
    P <- gec_pmf_cpp(lam, fit$delta, as.integer(km), as.integer(ms))
    k <- seq_len(ncol(P)) - 1
    m1 <- as.numeric(P %*% k); m2 <- as.numeric(P %*% (k * k))
    list(mean = m1, var = pmax(m2 - m1^2, 0), p0 = P[, 1], truncated = isTRUE(fit$truncated))
  } else if (inherits(fit, "count_reg")) {
    fam <- .count_fam(fit$family); mu <- fit$mu                # natural parameter
    list(mean = fam$meanfun(mu, fit$theta), var = fam$varfun(mu, fit$theta),
         p0 = fam$p0(mu, fit$theta), truncated = isTRUE(fit$truncated))
  } else stop("dispersion_profile() takes single-equation fits: cpb, cpb_fe, gec, gec_fe, or count_reg.")
}

## conditional (zero-truncated) mean and variance from the untruncated moments
.disp_condition <- function(m) {
  if (!m$truncated) return(list(mean = m$mean, var = m$var))
  q <- pmax(1 - m$p0, 1e-12)
  e1 <- m$mean / q; e2 <- (m$var + m$mean^2) / q
  list(mean = e1, var = e2 - e1^2)
}

#' Dispersion profile: variance-to-mean ratio against the fitted mean
#'
#' Bins the observations by the fitted mean of the first model, computes the
#' empirical conditional variance-to-mean ratio in each bin (the mean of the
#' Pearson contributions \eqn{(y_i - \hat\mu_i)^2/\hat\mu_i}), and sets it
#' against the ratio each fitted family implies at those means. The families
#' differ in how their dispersion moves with the mean: the CPB and the Katz
#' family hold it essentially constant (at `alpha` and `delta`; the implied
#' curves are the exact moments of the fitted pmf, so the renormalized CPB's
#' ratio dips a little below `alpha` where the ceiling is only two or three),
#' the COM-Poisson, gamma-count, and double Poisson approach a constant only as
#' the mean grows and rise toward the Poisson at small means, and the negative
#' binomial's ratio rises linearly. The profile therefore
#' shows which mechanism the data follow, where a family's variance function
#' fails, and whether the observed underdispersion is confined to a range of
#' means. Zero-truncated fits are profiled on their conditional moments. With
#' frequency weights (those of the first model) the bins hold equal weight
#' rather than equal numbers of rows, each row kept whole, every bin statistic
#' is a weighted mean, and `n` is the bin's weight total, so a weighted fit is
#' profiled as the data it stands for.
#'
#' @param object A fitted single-equation model (`cpb`, `cpb_fe`, `gec`,
#'   `gec_fe`, or `count_reg`); its fitted means define the bins and the
#'   empirical ratio.
#' @param ... Further fitted models on the same data, whose implied ratios are
#'   added as columns (name them for readable labels).
#' @param bins Number of equal-frequency bins by fitted mean (default 10).
#' @param plot Draw the profile (default `TRUE`); the table is returned
#'   invisibly either way.
#' @param ylim,main,xlab,ylab Plot settings.
#' @return A data frame of class `"dispersion_profile"` with one row per bin:
#'   `bin`, `n`, `mean_fitted`, `mean_y`, `ratio_empirical`, and one
#'   `ratio_<model>` column per fitted model.
#' @examples
#' set.seed(4); n <- 400; x <- rnorm(n)
#' d <- data.frame(y = rgammacount(n, exp(0.8 + 0.6 * x), alpha = 2.5), x = x)
#' m_gc <- count_reg(y ~ x, d, family = "gammacount")
#' m_nb <- count_reg(y ~ x, d, family = "negbin")
#' m_cp <- cpb(y ~ x, d, truncated = FALSE, se = "none")
#' dispersion_profile(gammacount = m_gc, negbin = m_nb, cpb = m_cp)
#' @seealso [dispersion_test()], [compare_models()], [rootogram()]
#' @export
dispersion_profile <- function(object, ..., bins = 10, plot = TRUE, ylim = NULL,
                               main = "Dispersion profile", xlab = "Fitted mean",
                               ylab = "Conditional Var / Mean") {
  ## fits may all arrive by name (then `object` is missing and they sit in `...`)
  fits <- c(if (!missing(object)) list(object), list(...))
  if (!length(fits)) stop("Supply at least one fitted model.")
  nm <- names(fits); if (is.null(nm)) nm <- rep("", length(fits))
  nm <- ifelse(nzchar(nm), nm, vapply(fits, function(f) if (inherits(f, "count_reg")) f$family else class(f)[1L], character(1)))
  nm <- make.unique(nm)
  y <- .obs_counts(fits[[1L]])
  m1 <- .disp_condition(.disp_moments(fits[[1L]]))
  if (any(vapply(fits[-1L], function(f) length(.obs_counts(f)) != length(y), logical(1))))
    stop("All models must be fit on the same observations.")
  mu <- m1$mean
  wt <- .ud_fit_w(fits[[1L]], length(y))           # frequency weights of the first model
  ord <- order(mu)
  bin <- integer(length(mu))
  bin[ord] <- if (is.null(wt)) ceiling(seq_along(mu) * bins / length(mu))
              else { cw <- cumsum(wt[ord]); ceiling(cw * bins / cw[length(cw)]) }   # equal weight, rows whole
  bin <- pmin(pmax(bin, 1L), bins)
  fb <- factor(bin, levels = seq_len(bins))
  bmean <- function(v) if (is.null(wt)) as.numeric(tapply(v, fb, mean))
                       else as.numeric(tapply(wt * v, fb, sum) / tapply(wt, fb, sum))
  emp <- (y - mu)^2 / pmax(mu, 1e-12)
  nb <- if (is.null(wt)) as.integer(tabulate(bin, bins))
        else { s <- as.numeric(tapply(wt, fb, sum)); s[is.na(s)] <- 0; s }
  out <- data.frame(bin = seq_len(bins), n = nb, mean_fitted = bmean(mu), mean_y = bmean(y),
                    ratio_empirical = bmean(emp))
  for (j in seq_along(fits)) {
    mj <- .disp_condition(.disp_moments(fits[[j]]))
    out[[paste0("ratio_", nm[j])]] <- bmean(mj$var / pmax(mj$mean, 1e-12))
  }
  out <- out[out$n > 0, , drop = FALSE]
  class(out) <- c("dispersion_profile", "data.frame")
  attr(out, "models") <- nm
  if (plot) plot(out, ylim = ylim, main = main, xlab = xlab, ylab = ylab)
  invisible(out)
}

#' @rdname dispersion_profile
#' @param x A `"dispersion_profile"` table.
#' @method plot dispersion_profile
#' @export
plot.dispersion_profile <- function(x, ylim = NULL, main = "Dispersion profile", xlab = "Fitted mean",
                                    ylab = "Conditional Var / Mean", ...) {
  nm <- attr(x, "models"); cols <- paste0("ratio_", nm)
  vals <- c(x$ratio_empirical, unlist(x[cols]), 1)
  if (is.null(ylim)) ylim <- range(vals[is.finite(vals)])
  pal <- c("firebrick", "steelblue", "darkgreen", "darkorange", "purple", "grey40")
  graphics::plot(x$mean_fitted, x$ratio_empirical, pch = 19, ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)
  graphics::abline(h = 1, col = "grey60", lty = 3)
  for (j in seq_along(cols))
    graphics::lines(x$mean_fitted, x[[cols[j]]], col = pal[(j - 1L) %% length(pal) + 1L], lwd = 2)
  graphics::legend("topright", bty = "n", legend = c("empirical", nm), pch = c(19, rep(NA, length(nm))),
                   lty = c(NA, rep(1, length(nm))), lwd = c(NA, rep(2, length(nm))),
                   col = c("black", pal[(seq_along(nm) - 1L) %% length(pal) + 1L]))
  invisible(x)
}

#' @export
print.dispersion_profile <- function(x, digits = 3, ...) {
  y <- as.data.frame(x); num <- vapply(y, is.numeric, logical(1)); y[num] <- lapply(y[num], round, digits)
  cat("Dispersion profile (", nrow(y), " bins by fitted mean; models: ", paste(attr(x, "models"), collapse = ", "), ")\n", sep = "")
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}
