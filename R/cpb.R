## ---------------------------------------------------------------------------
## Continuous Parameter Binomial (CPB) regression: estimator + internal helpers
## ---------------------------------------------------------------------------

## Factor levels of a model frame (for safe prediction on newdata), avoiding the
## non-exported stats:::.getXlevels.
.cpb_xlevels <- function(mf) {
  fac <- vapply(mf, is.factor, logical(1))
  if (!any(fac)) return(NULL)
  lapply(mf[fac], levels)
}

## Negative log-likelihood via the compiled C++ backend (see src/cpb_ll.cpp).
## `offset` is a length-n vector added to the log-mean (an exposure); 0 gives the
## plain model.
.cpb_nll <- function(params, X, Y, offset, max.support, truncated) {
  cpb_nll_cpp(params, X, as.integer(Y), offset, as.integer(max.support), truncated)
}

## Single Nelder-Mead fit from one starting value.
.cpb_fit_one <- function(X, Y, start, offset, max.support, truncated, maxit, reltol) {
  optim(start, .cpb_nll, X = X, Y = Y, offset = offset, max.support = max.support, truncated = truncated,
        method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
}

## Multi-start fit (the CPB surface is non-smooth at the feasibility boundary, so a
## single start can stick at a poor local optimum). Returns the lowest-objective fit.
.cpb_fit_ms <- function(X, Y, offset, max.support, truncated, alpha.starts, maxit, reltol) {
  bs <- tryCatch({ v <- glm.fit(X, Y, offset = offset, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  cand <- lapply(alpha.starts, function(a0)
    tryCatch(.cpb_fit_one(X, Y, c(bs, qlogis(a0)), offset, max.support, truncated, maxit, reltol),
             error = function(e) NULL))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) return(NULL)
  cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
}

## Baseline log-likelihood: (zero-truncated) Poisson MLE, the correct H0: alpha = 1.
.cpb_baseline_loglik <- function(X, Y, offset, truncated) {
  nll <- if (truncated)
    function(b) { l <- as.numeric(exp(offset + X %*% b))
      v <- sum(dpois(Y, l, log = TRUE) - log(1 - exp(-l))); if (!is.finite(v)) 1e10 else -v }
  else
    function(b) { l <- as.numeric(exp(offset + X %*% b))
      v <- sum(dpois(Y, l, log = TRUE)); if (!is.finite(v)) 1e10 else -v }
  bs <- tryCatch({ v <- glm.fit(X, Y, offset = offset, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  o <- tryCatch(optim(bs, nll, method = "BFGS"), error = function(e) NULL)
  if (is.null(o)) NA_real_ else -o$value
}

#' Fit a Continuous Parameter Binomial (CPB) regression
#'
#' Fits the underdispersed continuous parameter binomial model of King (1989), in
#' which the conditional variance is a fraction of the conditional mean,
#' \eqn{\mathrm{Var}(Y\mid x) = \alpha\,\mathrm{E}(Y\mid x)} with
#' \eqn{0 < \alpha < 1}, and each observation has an endogenous ceiling
#' \eqn{\lambda_i/(1-\alpha)}. A zero-truncated variant (the default) conditions on
#' \eqn{Y \ge 1}, appropriate when underdispersion lives among the positive counts of
#' an otherwise zero-inflated outcome.
#'
#' The mean is modelled log-linearly, \eqn{\lambda_i = \exp(x_i'\beta)}. Because the
#' support depends on the parameters, the log-likelihood is non-smooth at the
#' feasibility boundary and the numerical Hessian is unreliable; inference therefore
#' uses a cold-multistart bootstrap for the coefficients (validated to nominal
#' coverage) and a profile-likelihood interval for \eqn{\alpha} (see [confint.cpb()]).
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param truncated Logical; if `TRUE` (default) fit the zero-truncated CPB
#'   (requires all `Y >= 1`); if `FALSE` fit the untruncated CPB on `Y >= 0`.
#'   **Note the sibling default differs:** [cpb_fe()] defaults to
#'   `truncated = FALSE`, because its typical call fits a whole panel including
#'   zeros, while `cpb()`'s typical call fits the positive counts of a
#'   zero-inflated outcome. State `truncated` explicitly when moving a
#'   specification between the two.
#' @param se Inference method: `"none"` (default; fast, no standard errors) or
#'   `"bootstrap"` (cold-multistart pairs/cluster bootstrap; slower).
#' @param B Number of bootstrap resamples (default 500).
#' @param cluster Optional cluster identifier for cluster-robust inference: a
#'   column name in `data` or a vector aligned to its rows. When supplied, the
#'   bootstrap resamples whole clusters (a block bootstrap), giving
#'   cluster-robust standard errors and intervals; when `NULL` (default) it
#'   resamples observations, giving heteroskedasticity-robust errors. This is
#'   the appropriate route to robust inference here because the CPB's
#'   parameter-dependent support makes a Hessian-based sandwich unreliable.
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`, making the model a rate model.
#' @param alpha.start Starting value for the dispersion parameter (default 0.5); the
#'   fit also multi-starts over a spread of alpha values.
#' @param max.support Guard on the maximum evaluated support (default 500); fits with
#'   an implied ceiling above this are treated as infeasible (alpha near 1).
#' @param maxit,reltol Optimizer controls passed to [stats::optim()].
#'
#' @return An object of class `"cpb"`: a list with `coefficients`, `alpha`, bootstrap
#'   `se.beta`/`vcov`/`ci.beta`, `loglik`, `loglik.null`, `fitted.values`, `ceiling`
#'   (the observation-specific implied ceiling), the model frame pieces, and the
#'   original `call`.
#'
#' @references King, G. (1989). Variance specification in event count models.
#'   \emph{American Journal of Political Science}, 33(3), 762-784.
#'
#' @examples
#' set.seed(1)
#' n <- 300; x <- rnorm(n)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
#' y <- rbinom(n, N, 0.5)                 # underdispersed (var/mean approx 0.5)
#' d <- data.frame(y = y, x = x)
#' fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
#' summary(fit)
#' @seealso [ud_screen()], [confint.cpb()], [predict.cpb()], [implied_ceiling()]
#' @export
cpb <- function(formula, data, truncated = TRUE, se = c("none", "bootstrap"),
                B = 500, cluster = NULL, offset = NULL, alpha.start = 0.5, max.support = 500,
                maxit = 20000, reltol = 1e-8) {
  se <- match.arg(se)
  cl <- match.call()
  mf <- model.frame(formula, data, na.action = na.omit)
  Y  <- model.response(mf); X <- model.matrix(formula, mf)
  n  <- length(Y); p <- ncol(X)
  off <- rep_len(0, n)                                          # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    off <- as.numeric(ov[match(rownames(mf), rownames(data))])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }

  ## optional clustering variable, aligned to the (possibly na.omit-reduced) model frame
  clid <- NULL
  if (!is.null(cluster)) {
    cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    if (is.null(cv)) stop("'cluster' must be a column name in 'data' or a vector.")
    clid <- cv[match(rownames(mf), rownames(data))]
    if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
  }
  if (any(Y < 0) || any(Y != floor(Y)))
    stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1))
    stop("truncated = TRUE requires all Y >= 1; use truncated = FALSE for zero-inclusive data.")

  astarts <- unique(c(alpha.start, 0.15, 0.4, 0.65, 0.9))
  fit <- .cpb_fit_ms(X, Y, off, max.support, truncated, astarts, maxit, reltol)
  if (is.null(fit)) stop("All optimization starts failed.")
  if (fit$convergence != 0)
    warning("Optimization did not converge (code ", fit$convergence, ").")

  beta  <- fit$par[1:p]; names(beta) <- colnames(X)
  alpha <- plogis(fit$par[p + 1]); loglik <- -fit$value
  lam   <- as.numeric(exp(off + X %*% beta))
  loglik0 <- .cpb_baseline_loglik(X, Y, off, truncated)

  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL
  boot <- NULL; ci.beta <- NULL
  if (se == "bootstrap") {
    boot <- matrix(NA_real_, B, p + 1)
    grp <- if (!is.null(clid)) split(seq_len(n), clid) else NULL   # NULL = ordinary bootstrap
    for (b in seq_len(B)) {
      idx <- if (is.null(grp)) sample.int(n, n, replace = TRUE)     # resample rows, or
             else unlist(grp[sample.int(length(grp), length(grp), replace = TRUE)],
                         use.names = FALSE)                          # resample whole clusters
      fb <- .cpb_fit_ms(X[idx, , drop = FALSE], Y[idx], off[idx], max.support, truncated,
                        c(0.2, 0.5, 0.8), maxit, 1e-6)
      if (!is.null(fb) && fb$convergence == 0 && fb$value < 1e9) {
        bp <- fb$par; bp[p + 1] <- plogis(bp[p + 1]); boot[b, ] <- bp
      }
    }
    ok <- boot[complete.cases(boot), , drop = FALSE]
    if (nrow(ok) < 2) warning("Too few bootstrap resamples converged for stable inference.")
    se.beta <- setNames(apply(ok[, 1:p, drop = FALSE], 2, sd), colnames(X))
    vcv <- cov(ok[, 1:p, drop = FALSE]); dimnames(vcv) <- list(colnames(X), colnames(X))
    ci.beta <- t(apply(ok[, 1:p, drop = FALSE], 2, quantile, c(.025, .975)))
    rownames(ci.beta) <- colnames(X)
    attr(boot, "nboot_ok") <- nrow(ok)
  }

  structure(list(
    coefficients = beta, alpha = alpha, se.beta = se.beta, vcov = vcv,
    ci.beta = ci.beta, boot = boot, loglik = loglik, loglik.null = loglik0,
    fitted.values = lam, linear.predictors = as.numeric(off + X %*% beta), offset = off,
    ceiling = lam / (1 - alpha), residuals = Y - lam,
    n = n, df = p + 1, p = p, truncated = truncated, max.support = max.support,
    se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = (fit$convergence == 0),
    formula = formula, terms = attr(mf, "terms"), levels = .cpb_xlevels(mf),
    contrasts = attr(X, "contrasts"), X = X, Y = Y, call = cl),
    class = "cpb")
}

## Profile-likelihood interval for alpha: invert the LR test by fixing alpha and
## re-optimizing beta. Well calibrated except at strong underdispersion, where alpha
## sits at the feasibility boundary and the interval becomes one-sided (see paper).
.cpb_alpha_profile_ci <- function(object, level = 0.95) {
  X <- object$X; Y <- object$Y; ms <- object$max.support; tr <- object$truncated
  off <- if (!is.null(object$offset)) object$offset else rep_len(0, length(Y))
  llmax <- object$loglik; cut <- llmax - qchisq(level, 1) / 2
  bstart <- object$coefficients
  llprof <- function(a) {
    o <- optim(bstart, function(bb) cpb_nll_cpp(c(bb, qlogis(a)), X, as.integer(Y),
                                                off, as.integer(ms), tr),
               method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-8))
    -o$value
  }
  ah <- object$alpha
  lo <- tryCatch(uniroot(function(a) llprof(a) - cut, c(max(1e-3, ah / 6), ah - 1e-4))$root,
                 error = function(e) NA_real_)
  hi <- tryCatch(uniroot(function(a) llprof(a) - cut, c(ah + 1e-4, 0.999))$root,
                 error = function(e) NA_real_)
  boundary <- is.na(lo)
  c(lower = if (is.na(lo)) ah else unname(lo),
    upper = if (is.na(hi)) 0.999 else unname(hi),
    boundary = as.numeric(boundary))
}
