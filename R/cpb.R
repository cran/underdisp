## ---------------------------------------------------------------------------
## Continuous Parameter Binomial (CPB) regression: estimator + internal helpers
## ---------------------------------------------------------------------------

## The levels of a model frame's factor and character columns, so that new data
## are coded as the fitting data were (a character column's levels are its
## sorted values, as in stats::.getXlevels(); without them a newdata column
## holding one value would become a one-level factor).
.cpb_xlevels <- function(mf) {
  fc <- vapply(mf, function(v) is.factor(v) || is.character(v), logical(1))
  if (!any(fc)) return(NULL)
  lapply(mf[fc], function(v) if (is.factor(v)) levels(v) else sort(unique(v)))
}

## Negative log-likelihood via the compiled C++ backend (see src/cpb_ll.cpp).
## `offset` is a length-n vector added to the log-mean (an exposure); 0 gives the
## plain model. With frequency weights `w` the weighted form is used.
.cpb_nll <- function(params, X, Y, offset, max.support, truncated, w = NULL, support) {
  if (is.null(w)) cpb_nll_cpp(params, X, as.integer(Y), offset, as.integer(max.support), truncated, .cpb_fl(support))
  else cpb_wnll_cpp(params, X, as.integer(Y), w, offset, as.integer(max.support), truncated, .cpb_fl(support))
}

## Single Nelder-Mead fit from one starting value on the (column-scaled) design.
.cpb_fit_one <- function(Xs, Y, start, offset, max.support, truncated, maxit, reltol, w = NULL, support) {
  optim(start, .cpb_nll, X = Xs, Y = Y, offset = offset, max.support = max.support, truncated = truncated,
        w = w, support = support, method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol))
}

## The largest intercept shift any count needs to be feasible at dispersion a
## and linear predictor eta (-Inf when none binds): the floor support needs
## y <= n, King's y - 1 < n, so on King's support a count of one is always
## feasible. A start shifted by this plus a margin is never on the infeasible
## plateau.
.cpb_feas_shift <- function(Y, eta, a, support) {
  if (.cpb_fl(support)) {
    pos <- Y > 0
    if (!any(pos)) return(-Inf)
    max(log(Y[pos] * (1 - a)) - eta[pos])
  } else {
    pos <- Y > 1
    if (!any(pos)) return(-Inf)
    max(log((Y[pos] - 1) * (1 - a)) - eta[pos])
  }
}

## King's CPB approaches the Poisson continuously as alpha -> 1, so a fit to
## data without underdispersion runs to that limit (the floor support's guard
## warning covers the same case there).
## a CPB fit at that limit (King's support only; the floor support's guard is
## reported by $support_binding)
.cpb_at_limit <- function(object)
  inherits(object, c("cpb", "cpb_fe")) && !.cpb_fl(.cpb_support(object)) && isTRUE(object$alpha > 0.999)
.cpb_warn_limit <- function(alpha, support) {
  if (!.cpb_fl(support) && alpha > 0.999)
    warning("alpha-hat (", signif(alpha, 6), ") is at the CPB's equidispersed limit: the fit is the Poisson, ",
            "and the data show no underdispersion (gec() lets the dispersion go either way).", call. = FALSE)
  invisible(NULL)
}
## The other end: as alpha -> 0 the pmf at each rate concentrates on at most two
## adjacent counts (under King's support the top two points, split by where the
## ceiling sits above an integer), which a series taking only two neighbouring
## values fits exactly. Such a maximum is legitimate, so the print methods note
## it rather than warn; the interval and tests for alpha are then at the boundary.
.cpb_lower_note <- function(alpha) {
  if (is.finite(alpha) && alpha < 1e-3)
    cat("Note: alpha-hat is at the CPB's lower limit (alpha -> 0): the fitted distribution puts its mass on at most\n",
        "two adjacent counts at each rate, and the interval and tests for alpha are at that boundary.\n", sep = "")
  invisible(NULL)
}

## Multi-start fit. On the floor support the surface is piecewise smooth: every
## time a parameter move shifts an observation's support floor, the normalizing
## constant jumps, so a simplex can settle on a local tooth; on King's support it
## is continuous with kinks. The fit multi-starts over the dispersion parameter,
## polishes the incumbent by Nelder-Mead restarts, and then restarts from
## deterministic perturbations of the incumbent whose intercept is raised to
## keep every count feasible, keeping the best (on King's support a final BFGS
## pass follows). The design columns are scaled to unit standard deviation for
## the optimizer and the coefficients mapped back, so the fit is invariant to
## the covariates' units. Returns the lowest-objective fit with `par` on the
## original scale.
.cpb_fit_ms <- function(X, Y, offset, max.support, truncated, alpha.starts, maxit, reltol, w = NULL,
                        restarts = 2L, start = NULL, support) {
  wv <- .ud_w1(w, nrow(X))
  s <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/"); p <- ncol(X)
  fn <- function(par) .cpb_nll(par, Xs, Y, offset, max.support, truncated, w, support = support)
  bs <- if (!is.null(start)) start * s                       # slopes given on the original scale
        else tryCatch({ v <- glm.fit(Xs, Y, weights = wv, offset = offset, family = poisson())$coefficients
                        v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  icol <- match("(Intercept)", colnames(X)); pos <- Y > 0
  ## feasibility repair: raise the intercept until every count is feasible,
  ## so that a start (and a restart) is never on the infeasible plateau
  repair <- function(st) {
    if (is.na(icol) || !any(pos)) return(st)
    a <- plogis(st[p + 1L]); eta <- offset + as.numeric(Xs %*% st[seq_len(p)])
    need <- .cpb_feas_shift(Y, eta, a, support) + 1e-6
    if (need > 0) st[icol] <- st[icol] + need
    st
  }
  cand <- lapply(alpha.starts, function(a0)
    tryCatch(.cpb_fit_one(Xs, Y, repair(c(bs, qlogis(a0))), offset, max.support, truncated, maxit, reltol, w,
                          support = support),
             error = function(e) NULL))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) return(NULL)
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  best <- .ud_nm_polish(best, fn, maxit, reltol, restarts = 1L)
  if (restarts > 0L && best$value < 1e9) {
    pats <- list(rep(0.15, p + 1L), rep(-0.15, p + 1L), 0.15 * (-1)^seq_len(p + 1L), -0.15 * (-1)^seq_len(p + 1L))
    for (dlt in pats[seq_len(min(restarts, 4L))]) {
      st <- repair(best$par + dlt)
      o2 <- tryCatch(.cpb_fit_one(Xs, Y, st, offset, max.support, truncated, maxit, reltol, w, support = support),
                     error = function(e) NULL)
      if (!is.null(o2) && o2$value < best$value - 1e-8) best <- .ud_nm_polish(o2, fn, maxit, reltol)
    }
  }
  if (!.cpb_fl(support) && best$value < 1e9) best <- .ud_bfgs_pass(best, fn)
  best$par[seq_len(p)] <- best$par[seq_len(p)] / s
  best
}

## Baseline log-likelihood: (zero-truncated) Poisson MLE, the correct H0: alpha = 1.
.cpb_baseline_loglik <- function(X, Y, offset, truncated, w = NULL) {
  wv <- .ud_w1(w, nrow(X)); s <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/")
  nll <- if (truncated)
    function(b) { l <- as.numeric(exp(offset + Xs %*% b))
      v <- sum(wv * (dpois(Y, l, log = TRUE) - log(1 - exp(-l)))); if (!is.finite(v)) 1e10 else -v }
  else
    function(b) { l <- as.numeric(exp(offset + Xs %*% b))
      v <- sum(wv * dpois(Y, l, log = TRUE)); if (!is.finite(v)) 1e10 else -v }
  bs <- tryCatch({ v <- glm.fit(Xs, Y, weights = wv, offset = offset, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  o <- tryCatch(optim(bs, nll, method = "BFGS"), error = function(e) NULL)
  if (is.null(o)) NA_real_ else -o$value
}

#' Fit a Continuous Parameter Binomial (CPB) regression
#'
#' Fits the underdispersed continuous parameter binomial model of King (1989): a
#' binomial count with a continuous number of trials, \eqn{\lambda_i/(1-\alpha)},
#' and success probability \eqn{1-\alpha}, so that the conditional variance is
#' close to a fraction of the conditional mean,
#' \eqn{\mathrm{Var}(Y\mid x) \approx \alpha\,\mathrm{E}(Y\mid x)} with
#' \eqn{0 < \alpha < 1}, and each observation has an endogenous ceiling,
#' \eqn{\lceil \lambda_i/(1-\alpha) \rceil}. A zero-truncated variant (the default) conditions on
#' \eqn{Y \ge 1}, appropriate when underdispersion lives among the positive counts of
#' an otherwise zero-inflated outcome.
#'
#' The rate is modeled log-linearly, \eqn{\lambda_i = \exp(x_i'\beta)}. The
#' probabilities follow the Katz recursion
#' \eqn{P(Y = k+1)/P(Y = k) = (\lambda_i - k(1-\alpha))/(\alpha(k+1))} on the
#' support \eqn{0, \ldots, \lceil \lambda_i/(1-\alpha) \rceil}, King's
#' normalization as Winkelmann, Signorino and King (1995) give it. The CPB is
#' therefore the generalized event count model of [gec()] at
#' \eqn{\delta = \alpha}, and the two functions return the same fit wherever the
#' GEC's estimate of \eqn{\delta} is below one, given the same `truncated`
#' (`cpb()` defaults to the zero-truncated model, `gec()` to the untruncated
#' one). Where \eqn{\lambda_i/(1-\alpha)}
#' is an integer the distribution is the binomial with that many trials, with mean
#' \eqn{\lambda_i} and variance \eqn{\alpha\lambda_i}; elsewhere its mean is at
#' least \eqn{\lambda_i}. In King's notation `alpha` is \eqn{\sigma^2}, the
#' variance-to-mean ratio at an integer ceiling; Winkelmann, Signorino and King
#' write the ceiling \eqn{\lambda_i/(1-\alpha)} as \eqn{\alpha} and \eqn{1-\alpha}
#' as \eqn{p}. Because the support moves with the parameters, the
#' log-likelihood is continuous with kinks wherever an observation's ceiling
#' crosses an integer, and it tends to the Poisson log-likelihood as
#' \eqn{\alpha \to 1}; a fit that reaches that limit (\eqn{\hat\alpha > 0.999})
#' warns, and [dispersion_test()] gives it the boundary value of the test. The
#' maximizer is found by a deterministic multistart, Nelder-Mead from the
#' Poisson solution at `alpha.start` and seven further dispersions, the best
#' polished by Nelder-Mead with feasibility-repaired restarts and a final BFGS
#' pass, on covariates scaled to
#' unit standard deviation (the coefficients are mapped back, so the fit does not
#' depend on the covariates' units), and the reported maximum is the maximum of
#' the fit's own profile in `alpha`: the profile is traced by continuation on both
#' sides of the estimate (the slopes re-maximized at each step from the
#' neighboring solution) until it has dropped four log-likelihood units, and a
#' trace point above the multistart's value restarts the fit from there. The same
#' trace gives the profile interval of [confint.cpb()]. Unit intercepts belong in
#' the fixed-effects estimator [cpb_fe()], which solves each of them exactly (a
#' pooled fit with many dummy columns can stop short of it).
#'
#' The likelihood is not regular at the crossings. Just above an integer the
#' top support point carries probability proportional to the ceiling's distance
#' \eqn{\epsilon} past it, so the information an observation carries about its
#' rate grows like \eqn{1/\epsilon} (about \eqn{0.28/\epsilon} for
#' \eqn{\log\lambda} at \eqn{\alpha = 0.5} and a ceiling of \eqn{3 + \epsilon},
#' about \eqn{1.8/\epsilon} at \eqn{\alpha = 0.3}; King 1989, footnote 7, notes
#' that asymptotic normality needs non-integer ceilings). With continuous
#' covariates the ceilings' fractional parts are spread and the expected
#' information per observation is unbounded, so standard errors from the
#' numerical Hessian are unreliable (at small `alpha` they fall well below the
#' sampling variability of the estimates), and inference uses a
#' cold-multistart bootstrap for the coefficients and a profile-likelihood
#' interval for \eqn{\alpha}. In simulations from King's CPB at n = 300, the
#' bootstrap's 95% interval for a slope covers 0.91 at alpha 0.3, 0.95 at 0.5
#' and 0.94 at 0.7 (0.92 at alpha 0.3 with n = 1,000); on the floor support it
#' covers 0.95. The interval for
#' \eqn{\alpha} is first-order by default and calibrated by parametric bootstrap
#' after [calibrate_alpha()] (see [alpha_confint()]).
#'
#' `support = "floor"` gives the package's definition before version 0.2.0, on
#' \eqn{0, \ldots, \lfloor \lambda_i/(1-\alpha) \rfloor}, whose mean is at most
#' \eqn{\lambda_i}. Its log-likelihood is discontinuous: whenever a parameter move
#' carries an observation's ceiling across an integer its normalizing constant
#' jumps, by about \eqn{(1-\alpha)^{k}} at ceiling \eqn{k}. The fit then adds a
#' scan of the breakpoints in `alpha` next to the estimate; two parameterizations
#' of the same design (say, treatment and sum contrasts) can settle on different
#' teeth of the surface, with log-likelihoods differing by the order of the
#' jumps; and ceilings above `max.support` are infeasible, so `alpha` can be
#' bounded by that guard rather than by the data.
#'
#' Runtime: a fit with `se = "none"` takes about a second at 500 rows and four
#' at 2,000 on one core, the multistart and the profile trace included (about
#' four and fifteen on the floor support).
#' The bootstrap replicates are maximized by the multistart without the profile
#' trace (about 0.1 seconds each at 500 rows; on the floor support a
#' replicate can therefore sit on a different tooth from the point estimate, by
#' the order of the jumps); `cores` runs them in parallel.
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
#'   resamples observations, giving heteroskedasticity-robust errors (a
#'   Hessian-based sandwich is unreliable here; see Details). Robust errors
#'   protect inference about the best-fitting CPB's parameters. The mean effects
#'   computed from them are the data's only when the CPB is the data's process:
#'   on underdispersed counts from another process (binomial counts with
#'   varying trials, a COM-Poisson) the CPB's average marginal effects can be
#'   off by several percent at any sample size, while the robust Poisson's
#'   ([count_reg()] with `se = "robust"`) are consistent whenever the mean is
#'   log-linear; the two fits' [ame()] side by side is the check. Nor do robust errors protect
#'   the CPB's distributional claims (its ceiling, its predicted probabilities)
#'   when the CPB is wrong: King and Roberts (2015) read a gap between robust
#'   and model-based errors as a sign of misspecification to fix.
#'   [ceiling_exceedance()], [rootogram()], [pit_hist()] and [compare_models()]
#'   check those claims.
#' @param offset Optional offset on the log-mean scale (an exposure): a numeric
#'   vector or the name of a column in `data`, making the model a rate model.
#' @param weights Optional frequency weights: a numeric vector or the name of a
#'   column in `data`. A weight of \eqn{w_i} is equivalent to \eqn{w_i} copies of
#'   row \eqn{i}; the bootstrap resamples rows together with their weights, and
#'   `nobs()` returns the weight total.
#' @param cores Number of worker processes for the bootstrap (default 1). With
#'   `cores > 1` the replicates run on a socket cluster whose random streams are
#'   seeded from the calling session, so `set.seed()` before the call
#'   reproduces a parallel run (parallel and serial draws differ).
#' @param alpha.start Starting value for the dispersion parameter (default 0.5); the
#'   fit also multi-starts over a spread of alpha values.
#' @param max.support Guard on the evaluated support. On King's support the pmf
#'   is summed over the points within 36 log units of its mode, and a
#'   distribution with more such points than `max.support` is treated as
#'   infeasible (the guard caps the work for very large rates, binding at 500
#'   where \eqn{\alpha\lambda} exceeds about 880; it does not bind as `alpha`
#'   approaches 1). On the floor support parameter values whose
#'   implied ceiling exceeds it are infeasible (they lie in the `alpha` near 1
#'   region where that support explodes), and a fit whose largest ceiling
#'   reaches the guard is flagged with a warning and in `$support_binding`,
#'   because `alpha` is then bounded by the guard rather than by the data. The
#'   default `NULL` sets it to `max(500, 10 * max(y))`.
#' @param maxit,reltol Optimizer controls passed to [stats::optim()].
#' @param support The CPB's support. `"ceiling"` (the default) is King's (1989)
#'   CPB as Winkelmann, Signorino and King (1995) normalize it, on
#'   \eqn{0, \ldots, \lceil \lambda_i/(1-\alpha) \rceil}; `"floor"` is the
#'   package's definition before version 0.2.0, on
#'   \eqn{0, \ldots, \lfloor \lambda_i/(1-\alpha) \rfloor}, kept so that earlier
#'   results reproduce. The two coincide when the ceiling is an integer.
#'
#' @return An object of class `"cpb"`: a list with `coefficients`, `alpha`, bootstrap
#'   `se.beta`/`vcov`/`ci.beta`, `loglik`, `loglik.null`, `fitted.values`, `ceiling`
#'   (the observation-specific implied ceiling), the model frame pieces, and the
#'   original `call`.
#'
#' @references King, G. (1989). Variance specification in event count models:
#'   From restrictive assumptions to a generalized estimator.
#'   \emph{American Journal of Political Science}, 33(3), 762-784.
#'   Winkelmann, R., Signorino, C. S. and King, G. (1995). A correction for an
#'   underdispersed event count probability distribution. \emph{Political
#'   Analysis}, 5, 215-228.
#'   King, G. and Roberts, M. E. (2015). How robust standard errors expose
#'   methodological problems they do not fix, and what to do about it.
#'   \emph{Political Analysis}, 23(2), 159-179.
#'
#' @examples
#' set.seed(1)
#' n <- 200; x <- rnorm(n)
#' N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
#' y <- rbinom(n, N, 0.5)                 # underdispersed (var/mean approx 0.5)
#' d <- data.frame(y = y, x = x)
#' fit <- cpb(y ~ x, data = d[d$y > 0, ], se = "none")
#' summary(fit)
#' @seealso [ud_screen()], [confint.cpb()], [predict.cpb()], [implied_ceiling()],
#'   [dispersion_test()]
#' @export
cpb <- function(formula, data, truncated = TRUE, se = c("none", "bootstrap"),
                B = 500, cluster = NULL, offset = NULL, weights = NULL, cores = 1L,
                alpha.start = 0.5, max.support = NULL, maxit = 20000, reltol = 1e-8,
                support = c("ceiling", "floor")) {
  .ud_no_formula_offset(formula)
  se <- match.arg(se); support <- match.arg(support)
  cl <- match.call()
  .ud_check_whole(max.support, "max.support", null_ok = TRUE)
  if (!is.numeric(alpha.start) || !length(alpha.start) || any(!is.finite(alpha.start)) ||
      any(alpha.start <= 0 | alpha.start >= 1))
    stop("'alpha.start' must be numeric, with every value strictly between 0 and 1.", call. = FALSE)
  if (se == "bootstrap") .ud_check_whole(B, "B", lower = 2)
  mf <- model.frame(formula, data, na.action = na.omit, drop.unused.levels = TRUE)
  Y  <- model.response(mf); X <- model.matrix(formula, mf)
  n  <- length(Y); p <- ncol(X); rows <- .ud_kept_rows(mf, data)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (n == 0L) stop("No complete cases on the model variables.")
  off <- rep_len(0, n)                                          # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    if (is.null(ov)) stop("'offset' must name a column of 'data' or be a numeric vector.")
    if (length(ov) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    off <- as.numeric(ov[rows])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  w <- .ud_weights(weights, data, rows)

  ## optional clustering variable, aligned to the (possibly na.omit-reduced) model frame
  clid <- NULL
  if (!is.null(cluster)) {
    if (se != "bootstrap") {
      warning("'cluster' only affects the bootstrap; it is ignored with se = \"none\" (use se = \"bootstrap\").")
    } else {
      cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
      if (is.null(cv)) stop("'cluster' must be a column name in 'data' or a vector.")
      if (length(cv) != nrow(data)) stop("'cluster' must have one value per row of 'data'.")
      clid <- cv[rows]
      if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
    }
  }
  if (any(Y < 0) || any(Y != floor(Y)))
    stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1))
    stop("truncated = TRUE requires all Y >= 1; use truncated = FALSE for zero-inclusive data.")
  .ud_warn_all_zero(Y)
  .ud_rank_check(X)
  if (is.null(max.support)) max.support <- max(500L, 10L * max(Y))

  astarts <- unique(c(alpha.start, 0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95))   # the floor support's teeth in alpha are thin
  fit <- .cpb_fit_ms(X, Y, off, max.support, truncated, astarts, maxit, reltol, w = w, support = support)
  if (is.null(fit)) stop("All optimization starts failed.")
  if (fit$value >= 1e9)
    stop("no feasible fit: every start left some count above its implied ceiling (raise max.support, ",
         "currently ", max.support, ", or check the data).")
  ## The reported maximum is the maximum of the fit's own alpha profile: the
  ## profile is traced by continuation on both sides of alpha-hat (the slopes
  ## re-maximized at each alpha from the neighbouring solution), and when a
  ## point of the trace exceeds the multistart's value the fit is re-run from
  ## it with alpha free. The trace is stored for the profile interval. On the
  ## floor support a breakpoint scan in alpha first resolves the teeth next to
  ## the estimate; King's likelihood has no teeth.
  for (round in 1:2) {
    if (.cpb_fl(support)) fit <- .cpb_alpha_scan(X, Y, off, max.support, truncated, w, fit, maxit, reltol)
    tr <- .cpb_profile_trace(X, Y, off, max.support, truncated, w, fit, support = support)
    k <- which.max(tr$ll)
    if (tr$ll[k] <= -fit$value + 1e-8) break
    fit <- .cpb_fit_ms(X, Y, off, max.support, truncated, tr$alpha[k], maxit, reltol, w = w,
                       start = tr$b[[k]], support = support)
  }
  if (fit$convergence != 0)
    warning("Optimization did not converge (code ", fit$convergence, ").")

  beta  <- fit$par[1:p]; names(beta) <- colnames(X)
  alpha <- plogis(fit$par[p + 1]); loglik <- -fit$value
  profile <- tr[c("alpha", "ll", "b")]                      # the alpha-profile trace the fit is the maximum of
  lam   <- as.numeric(exp(off + X %*% beta))
  fitted <- .cpb_mean(lam, alpha, truncated, support)       # exact mean of the fitted pmf
  ## the floor support treats a ceiling above max.support as infeasible, so the
  ## guard can bound alpha; King's pmf does not form the ceiling, and its fit
  ## approaches the Poisson continuously as alpha -> 1
  binding <- .cpb_fl(support) && max(lam / (1 - alpha)) >= 0.95 * max.support
  if (binding)
    warning("the fitted ceiling reaches max.support (", max.support, "): alpha is bounded by the guard, not the data. ",
            "The data may not be underdispersed (compare gec(), which lets the dispersion go either way), ",
            "or raise max.support.")
  .cpb_warn_limit(alpha, support)
  loglik0 <- .cpb_baseline_loglik(X, Y, off, truncated, w = w)

  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL
  boot <- NULL; ci.beta <- NULL
  if (se == "bootstrap") {                  # rows, or whole clusters; cold multistart in each replicate
    boot <- .ud_boot(n, clid, function(idx, blk) {
      fb <- .cpb_fit_ms(X[idx, , drop = FALSE], Y[idx], off[idx], max.support, truncated,
                        c(0.2, 0.5, 0.8), maxit, 1e-6, w = if (is.null(w)) NULL else w[idx], restarts = 0L,
                        support = support)
      if (!is.null(fb) && fb$convergence == 0 && fb$value < 1e9) c(fb$par[1:p], plogis(fb$par[p + 1])) else NULL
    }, c(colnames(X), "alpha"), B, cores, label = if (is.null(clid)) "obs" else .ud_cluster_label(cluster))
    bs <- .ud_boot_stats(boot, colnames(X), beta)
    se.beta <- bs$se; vcv <- bs$vcov; ci.beta <- bs$ci
  }

  structure(list(
    coefficients = beta, alpha = alpha, se.beta = se.beta, vcov = vcv,
    ci.beta = ci.beta, boot = boot, loglik = loglik, loglik.null = loglik0,
    fitted.values = fitted, rate = lam, linear.predictors = as.numeric(off + X %*% beta), offset = off, offset_name = .ud_offset_name(offset),
    weights = w, nobs_weighted = if (is.null(w)) n else sum(w),
    ceiling = lam / (1 - alpha), residuals = Y - fitted,
    n = n, df = p + 1, p = p, truncated = truncated, support = support,
    max.support = max.support, support_binding = binding,
    profile = profile,
    se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = (fit$convergence == 0),
    formula = formula, terms = attr(mf, "terms"), levels = .cpb_xlevels(mf),
    contrasts = attr(X, "contrasts"), X = X, Y = Y, rows = rows, n_dropped = nrow(data) - n, call = cl),
    class = "cpb")
}

## Profile-likelihood interval for alpha: invert the LR test by fixing alpha and
## re-optimizing beta (on the column-scaled design). The interval is first-order
## and model-based. On King's support the probabilities are continuous in the
## parameters and the chi-square cut covers 0.940 to 0.956 in simulations from
## the CPB. On the floor support, floor(lambda_i/(1-alpha)) moves with the
## parameters, so the log-ceiling coefficients (the intercept less log(1-alpha),
## and the slopes) are endpoint-type parameters: the likelihood gains whenever an
## unobserved top support point is shed, the maximum stops at the first count on
## its ceiling, and alpha-hat inherits a bias toward zero from them (given the
## ceilings, alpha is a regular exponential-family parameter). With the
## chi-square cut the interval covers about 0.90 to 0.94 there, most misses on
## the upper side.
## The profile is traced by continuation: each alpha starts from the solution
## at the previous alpha, and the start is first repaired to feasibility (a
## smaller alpha shrinks every ceiling lambda_i/(1-alpha), so the previous
## coefficients can leave observations above their ceiling; raising the
## intercept restores feasibility) and then re-maximized over beta with
## restarts. Each limit is the OUTERMOST crossing of its cut on the trace (the
## profile is a saw-tooth, so the set above a cut need not be connected),
## located by uniroot between the last trace point above the cut and the next
## one out. The interval is one-sided (`boundary`) only when the profile has not
## dropped to the cut by the trace's lower end (alpha = 1e-6 on King's support, 0.005 on the floor).
## The profile machinery of the pooled fit. `.cpb_prof_at()` maximizes the
## slopes at a fixed alpha with the fit's own effort (BFGS, Nelder-Mead polish,
## feasibility-repaired perturbation restarts, since a single simplex can stop
## on a lower tooth); `.cpb_profile_trace()` traces the profile by continuation
## on both sides of alpha-hat, in steps of max(0.01, alpha-hat/20), until it has
## dropped `drop` log-likelihood units below its running maximum (or reached
## the parameter bounds); `.cpb_alpha_profile_ci()` reads the crossing of the
## likelihood-ratio cut off the trace and refines it by root finding.
.cpb_prof_at <- function(a, b0, Xs, Y, off, ms, tr, w, s, icol, pos, b_cold = NULL, support) {
  repair <- function(b) {                        # shift the intercept so every count is feasible
    if (is.na(icol) || !any(pos)) return(b)
    eta <- off + as.numeric(Xs %*% b)
    need <- .cpb_feas_shift(Y, eta, a, support) + 1e-6
    if (need > 0) b[icol] <- b[icol] + need * s[icol]
    b
  }
  fn <- function(bb) .cpb_nll(c(bb, qlogis(a)), Xs, Y, off, ms, tr, w, support = support)
  ## from the neighbouring solution (continuation) and, when given, from the
  ## Poisson slopes (a cold start, as the multistart's alpha.start values use)
  run <- function(st) {
    st <- repair(st)
    o <- tryCatch(optim(st, fn, method = "BFGS", control = list(maxit = 3000, reltol = 1e-8)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value) || o$value >= 1e9)
      o <- optim(st, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-8))
    o
  }
  o <- run(b0)
  if (!is.null(b_cold)) { oc <- run(b_cold); if (oc$value < o$value) o <- oc }
  o <- .ud_nm_polish(o, fn, 3000, 1e-8, restarts = 1L)
  p <- length(o$par)
  for (dlt in if (p > 1L) list(rep(0.15, p), -0.15 * (-1)^seq_len(p)) else list()) {
    st <- repair(o$par + dlt)
    o2 <- tryCatch(optim(st, fn, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-8)),
                   error = function(e) NULL)
    if (!is.null(o2) && o2$value < o$value - 1e-8) o <- .ud_nm_polish(o2, fn, 3000, 1e-8, restarts = 1L)
  }
  list(ll = -o$value, b = o$par)
}
## An exact one-dimensional search in alpha at the fitted slopes. At fixed
## beta the log-likelihood is a saw-tooth in alpha: observation i's ceiling
## reaches k at alpha = 1 - lambda_i/k, where its normalizer grows and the log
## pmf drops by about (1-alpha)^k, so the supremum of a tooth is its interior
## critical point or its right end (the left limit of the next breakpoint).
## Every breakpoint within `window` of the estimate whose jump exceeds 1e-7 is
## enumerated, the objective evaluated at its left limit, the best tooth
## refined by optimize(), and, when the winner beats the incumbent, the fit is
## re-run from it with alpha free (up to three rounds). The profile trace sees
## the large-scale shape of the profile; this step resolves the teeth next to
## the estimate, which are narrower than any trace step.
.cpb_alpha_scan <- function(X, Y, off, ms, tr, w, fit, maxit, reltol, window = 0.02) {   # floor support only
  p <- ncol(X); wv <- .ud_w1(w, length(Y))
  s <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/")
  fn <- function(par) .cpb_nll(par, Xs, Y, off, ms, tr, w, support = "floor")
  par <- fit$par; par[seq_len(p)] <- par[seq_len(p)] * s
  for (round in 1:3) {
    ah <- plogis(par[p + 1L]); lam <- exp(off + as.numeric(Xs %*% par[seq_len(p)]))
    lo <- max(0.005, ah - window); hi <- min(0.995, ah + window)
    kcap <- floor(log(1e-7) / log1p(-hi))
    bp <- unlist(lapply(seq_along(Y), function(i) {
      k1 <- max(Y[i], ceiling(lam[i] / (1 - lo))); k2 <- min(floor(lam[i] / (1 - hi)), kcap)
      if (k2 < k1 || k1 < 1) numeric(0) else 1 - lam[i] / seq(k1, k2)
    }))
    bp <- sort(unique(round(bp[bp > lo & bp < hi], 12)))
    if (!length(bp)) break
    cand <- bp - 1e-9                                        # left limits
    vals <- vapply(cand, function(a) fn(c(par[seq_len(p)], qlogis(a))), numeric(1))
    j <- which.min(vals); best <- list(par = c(par[seq_len(p)], qlogis(cand[j])), value = vals[j])
    if (j > 1L) {                                            # the interior of the winning tooth
      o <- stats::optimize(function(a) fn(c(par[seq_len(p)], qlogis(a))), c(bp[j - 1L] + 1e-9, cand[j]), tol = 1e-10)
      if (o$objective < best$value) best <- list(par = c(par[seq_len(p)], qlogis(o$minimum)), value = o$objective)
    }
    if (best$value >= fit$value - 1e-8) break
    o2 <- tryCatch(.cpb_fit_one(Xs, Y, best$par, off, ms, tr, maxit, reltol, w, support = "floor"),
                   error = function(e) NULL)
    if (!is.null(o2) && o2$value < best$value) best <- o2
    best <- .ud_nm_polish(best, fn, maxit, reltol, restarts = 1L)
    if (best$value >= fit$value - 1e-8) break
    if (is.null(best$convergence)) best$convergence <- 0L      # a scan point carries no optim() fields
    if (is.null(best$counts)) best$counts <- c(0L, 0L)
    fit <- best; fit$par[seq_len(p)] <- fit$par[seq_len(p)] / s; par <- best$par
  }
  fit
}
.cpb_profile_trace <- function(X, Y, off, ms, tr, w, fit, drop = 4, support) {
  p <- ncol(X); wv <- .ud_w1(w, length(Y))
  s <- .ud_colscale(X, wv); Xs <- sweep(X, 2, s, "/")
  icol <- match("(Intercept)", colnames(X)); pos <- Y > 0
  bs <- tryCatch({ v <- glm.fit(Xs, Y, weights = wv, offset = off, family = poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  ah <- plogis(fit$par[p + 1L]); bh <- fit$par[seq_len(p)] * s; llmax <- -fit$value
  step <- max(0.01, ah / 20)
  ## King's CPB has the two-point limit as alpha -> 0, so its trace continues below 0.01 in geometric steps down to
  ## 1e-6; the floor support stops at 0.005, as before version 0.2.0
  fl <- .cpb_fl(support); lb <- .cpb_alpha_floor(support)
  alpha <- ah; ll <- llmax; b <- list(bh / s)
  for (dir in c(-1, 1)) {
    a <- ah; bp <- bh; top <- llmax
    repeat {
      a <- if (!fl && dir < 0 && a - step < 0.01) a * 0.6 else a + dir * step
      if (a < lb || a > 0.995) break
      r <- .cpb_prof_at(a, bp, Xs, Y, off, ms, tr, w, s, icol, pos, b_cold = bs, support = support)
      alpha <- c(alpha, a); ll <- c(ll, r$ll); b <- c(b, list(r$b / s))
      if (r$ll > top) top <- r$ll
      if (r$ll < top - drop) break
      bp <- r$b
    }
  }
  ## a fine pass next to the estimate: the teeth of the profile in alpha can be
  ## narrower than the step, so the neighbourhood is sampled at a fifth of it
  for (a in ah + step * c(-0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8)) {
    if (a < lb || a > 0.995) next
    r <- .cpb_prof_at(a, bh, Xs, Y, off, ms, tr, w, s, icol, pos, b_cold = bs, support = support)
    alpha <- c(alpha, a); ll <- c(ll, r$ll); b <- c(b, list(r$b / s))
  }
  o <- order(alpha)
  list(alpha = alpha[o], ll = ll[o], b = b[o])
}
## the smallest alpha the profile visits: 0.005 on the floor support (as before version 0.2.0), 1e-6 on King's
.cpb_alpha_floor <- function(support) if (.cpb_fl(support)) 0.005 else 1e-6
.cpb_alpha_profile_ci <- function(object, level = 0.95, drops = NULL) {
  X <- object$X; Y <- object$Y; ms <- object$max.support; tr <- object$truncated; sup <- .cpb_support(object)
  off <- if (!is.null(object$offset)) object$offset else rep_len(0, length(Y))
  w <- object$weights
  s <- .ud_colscale(X, .ud_w1(w, length(Y))); Xs <- sweep(X, 2, s, "/")
  icol <- match("(Intercept)", colnames(X)); pos <- Y > 0
  llmax <- object$loglik
  ## `drops`: how far below the maximum each limit sits, c(lower, upper) in
  ## log-likelihood units; the chi-square cut by default, side-specific values
  ## when the signed root has been calibrated
  if (is.null(drops)) drops <- rep(qchisq(level, 1) / 2, 2L)
  drops <- pmax(as.numeric(drops), 0); cut <- llmax - drops
  step <- max(0.01, object$alpha / 20); fl <- .cpb_fl(sup); lb <- .cpb_alpha_floor(sup)
  ## a side of the trace is usable when it falls below its cut or ends at the
  ## parameter bound; the stored trace stops 4 units below the maximum, so a
  ## deeper cut re-runs it (the fit is not touched)
  at_lb <- function(a) if (fl) a - step < lb else a * 0.6 < lb   # the trace's next step down would pass the bound
  reaches <- function(pr) {
    if (is.null(pr)) return(FALSE)
    lo <- pr$alpha < object$alpha - 1e-12; hi <- pr$alpha > object$alpha + 1e-12
    (any(pr$ll[lo] < cut[1L]) || (any(lo) && at_lb(min(pr$alpha[lo]))) || (!fl && object$alpha * 0.6 < lb)) &&
      (any(pr$ll[hi] < cut[2L]) || (any(hi) && max(pr$alpha[hi]) + step > 0.995))
  }
  pr <- object$profile
  if (is.null(pr) || max(pr$ll) > llmax + 1e-8 || !reaches(pr)) {        # an older object, a stale trace, or a cut below it
    fit <- list(par = c(object$coefficients, qlogis(object$alpha)), value = -llmax)
    pr <- .cpb_profile_trace(X, Y, off, ms, tr, w, fit, drop = max(drops, qchisq(level, 1) / 2) + 2, support = sup)
  }
  prof <- function(a, b0) .cpb_prof_at(a, b0, Xs, Y, off, ms, tr, w, s, icol, pos, support = sup)$ll
  k <- which(abs(pr$alpha - object$alpha) < 1e-12)[1L]; if (is.na(k)) k <- which.max(pr$ll)
  cross <- function(idx, cut) {                              # idx runs outward from alpha-hat
    if (!length(idx)) return(NA_real_)
    if (cut >= llmax - 1e-10) return(object$alpha)           # a zero drop: the limit is the estimate itself
    above <- which(pr$ll[idx] >= cut)
    m <- if (length(above)) max(above) else 0L               # the outermost trace point still above the cut
    if (m == length(idx)) return(NA_real_)                   # the trace ends above the cut, at a parameter bound
    i <- if (m == 0L) k else idx[m]; j <- idx[m + 1L]
    f <- function(aa) prof(aa, pr$b[[i]] * s) - cut
    ## the root to 1e-4 on the floor support (as before version 0.2.0), to a thousandth of the bracket's lower end
    ## on King's, whose limits can lie near zero
    tl <- if (fl) 1e-4 else max(1e-9, 1e-3 * min(pr$alpha[i], pr$alpha[j]))
    tryCatch(uniroot(f, sort(c(pr$alpha[i], pr$alpha[j])), tol = tl)$root, error = function(e) pr$alpha[j])
  }
  lo <- cross(if (k > 1L) rev(seq_len(k - 1L)) else integer(0), cut[1L])
  hi <- cross(if (k < length(pr$ll)) seq(k + 1L, length(pr$ll)) else integer(0), cut[2L])
  ## a side that never reaches its cut ends at the parameter bound: 0.005 or
  ## 0.999 on the floor support; on King's the lower side reaches alpha's limit,
  ## 0 (the two-point limit), and the upper side 0.999 or the estimate beyond it
  bl <- if (fl) 0.005 else 0
  bu <- if (.cpb_fl(sup)) 0.999 else max(0.999, object$alpha)
  c(lower = if (is.na(lo)) bl else unname(lo),
    upper = if (is.na(hi)) bu else unname(hi),
    boundary = as.numeric(is.na(lo)))
}
