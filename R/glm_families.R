## ---------------------------------------------------------------------------
## Matched count families through one interface: single-equation, fixed-effects,
## zero-truncated, hurdle, and zero-inflated regressions for every family in
## R/count_families.R (Poisson, negative binomial, the classical and the
## mean-parameterized COM-Poisson, generalized Poisson, gamma-count, double
## Poisson), all fit and scored with identical df, log-likelihood, proper-score,
## weight, and robust/clustered-SE accounting, so compare_models() can place a
## cpb_fe next to a gamma-count FE fit or a zero-inflated generalized Poisson on
## exactly the same footing.
##
## Design: everything reduces to a per-observation log-likelihood in the natural
## parameter mu (log link) plus an optional shape parameter carried on the
## family's own unconstrained scale (fam$shape_link / fam$shape_inv). A single
## sandwich (.count_vcov) turns per-observation scores into analytic, HC-robust,
## or cluster-robust covariance, uniformly across every variant, with frequency
## weights entering the log-likelihood, the information, and the score products
## exactly as duplicated rows would.
## ---------------------------------------------------------------------------

## An optimizer excursion: a fitted mean outside [1e-12, 1e8]. The bound applies on the family's mean
## scale, which for the rate-parameterized COM-Poisson is lambda^(1/nu), its mean to leading order. On a
## nearly degenerate series the COM-Poisson likelihood rises along a ridge on which nu grows and lambda,
## about the mode to the power nu, grows without bound, so a bound on lambda itself would stop the fit
## short of the maximum.
.count_excursion <- function(eta, theta, fam) {
  le <- if (is.null(fam$log_mean)) eta else fam$log_mean(eta, theta)
  any(!is.finite(le)) || any(le > log(1e8)) || any(le < log(1e-12))
}

## per-observation log-likelihood in (beta, shape) for a count family, with
## optional zero-truncation. par = c(beta, shape_link(theta)) when the family
## has a shape.
.count_llik_i <- function(par, X, Y, fam, truncated, offset = 0) {
  p <- ncol(X)
  eta <- offset + as.numeric(X %*% par[seq_len(p)])       # offset enters the log-mean, exp(offset + x'b)
  theta <- if (fam$nshape) fam$shape_inv(par[p + seq_len(fam$nshape)]) else NA_real_
  if (.count_excursion(eta, theta, fam)) return(rep(-Inf, length(Y)))
  mu <- exp(eta)
  li <- fam$logpmf(Y, mu, theta)
  if (truncated) li <- li - log1p(-fam$p0(mu, theta))   # condition on Y >= 1
  li
}

## robust/clustered sandwich from a per-observation log-likelihood `ll_i(par)`.
##   type "analytic" -> inverse observed information (model-based)
##   type "robust"   -> HC sandwich  bread %*% meat %*% bread
##   type "cluster"  -> clustered meat (scores summed within cluster)
## Frequency weights w enter as duplicated rows: the information is that of
## sum(w * ll_i), the HC meat is sum(w_i s_i s_i'), and the cluster meat sums
## w_i s_i within clusters before the outer product.
## The derivatives are taken in the optimizer's scaled parameterization
## (par_s = par * scale, so a coefficient on a covariate in cents and one on a
## covariate in dollars are differenced at comparable steps) and mapped back,
## so the covariance does not depend on the covariates' units. The cluster
## meat carries the usual G/(G-1) small-sample factor (sandwich::vcovCL's
## default, and Stata's).
## `info`: the expected information in the scaled parameterization, which
## replaces the numerical Hessian for a family whose log-likelihood is kinked
## (.count_expected_info); its scores are then central differences at a small
## step, which a kink affects only for an observation within the step of it.
.count_vcov <- function(par, ll_i, type, cluster = NULL, nm, w = NULL, scale = NULL, info = NULL) {
  n <- length(ll_i(par)); w <- .ud_w1(w, n)
  k <- length(par); sc <- if (is.null(scale)) rep(1, k) else scale
  ps <- par * sc; to_par <- function(v) v / sc
  ## differenced at an absolute step of 1e-3 in the scaled parameterization
  ## (numDeriv's default is a step of a tenth of the parameter, which at large
  ## fitted means leaves the family's support and returns NA)
  ma <- list(eps = 1e-3, d = 0.1, zero.tol = 1e-8, r = 4, v = 2, show.details = FALSE)
  H <- if (!is.null(info)) -info
       else numDeriv::hessian(function(u) sum(w * ll_i(to_par(ps + u))), rep(0, k), method.args = ma)
  bread <- tryCatch(solve(-H), error = function(e) {
    ## the parameters the flattest direction loads on (a dispersion at its limit, say the negative binomial's
    ## log(theta) at the Poisson, or a dummy the data do not identify)
    ev <- eigen(-H, symmetric = TRUE); v <- ev$vectors[, which.min(abs(ev$values))]
    warning("the information is singular (flat in the direction of ", paste(nm[abs(v) > 0.3], collapse = ", "),
            "), so the covariance uses its generalized inverse; the standard errors of the parameters in that ",
            "direction are not meaningful.", call. = FALSE)
    MASS::ginv(-H)
  })
  V <- if (type == "analytic") {
    bread
  } else {
    S <- if (!is.null(info)) .count_cd_scores(ll_i, ps, to_par)
         else numDeriv::jacobian(function(u) ll_i(to_par(ps + u)), rep(0, k), method.args = ma)   # n x length(par) scores
    meat <- if (type == "cluster") {
      if (is.null(cluster)) stop("cluster-robust SE requires a 'cluster'.")
      G <- length(unique(cluster))
      (G / (G - 1)) * crossprod(as.matrix(rowsum(S * w, group = cluster)))   # sum weighted scores within cluster
    } else crossprod(S * sqrt(w))                            # HC (per-observation)
    bread %*% meat %*% bread
  }
  V <- V * tcrossprod(1 / sc)                                # back to the original parameterization
  dimnames(V) <- list(nm, nm)
  if (type != "analytic") attr(V, "model") <- structure(bread * tcrossprod(1 / sc), dimnames = list(nm, nm))
  V
}

## per-observation scores by central differences at the step `h` in the scaled parameterization (n x k)
.count_cd_scores <- function(ll_i, ps, to_par, h = 1e-6) {
  m <- length(ll_i(to_par(ps)))
  vapply(seq_along(ps), function(j) { e <- replace(numeric(length(ps)), j, h)
    (ll_i(to_par(ps + e)) - ll_i(to_par(ps - e))) / (2 * h) }, numeric(m))
}

## The expected information of a likelihood that is continuous but kinked in its parameters, in the scaled
## parameterization: sum_i w_i sum_y f_i(y) s_iy s_iy', the scores of every support point y of every observation
## differenced at a small step (a kink sits at one value of the parameters for a given y, so the step straddles
## it for almost no support point). A finite-difference Hessian of the summed log-likelihood instead straddles
## the kinks of the observations near them and measures only the local curvature of the pieces, which gives
## erratic and mostly too small standard errors. `llik(par, idx, y)` is the log-likelihood of observation idx[j]
## at the count y[j]; `ends` the last support point of each observation; `ymin` the first (1 when zero-truncated).
.count_expected_info <- function(par, llik, ends, ymin, w, scale) {
  n <- length(ends); w <- .ud_w1(w, n); ps <- par * scale; to_par <- function(v) v / scale
  len <- pmax(ends - ymin + 1L, 0L)
  idx <- rep(seq_len(n), len); yy <- unlist(lapply(seq_len(n), function(i) seq.int(ymin, length.out = len[i])))
  f <- exp(llik(par, idx, yy)); f[!is.finite(f)] <- 0
  S <- .count_cd_scores(function(p) llik(p, idx, yy), ps, to_par)
  S[!is.finite(S)] <- 0                                       # support points without mass
  crossprod(S * sqrt(f * w[idx]))
}
## the last support point of each fitted observation of a count family (its `kmax` rule, one per mean)
.count_ends <- function(fam, mu, theta) vapply(mu, function(m) as.integer(fam$kmax(m, theta)), integer(1))

## The distinct rows of a count fit: the rows with the same response, covariates, offset and (for a clustered
## covariance) cluster, as one row carrying the group's total weight. The log-likelihood, the information and the
## score sums over them are those of the full data up to the order of summation, so a design with repeated rows (an
## intercept-only fit, a few categorical covariates, tabulated data) is evaluated once per distinct row: the EWP's
## intercept-only fit to 5,414 rows of seven distinct counts takes about a second instead of minutes. NULL when every
## row is distinct. Rows are compared exactly, after sorting.
.count_collapse <- function(X, Y, off, w, cl = NULL) {
  n <- length(Y)
  if (n < 2L) return(NULL)
  cols <- c(list(as.numeric(Y), rep_len(as.numeric(off), n)), lapply(seq_len(ncol(X)), function(j) X[, j]),
            if (!is.null(cl)) list(as.character(cl)))
  o <- do.call(order, unname(cols))
  new <- c(TRUE, rep(FALSE, n - 1L))
  for (v in cols) { s <- v[o]; new[-1L] <- new[-1L] | (s[-1L] != s[-n]) }
  if (all(new)) return(NULL)
  id <- integer(n); id[o] <- cumsum(new)
  list(first = o[new], w = as.numeric(rowsum(.ud_w1(w, n), id, reorder = TRUE)))
}

## The saturated log-likelihood: the empirical distribution of the counts within each covariate pattern (the same
## covariates and offset, compared exactly), an upper bound of any count model's log-likelihood on these rows
.count_saturated <- function(X, Y, off, w) {
  n <- length(Y); w <- .ud_w1(w, n); off <- rep_len(as.numeric(off), n)
  grp <- function(cols) {                                    # group ids of the rows, by exact comparison after sorting
    o <- do.call(order, unname(cols)); new <- c(TRUE, rep(FALSE, n - 1L))
    for (v in cols) { s <- v[o]; new[-1L] <- new[-1L] | (s[-1L] != s[-n]) }
    id <- integer(n); id[o] <- cumsum(new); id
  }
  xc <- c(list(off), lapply(seq_len(ncol(X)), function(j) X[, j]))
  pat <- grp(xc); cell <- grp(c(xc, list(as.numeric(Y))))
  wc <- as.numeric(rowsum(w, cell, reorder = TRUE)); pc <- pat[match(seq_along(wc), cell)]
  wp <- as.numeric(rowsum(w, pat, reorder = TRUE))
  sum(wc[wc > 0] * log(wc[wc > 0] / wp[pc[wc > 0]]))
}

## core estimator shared by count_reg() and the intensity of hurdle/zi; optimized over the distinct rows
## (.count_collapse()), with the fitted means of every row.
.count_fit <- function(X, Y, fam, truncated, offset = 0, start = NULL, w = NULL) {
  p <- ncol(X); off <- rep_len(offset, nrow(X)); w <- .ud_w1(w, nrow(X))
  Xall <- X; off_all <- off
  cc <- .count_collapse(X, Y, off, w)
  if (!is.null(cc)) { X <- X[cc$first, , drop = FALSE]; Y <- Y[cc$first]; off <- off[cc$first]; w <- cc$w }
  s <- .ud_colscale(X, w); Xs <- sweep(X, 2, s, "/")           # optimize on unit-SD (weighted) columns
  sc <- c(s, rep(1, fam$nshape))
  if (!is.null(start)) start <- start * sc
  b0 <- tryCatch({ v <- stats::glm.fit(Xs, Y, weights = w, offset = off, family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  nll0 <- function(th) { v <- sum(w * .count_llik_i(th, Xs, Y, fam, truncated, off)); if (!is.finite(v)) 1e10 else -v }
  idx <- p + seq_len(fam$nshape)                               # the shape's coordinates (unscaled)
  nll <- function(th) nll0(.shape_clamp_search(th, fam, idx))   # a bounded shape held in the search range
  ## Shape families start at their Poisson value and at one underdispersed (or,
  ## for the negative binomial, near-Poisson) value, and keep the better
  ## optimum: on underdispersed data the negative-binomial likelihood is
  ## monotone toward its Poisson limit, and the COM-Poisson's BFGS can stall
  ## exactly at nu = 1 on real data whose optimum lies at nu > 1.
  starts <- if (!is.null(start)) list(start)
            else if (fam$nshape) lapply(fam$shape_starts, function(s) c(b0, fam$shape_link(s)))
            else list(b0)
  o <- NULL
  for (st in starts) {
    oi <- tryCatch(stats::optim(st, nll, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)),
                   error = function(e) NULL)
    if (!is.null(oi) && (is.null(o) || oi$value < o$value)) o <- oi
  }
  if (is.null(o)) stop("count model optimization failed from all starting values.")
  ## a Nelder-Mead polish from the BFGS solution, kept only if it improves; the fit is converged when either stage
  ## converged (a polish of a converged BFGS solution, or a polish that converged where BFGS ran out of iterations,
  ## as on the COM-Poisson's ridge to its supremum)
  if (length(o$par) > 1L) {
    o2 <- tryCatch(stats::optim(o$par, nll, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(o2) && o2$value < o$value) { o2$convergence <- min(o$convergence, o2$convergence); o <- o2 }
  }
  if (fam$nshape) {                                          # a shape that ran to the search's edge: the limit
    o$par <- .shape_clamp_search(o$par, fam, idx)
    o <- .shape_limit_candidate(o, nll0, fam, idx)
  }
  ## a fit at the saturated likelihood has nothing left to optimize, whatever the optimizers' stopping rules say (the
  ## COM-Poisson's ridge on a two-valued series approaches it as nu grows without bound, where neither BFGS nor the
  ## polish meets its relative tolerance)
  if (o$convergence != 0 && is.finite(o$value) && o$value < 1e9) {
    sat <- .count_saturated(X, Y, off, w)
    if (-o$value >= sat - 1e-6 * max(1, abs(sat))) o$convergence <- 0L
  }
  if (!is.finite(o$value) || o$value >= 1e9) {            # never left the infeasible plateau
    warning("the optimizer found no finite likelihood: a fitted mean lies outside [1e-12, 1e8] at every ",
            "point it visited. The fit did not converge; the coefficients are its last iterate.", call. = FALSE)
    o$value <- Inf; o$convergence <- 1L
  }
  par <- o$par / sc                                          # back to the original scale
  beta <- par[seq_len(p)]; names(beta) <- colnames(X)
  theta <- if (fam$nshape) unname(fam$shape_inv(par[p + seq_len(fam$nshape)])) else NA_real_
  list(par = par, beta = beta, theta = theta, loglik = -o$value, offset = off_all, scale = sc,
       mu = as.numeric(exp(off_all + Xall %*% beta)), converged = o$convergence == 0)
}

## resolve `cluster` to the estimation rows; NULL stays NULL
.count_cluster <- function(cluster, data, rows) {
  if (is.null(cluster)) return(NULL)
  cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
  if (is.null(cv)) stop("se = \"cluster\" needs 'cluster' (a column name in 'data' or a vector).")
  if (length(cv) != nrow(data)) stop("'cluster' must have one value per row of 'data'.")
  clid <- cv[rows]
  if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
  ng <- length(unique(clid))
  if (ng < 30L)
    warning(sprintf("only %d clusters; the cluster-robust SE carries only the G/(G - 1) small-sample correction and may be unreliable.", ng))
  clid
}

#' Matched count regressions: Poisson, negative binomial, COM-Poisson, generalized Poisson, gamma-count, double Poisson, hyper-Poisson, exponentially weighted Poisson
#'
#' Fits a count regression with a log link from any of the package's count
#' families, optionally zero-truncated and/or with unit fixed effects, returning
#' an object that [compare_models()], [score()], [dispersion_test()],
#' [dispersion_profile()], and the \pkg{broom} methods treat on the same footing
#' as a [cpb()] or [gec()] fit. The families share one interface, so a
#' gamma-count and a COM-Poisson can be compared to a CPB with identical
#' degrees-of-freedom, log-likelihood, proper-score, weight, and
#' robust-standard-error accounting rather than reconciled across packages.
#'
#' @section Families:
#' * `"poisson"`: the equidispersed baseline.
#' * `"negbin"`: negative binomial (overdispersion only); `$theta` is the size.
#' * `"compois"`: the classical Conway--Maxwell--Poisson in its rate
#'   parameterization, \eqn{\log\lambda = x'\beta}. The rate \eqn{\lambda} is
#'   not the mean (though `predict(type = "response")` and `fitted()` return the
#'   mean); `$theta` is the dispersion \eqn{\nu} (\eqn{\nu > 1} underdispersed,
#'   \eqn{\nu = 1} Poisson, \eqn{\nu < 1} overdispersed).
#' * `"mpcmp"`: the COM-Poisson in Huang's (2017) mean parameterization,
#'   \eqn{\log\mathrm{E}(Y) = x'\beta}, with the same dispersion \eqn{\nu}; the
#'   coefficients are effects on the log mean and rate ratios are exact.
#' * `"genpois"`: the Consul--Jain generalized Poisson with constant dispersion
#'   \eqn{\lambda \in (-1, 1)} and \eqn{\mu = \exp(x'\beta)}, the mean for
#'   \eqn{\lambda \ge 0}, with \eqn{\mathrm{Var}/\mathrm{Mean} = 1/(1-\lambda)^2};
#'   \eqn{\lambda < 0} is underdispersion, on a finite support that the pmf is
#'   renormalized over, and the mean then departs from \eqn{\mu} where that
#'   support is short (small \eqn{\mu} at strong underdispersion: about
#'   \eqn{1.19\mu} at \eqn{\lambda = -0.6}, \eqn{\mu = 0.5}); the exact mean is
#'   reported (see [genpois-distribution]).
#' * `"gammacount"`: Winkelmann's (1995) gamma-count renewal-process model,
#'   underdispersed when the waiting times between events are more regular than
#'   exponential (\eqn{\alpha > 1}); \eqn{\exp(x'\beta)} is the long-run event
#'   rate and the exact mean is reported (see [gammacount-distribution]).
#' * `"doublepois"`: Efron's (1986) double Poisson with the exact normalizing
#'   constant, \eqn{\theta > 1} underdispersed (see [doublepois-distribution]).
#' * `"hyperpois"`: the hyper-Poisson of Bardwell and Crow (1964) in the mean
#'   parameterization of \enc{Sáez-Castillo}{Saez-Castillo} and
#'   \enc{Conde-Sánchez}{Conde-Sanchez} (2013), mean \eqn{\mu = \exp(x'\beta)},
#'   \eqn{\gamma < 1} underdispersed, with one \eqn{\gamma} for all
#'   observations: their constant-dispersion model hP\eqn{(\mu_i, \gamma)}. Their
#'   model with covariates in the dispersion, hP\eqn{(\mu_i, \gamma_i)}, is fit by
#'   `DGLMExtPois::glm.hP()`. Its underdispersion is bounded: as
#'   \eqn{\gamma \to 0} the variance falls only to \eqn{\mu - 1} for
#'   \eqn{\mu \ge 1} and \eqn{\mu(1 - \mu)} below, so Var/Mean is at least
#'   0.5 at a mean of 2 and 0.9 at 10 (see [hyperpois-distribution]).
#' * `"ewp3"` (alias `"ewp"`): the three-parameter exponentially weighted
#'   Poisson of Ridout and Besbeas (2004), a Poisson with centrality
#'   \eqn{\lambda = \exp(x'\beta)} reweighted by a two-sided exponential in the
#'   distance from it; `$theta` holds the two dispersion parameters
#'   \eqn{(\beta_1, \beta_2)}, one per tail, positive values underdispersed
#'   (see [ewp-distribution]). `"ewp2"` is its two-parameter form,
#'   \eqn{\beta_1 = \beta_2 = \beta}. The coefficients are effects on the log
#'   centrality, not the log mean; the exact mean is reported.
#'
#' The Poisson is the equidispersed member of every family, so
#' [dispersion_test()] tests each family's dispersion parameter against it, and
#' [dispersion_profile()] compares the families' implied variance-to-mean curves
#' against the data.
#'
#' Runtime: the Poisson, negative binomial, gamma-count, and double Poisson and
#' generalized Poisson families fit in under a second at a few hundred rows; the
#' COM-Poisson families, the hyper-Poisson and the EWP evaluate a normalizing
#' sum per observation (and the mean parameterizations solve a root per
#' observation), so they take seconds at a few hundred rows and minutes at
#' several thousand. The EWP's sum lengthens with a heavy upper tail (a
#' negative upper weight), and its standard errors sum over each observation's
#' support. Rows that repeat (the same count, covariates and offset) are
#' evaluated once with their total weight, so an intercept-only or
#' few-category design costs as many evaluations as it has distinct rows.
#' `se = "robust"` and `"cluster"` add numerical scores at the same cost per
#' parameter.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param family One of `"poisson"`, `"negbin"`, `"compois"`, `"mpcmp"`,
#'   `"genpois"`, `"gammacount"`, `"doublepois"`, `"hyperpois"`, `"ewp3"`,
#'   `"ewp2"` (see Families).
#' @param truncated Logical; if `TRUE`, fit the zero-truncated form (requires all
#'   `Y >= 1`), the count analogue of the CPB's zero-truncated default.
#' @param fe Optional column name(s) for fixed effects, entered as factor dummies
#'   (so the degrees of freedom count each absorbed intercept, matching
#'   [cpb_fe()]). Pass a vector of two columns for two-way (e.g. unit and time)
#'   fixed effects.
#' @param offset Optional offset entered on the linear-predictor (log) scale --
#'   the log-mean for the mean-parameterized families, the log-rate
#'   (\eqn{\log\lambda}) for `"compois"`. A numeric vector or the name of a
#'   column in `data`, e.g. \eqn{\log(\text{exposure})} so the model becomes a
#'   rate model.
#' @param weights Optional frequency weights: a numeric vector or the name of a
#'   column in `data`. A weight of \eqn{w_i} is equivalent to \eqn{w_i} copies of
#'   row \eqn{i} in the log-likelihood, the information, and the score products,
#'   and `nobs()` returns the weight total. Survey (probability) weights call
#'   for `se = "robust"` or `se = "cluster"`.
#' @param se Standard errors: `"analytic"` (inverse information, default),
#'   `"robust"` (heteroskedasticity-consistent sandwich), `"cluster"`
#'   (cluster-robust; needs `cluster`), or `"none"`. A robust or clustered
#'   covariance protects inference about the parameters of the best-fitting
#'   member of the family; the quantities that rest on the whole distribution
#'   (predicted probabilities, the mean of a truncated or two-part model) are
#'   still those of the fitted family. King and Roberts (2015) read a large gap
#'   between the robust and the model-based standard errors as evidence of a
#'   misspecification to fix; [summary()] reports their ratio. For the Poisson
#'   the coefficients are consistent for the conditional mean whenever
#'   \eqn{E(Y \mid x) = \exp(x'\beta)}, whatever the distribution
#'   (Gourieroux, Monfort and Trognon 1984; Wooldridge 1999), and the robust or
#'   clustered errors are then valid for them: the gap shows that the variance
#'   is not the Poisson's, not that the mean is wrong. That makes the robust
#'   Poisson the benchmark for mean effects, against which the families' own
#'   mean effects ([ame()]) can be compared; a family that fits the
#'   distribution better estimates them consistently only if it is the true
#'   model.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector aligned to its rows) for `se = "cluster"`; supplying it selects
#'   `se = "cluster"` unless `se` is given explicitly.
#' @return An object of class `"count_reg"`. The component `$theta` holds the
#'   dispersion parameter on its natural scale (the negative-binomial size, the
#'   COM-Poisson \eqn{\nu}, the generalized-Poisson \eqn{\lambda}, the
#'   gamma-count \eqn{\alpha}, the double-Poisson \eqn{\theta}, the
#'   hyper-Poisson \eqn{\gamma}, the EWP's \eqn{(\beta_1, \beta_2)} or
#'   \eqn{\beta}; `NA` for the Poisson); `$fitted.values` and `$residuals` are on the mean scale (the
#'   conditional mean E(Y | Y >= 1) for a zero-truncated fit), as `fitted()` and
#'   `residuals()` return them, and `$mu` is the family's natural parameter
#'   `exp(offset + x'b)`.
#' @references King, G. and Roberts, M. E. (2015). How robust standard errors
#'   expose methodological problems they do not fix, and what to do about it.
#'   \emph{Political Analysis}, 23(2), 159-179.
#'   Huang, A. (2017). Mean-parametrized Conway--Maxwell--Poisson
#'   regression models for dispersed counts. \emph{Statistical Modelling},
#'   17(6), 359-380. Winkelmann, R. (1995). Duration dependence and dispersion
#'   in count-data models. \emph{Journal of Business & Economic Statistics},
#'   13(4), 467-474. Efron, B. (1986). Double exponential families and their
#'   use in generalized linear regression. \emph{Journal of the American
#'   Statistical Association}, 81(395), 709-721. Consul, P. C. and Famoye, F.
#'   (1992). Generalized Poisson regression model. \emph{Communications in
#'   Statistics -- Theory and Methods}, 21(1), 89-109. Shmueli, G., Minka, T.
#'   P., Kadane, J. B., Borle, S. and Boatwright, P. (2005). A useful
#'   distribution for fitting discrete data: revival of the Conway--Maxwell--Poisson
#'   distribution. \emph{Journal of the Royal Statistical Society C}, 54(1),
#'   127-142. Sellers, K. F. and Shmueli, G. (2010). A flexible regression model
#'   for count data. \emph{Annals of Applied Statistics}, 4(2), 943-961.
#'   Bardwell, G. E. and Crow, E. L. (1964). A two-parameter family of
#'   hyper-Poisson distributions. \emph{Journal of the American Statistical
#'   Association}, 59(305), 133-141. \enc{Sáez-Castillo}{Saez-Castillo}, A. J.
#'   and \enc{Conde-Sánchez}{Conde-Sanchez}, A. (2013). A hyper-Poisson
#'   regression model for overdispersed and underdispersed count data.
#'   \emph{Computational Statistics & Data Analysis}, 61, 148-157. Ridout, M.
#'   S. and Besbeas, P. (2004). An empirical model for underdispersed count
#'   data. \emph{Statistical Modelling}, 4(1), 77-89. Gourieroux, C., Monfort,
#'   A. and Trognon, A. (1984). Pseudo maximum likelihood methods: applications
#'   to Poisson models. \emph{Econometrica}, 52(3), 701-720. Wooldridge, J. M.
#'   (1999). Distribution-free estimation of some nonlinear panel data models.
#'   \emph{Journal of Econometrics}, 90(1), 77-97.
#' @seealso [cpb()], [gec()], [compare_models()], [dispersion_test()],
#'   [dispersion_profile()], [hurdle_count()], [zi_count()]
#' @examples
#' set.seed(1); n <- 400; x <- rnorm(n)
#' y <- rgammacount(n, exp(1 + 0.5 * x), alpha = 2)     # regular event timing
#' d <- data.frame(y = y, x = x)
#' m <- count_reg(y ~ x, data = d, family = "gammacount")
#' m
#' compare_models(poisson = count_reg(y ~ x, d, family = "poisson"), gammacount = m)
#' @export
count_reg <- function(formula, data, family = c("poisson", "negbin", "compois", "mpcmp", "genpois",
                                                "gammacount", "doublepois", "hyperpois", "ewp3", "ewp2"),
                      truncated = FALSE, fe = NULL, offset = NULL, weights = NULL,
                      se = c("analytic", "robust", "cluster", "none"), cluster = NULL) {
  .ud_no_formula_offset(formula)
  data <- .ud_drop_na_fe(data, fe)
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  se_missing <- missing(se)
  fam <- .count_fam(family); family <- fam$tag; se <- match.arg(se)
  ## uniform cluster semantics across the family: supplying `cluster` implies
  ## clustered inference unless the user explicitly chose another se type
  if (!is.null(cluster) && se_missing) se <- "cluster"
  cl <- match.call()
  if (!is.null(fe)) {
    miss <- setdiff(fe, names(data))
    if (length(miss)) stop("'fe' names column(s) not in 'data': ", paste(miss, collapse = ", "), ".")
    if (is.character(cluster) && any(cluster %in% fe))
      warning("'cluster' is the fixed-effects variable: each unit dummy is identified within a single cluster, ",
              "so the cluster-robust covariance is degenerate for the dummies (the slopes' cluster-robust ",
              "standard errors remain valid).", call. = FALSE)
    formula <- .ud_add_fe_terms(formula, fe, data)               # fe may name >1 column (two-way FE)
  }
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit, drop.unused.levels = TRUE)
  Y <- stats::model.response(mf); X <- stats::model.matrix(formula, mf)
  n <- length(Y); p <- ncol(X); rows <- .ud_kept_rows(mf, data)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  .ud_warn_all_zero(Y)
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")
  .ud_rank_check(X)
  off <- 0                                                   # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    if (is.null(ov)) stop("'offset' must name a column of 'data' or be a numeric vector.")
    if (length(ov) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    off <- as.numeric(ov[rows])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  w <- .ud_weights(weights, data, rows)

  if (!is.null(cluster) && se != "cluster")
    warning("'cluster' is ignored unless se = \"cluster\".")
  clid <- if (se == "cluster") .count_cluster(cluster, data, rows) else NULL

  ft <- .count_fit(X, Y, fam, truncated, offset = off, w = w)
  df <- p + fam$nshape
  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL; vcv_full <- NULL; se_model <- NULL
  if (se != "none") {
    nm <- c(colnames(X), if (fam$nshape) fam$shape_name)
    ## over the distinct rows (within clusters for a clustered covariance; see .count_collapse())
    offn <- rep_len(off, n); cc <- .count_collapse(X, Y, offn, w, clid)
    r <- if (is.null(cc)) seq_len(n) else cc$first
    Xv <- X[r, , drop = FALSE]; Yv <- Y[r]; offv <- offn[r]; clv <- clid[r]
    wv <- if (is.null(cc)) w else cc$w
    info <- if (isTRUE(fam$kinked))
      .count_expected_info(ft$par, function(p, idx, y) .count_llik_i(p, Xv[idx, , drop = FALSE], y, fam, truncated, offv[idx]),
                           .count_ends(fam, ft$mu[r], ft$theta), if (truncated) 1L else 0L, wv, ft$scale)
    V <- tryCatch(.count_vcov(ft$par, function(th) .count_llik_i(th, Xv, Yv, fam, truncated, offv),
                              se, clv, nm, w = wv, scale = ft$scale, info = info), error = function(e) .count_vcov_failed(e))
    if (!is.null(V)) {
      vm <- attr(V, "model"); attr(V, "model") <- NULL
      vcv_full <- V
      vcv <- V[seq_len(p), seq_len(p), drop = FALSE]
      se.beta <- .count_se(vcv)
      if (!is.null(vm)) se_model <- sqrt(pmax(diag(vm), 0))[seq_len(p)]
    }
  }
  ## $mu is the family's natural parameter (the rate of the recursion); the
  ## stored fitted values and residuals are on the mean scale, as fitted() returns
  fv <- fam$meanfun(ft$mu, ft$theta)
  if (truncated) fv <- fv / pmax(1 - fam$p0(ft$mu, ft$theta), 1e-12)     # E(Y | Y >= 1)
  structure(list(
    coefficients = ft$beta, theta = ft$theta, family = family, truncated = truncated,
    se.beta = se.beta, se_model = se_model, vcov = vcv, vcov_full = vcv_full, loglik = ft$loglik, mu = ft$mu,
    fitted.values = fv,
    linear.predictors = as.numeric(rep_len(off, n) + X %*% ft$beta), residuals = Y - fv,
    offset = off, offset_name = .ud_offset_name(offset), weights = w, n = n, nobs_weighted = if (is.null(w)) n else sum(w),
    df = df, p = p, fe = fe, se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = ft$converged, formula = formula, terms = attr(mf, "terms"),
    levels = .cpb_xlevels(mf), contrasts = attr(X, "contrasts"), X = X, Y = Y, n_dropped = nrow(data) - n, call = cl),
    class = "count_reg")
}

## a covariance the numerical derivatives could not form: the fit says so and carries no standard errors
.count_vcov_failed <- function(e) {
  warning("the covariance could not be computed (", conditionMessage(e), "); the standard errors are NA.", call. = FALSE)
  NULL
}

## standard errors from a numerical covariance: a non-positive diagonal means
## the Hessian is not positive definite (the fit may sit on a feasibility
## boundary, with counts at their ceiling), and such an SE is NA, not 0
.count_se <- function(vcv) {
  d <- diag(vcv); bad <- !is.finite(d) | d <= 0
  if (any(bad)) {
    warning("the information is not positive definite for ", sum(bad), " coefficient(s) (",
            paste(names(d)[bad][seq_len(min(5L, sum(bad)))], collapse = ", "), if (sum(bad) > 5L) ", ..." else "",
            "): their standard errors are NA. This happens where the likelihood is flat in a parameter: a ",
            "fixed-effects unit with no positive count (its dummy at minus infinity), or a fit on a boundary of the ",
            "family's support or parameter space. The other standard errors are unaffected.", call. = FALSE)
    d[bad] <- NA_real_
  }
  sqrt(d)
}

## ---------------------------------------------------------------------------
## Hurdle count regression (participation logit + zero-truncated intensity).
## The likelihood factorizes, so the two margins are fit separately and each
## carries its own standard errors, exactly as hurdle_cpb().
## ---------------------------------------------------------------------------

#' Hurdle count regression for the matched count families
#'
#' Fits a participation model joined to a zero-truncated intensity from any of
#' the package's count families (see [count_reg()], Families), the count
#' analogue of [hurdle_cpb()]. Returned as a `"hurdle_count"` object that
#' [compare_models()] and [score()] accept, so a hurdle-CPB and a hurdle
#' gamma-count can be compared on one footing.
#'
#' @param formula Intensity formula (`y ~ x`).
#' @param data A data frame.
#' @param family The intensity family; see [count_reg()].
#' @param participation Optional one-sided formula for the participation model;
#'   defaults to the intensity right-hand side.
#' @param fe,part_fe Optional fixed-effects column names for the intensity and
#'   participation models (entered as factor dummies). The participation
#'   slopes carry the dummy logit's incidental-parameters bias in short panels
#'   (see [hurdle_cpb()]). The participation dummy of a unit that never
#'   participates is not counted in the degrees of freedom (see
#'   [compare_models()], Degrees of freedom). In new data such a unit, and one
#'   whose every count is positive, get the participation probability at their
#'   limits (zero and one; with one `part_fe` column), and a unit either
#'   equation has not seen (the intensity never sees a unit with no positive
#'   count) enters that equation at its reference unit, as in [hurdle_cpb()];
#'   the participation equation is checked for separation (`part_separation`).
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`. Positive participation coefficients raise the probability of
#'   a *positive count* (participation).
#' @param offset Optional offset for the intensity, on the linear-predictor (log)
#'   scale: a numeric vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name),
#'   applied to both margins; see [count_reg()].
#' @param se,cluster Standard-error type and optional cluster; see
#'   [count_reg()]. They apply to both margins: the participation equation's
#'   standard errors are model-based, robust or cluster-robust with the
#'   intensity's.
#' @return An object of class `"hurdle_count"`.
#' @examples
#' set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(0.4 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
#' hurdle_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
#'              participation = ~ z, se = "none")
#' @seealso [hurdle_cpb()], [zi_count()], [compare_models()]
#' @export
hurdle_count <- function(formula, data, family = c("poisson", "negbin", "compois", "mpcmp", "genpois",
                                                   "gammacount", "doublepois", "hyperpois", "ewp3", "ewp2"),
                         participation = NULL, fe = NULL, part_fe = NULL,
                         link = c("logit", "probit", "cloglog"), offset = NULL, weights = NULL,
                         se = c("analytic", "robust", "cluster", "none"), cluster = NULL) {
  .ud_no_formula_offset(formula, participation)
  se_missing <- missing(se)
  fam <- .count_fam(family); family <- fam$tag; se <- match.arg(se); link <- match.arg(link)
  if (!is.null(cluster) && se_missing) se <- "cluster"       # uniform cluster semantics
  ## reduce to complete cases on all model variables so the two margins stay aligned
  mv <- unique(c(all.vars(formula), all.vars(if (is.null(participation)) formula[-2L] else participation),
                 fe, part_fe))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  if (!is.null(offset) && !is.character(offset)) offset <- offset[keep]
  w_all <- .ud_weights(weights, data, seq_len(nrow(data)))
  if (!is.null(w_all)) w_all <- w_all[keep]
  cl_vals <- if (se == "cluster") .ud_cluster_values(cluster, data, keep) else NULL   # checked before the rows are reduced
  data <- data[keep, , drop = FALSE]
  y <- stats::model.response(stats::model.frame(formula, data))
  if (!is.numeric(y)) stop("Response must be a numeric count; got ", class(y)[1L], ".")
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  d <- as.integer(y > 0)
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d; pdata[[".w"]] <- if (is.null(w_all)) rep(1, nrow(data)) else w_all
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, weights = .w,
                     family = stats::binomial(link = link))
  pf <- .part_fields(pfit, part_fe, data, switch(se, cluster = "cluster", robust = "robust", "model"), cl_vals,
                     "hurdle_count")
  pos <- data[y > 0, , drop = FALSE]
  int_offset <- if (is.null(offset)) NULL                       # offset applies to the intensity
                else if (is.character(offset) && length(offset) == 1L) offset
                else offset[y > 0]
  ifit <- count_reg(formula, data = pos, family = family, truncated = TRUE,
                    fe = fe, offset = int_offset, weights = if (is.null(w_all)) NULL else w_all[y > 0],
                    se = se, cluster = if (is.null(cl_vals) || (is.character(cluster) && length(cluster) == 1L)) cluster
                                       else cl_vals[y > 0])
  ## per-observation intensity natural parameter for all units (missing FE levels
  ## -> reference), offset included
  beta <- ifit$coefficients
  mm <- stats::model.matrix(stats::delete.response(stats::terms(ifit$formula)), data)
  common <- intersect(names(beta), colnames(mm))
  off_all <- if (is.null(offset)) 0 else if (is.character(offset)) as.numeric(data[[offset]]) else as.numeric(offset)
  lambda_full <- exp(off_all + as.numeric(mm[, common, drop = FALSE] %*% beta[common]))
  llh <- as.numeric(stats::logLik(pfit)) + ifit$loglik
  zl <- .ud_fe_zero_levels(y, data, part_fe)   # never participate: not counted (.ud_zero_levels)
  ol <- .ud_fe_one_levels(y, data, part_fe)    # always participate: predicted at the limit
  dfh <- length(stats::coef(pfit)) - length(zl) + ifit$df
  structure(c(list(participation = pfit, intensity = ifit,
                 converged = isTRUE(pfit$converged) && isTRUE(ifit$converged), y = as.integer(y),
                 p_full = as.numeric(stats::fitted(pfit)), lambda_full = lambda_full,
                 theta = ifit$theta, family = family, loglik = llh, df = dfh, part_zero_units = zl,
                 part_one_units = ol, link = link,
                 int_fe_ref = mean(off_all),
                 int_beta = beta,
                 int_xref = .fe_dummy_ref(colMeans(mm[, common, drop = FALSE]), mm[, common, drop = FALSE], y, fe,
                                          .ident_trunc),
                 part_beta = stats::coef(pfit),
                 part_xref = .fe_dummy_ref(colMeans(stats::model.matrix(pfit)), stats::model.matrix(pfit), d, part_fe,
                                           .ident_binary),
                 n = length(y), nobs_weighted = if (is.null(w_all)) length(y) else sum(w_all),
                 weights = w_all, n_participate = sum(d), fe = fe, part_fe = part_fe, se.type = se,
                 formula = formula, part_formula = participation, call = match.call()), pf),
            class = "hurdle_count")
}

## ---------------------------------------------------------------------------
## Zero-inflated count regression (structural-zero mixture). Unlike the hurdle
## this does not factorize, so it is a single joint maximization.
## ---------------------------------------------------------------------------

## per-observation log-likelihood in (beta_count, gamma_zero, shape). The
## structural-zero probability uses an arbitrary link (logit/probit/cloglog).
.zi_llik_i <- function(par, Xc, Zz, Y, fam, linkinv = stats::plogis, offset = 0) {
  pc <- ncol(Xc); pz <- ncol(Zz)
  eta <- offset + as.numeric(Xc %*% par[seq_len(pc)])
  theta <- if (fam$nshape) fam$shape_inv(par[pc + pz + seq_len(fam$nshape)]) else NA_real_
  if (.count_excursion(eta, theta, fam)) return(rep(-Inf, length(Y)))
  mu <- exp(eta)
  pistar <- linkinv(as.numeric(Zz %*% par[pc + seq_len(pz)]))
  f0 <- fam$p0(mu, theta)
  ifelse(Y == 0, log(pistar + (1 - pistar) * f0), log1p(-pistar) + fam$logpmf(Y, mu, theta))
}

#' Zero-inflated count regression for the matched count families
#'
#' Fits a structural-zero mixture whose count component comes from any of the
#' package's count families (see [count_reg()], Families), the count analogue
#' of [zi_cpb()]. The count component uses a log link (on the mean for the
#' mean-parameterized families, on the rate \eqn{\lambda} for `"compois"`) and
#' the inflation probability a link set by `link`. Returned as a `"zi_count"`
#' object that [compare_models()] and [score()] accept.
#'
#' @param formula Count formula (`y ~ x`).
#' @param data A data frame.
#' @param family The count-component family; see [count_reg()].
#' @param zero Optional one-sided formula for the inflation (structural-zero)
#'   model; defaults to the count right-hand side.
#' @param fe Optional fixed-effects column for the count equation (factor dummies).
#' @param zero_fe Optional fixed-effects column for the inflation equation
#'   (factor dummies); zero-equation fixed effects are opt-in. In either
#'   equation the dummy of a unit with no positive count is not counted in the
#'   degrees of freedom (see [compare_models()], Degrees of freedom). In new
#'   data (with one column each) such a unit gets the structural-zero
#'   probability one through `zero_fe` and a count mean of zero through `fe`,
#'   and a unit whose every count is positive the structural-zero probability
#'   zero.
#' @param link Link for the inflation probability: `"logit"` (default),
#'   `"probit"`, or `"cloglog"`. Positive inflation coefficients raise the
#'   probability of a *structural zero*.
#' @param offset Optional offset for the count component, on the linear-predictor
#'   (log) scale: a numeric vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [count_reg()].
#' @param se,cluster Standard-error type and optional cluster; see [count_reg()].
#' @section Separation:
#' The fit checks its zero-inflation equation for separation as [zi_cpb()] does
#' (see its Separation section), under its own `link`. It warns when the data put
#' no useful bound on a covariate coefficient, and stores the check in the
#' component `separation`. Such a coefficient's standard error is not meaningful
#' either.
#' @return An object of class `"zi_count"`; `$fitted.values` is the marginal
#'   mean `(1 - pi) E(Y | count component)` that `fitted()` returns, and `$mu`
#'   the count component's natural parameter.
#' @examples
#' set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
#' zi_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
#'          zero = ~ z, se = "none")
#' @seealso [zi_cpb()], [hurdle_count()], [compare_models()]
#' @export
zi_count <- function(formula, data, family = c("poisson", "negbin", "compois", "mpcmp", "genpois",
                                               "gammacount", "doublepois", "hyperpois", "ewp3", "ewp2"),
                     zero = NULL, fe = NULL, zero_fe = NULL, link = c("logit", "probit", "cloglog"),
                     offset = NULL, weights = NULL, se = c("analytic", "robust", "cluster", "none"),
                     cluster = NULL) {
  .ud_no_formula_offset(formula, zero)
  se_missing <- missing(se)
  fam <- .count_fam(family); family <- fam$tag; se <- match.arg(se); link <- match.arg(link)
  if (!is.null(cluster) && se_missing) se <- "cluster"   # uniform cluster semantics
  linkinv <- stats::make.link(link)$linkinv
  ## a '.' expanded over the covariates: the columns other than the unit identifiers of either equation
  fexp <- .ud_add_fe_terms(formula, NULL, data[, setdiff(names(data), c(fe, zero_fe)), drop = FALSE])
  cform <- if (is.null(fe)) fexp else .ud_add_fe_terms(fexp, fe, data)
  zrhs <- if (is.null(zero)) fexp[-2L] else zero                 # the covariates, not the unit dummies
  if (!is.null(zero_fe))                                        # zero-equation FE is opt-in
    zrhs <- stats::reformulate(c(labels(stats::terms(zrhs, data = data)), paste0("factor(", zero_fe, ")")),
                               env = environment(zrhs))
  ## one model frame for both equations, with the count formula's own response expression
  allform <- stats::reformulate(unique(c(labels(stats::terms(cform)), labels(stats::terms(zrhs, data = data)))),
                                response = cform[[2L]], env = environment(cform))
  mf <- stats::model.frame(allform, data, na.action = stats::na.omit, drop.unused.levels = TRUE)
  Y <- stats::model.response(mf); Xc <- stats::model.matrix(cform, mf); Zz <- stats::model.matrix(zrhs, mf)
  n <- length(Y); pc <- ncol(Xc); pz <- ncol(Zz); rows <- .ud_kept_rows(mf, data)
  if (!is.numeric(Y)) stop("Response must be a numeric count; got ", class(Y)[1L], ".")
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  .ud_warn_all_zero(Y)
  .ud_rank_check(Xc, "count design"); .ud_rank_check(Zz, "inflation design")
  if (!is.null(cluster) && se != "cluster")
    warning("'cluster' is ignored unless se = \"cluster\".")
  clid <- if (se == "cluster") .count_cluster(cluster, data, rows) else NULL
  off <- 0                                                     # offset enters the count component
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    if (is.null(ov)) stop("'offset' must name a column of 'data' or be a numeric vector.")
    if (length(ov) != nrow(data)) stop("'offset' must have one value per row of 'data'.")
    off <- as.numeric(ov[rows])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  w <- .ud_weights(weights, data, rows); w1 <- .ud_w1(w, n)
  b0 <- tryCatch({ v <- stats::glm.fit(Xc, Y, weights = w1, offset = rep_len(off, n), family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, pc))
  g0 <- tryCatch({ v <- stats::glm.fit(Zz, as.integer(Y == 0), weights = w1, family = stats::binomial(link = link))$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, pz))
  ## pscl-style seeding: also size the count component to the POSITIVE counts via
  ## a zero-truncated fit, so heavy inflation cannot drag the starting rate toward
  ## a poor basin; keep the all-data Poisson start as a fallback and take the
  ## better optimum.
  st_tr <- tryCatch({
    pos <- Y > 0
    ft  <- .count_fit(Xc[pos, , drop = FALSE], Y[pos], fam, truncated = TRUE,
                      offset = rep_len(off, n)[pos], w = w1[pos])
    c(ft$beta, g0, if (fam$nshape) fam$shape_link(ft$theta))
  }, error = function(e) NULL)
  nll0 <- function(th) { v <- sum(w1 * .zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off)); if (!is.finite(v)) 1e10 else -v }
  idx <- pc + pz + seq_len(fam$nshape)
  nll <- function(th) nll0(.shape_clamp_search(th, fam, idx))
  ## the optimizer works on the covariates of both equations scaled to unit
  ## (weighted) standard deviation, and each start is polished by chained
  ## Nelder-Mead restarts; the coefficients are mapped back at the end
  sc <- c(.ud_colscale(Xc, w1), .ud_colscale(Zz, w1), rep(1, fam$nshape))
  nll_s <- function(v) nll(v / sc)
  ## starts that span the structural share of zeros, where the mixture
  ## likelihood's optima differ (see zi_cpb()): beside the two above, the
  ## zero-inflated Poisson's coefficients (logit link) and a low-inflation start,
  ## the all-data Poisson with a constant inflation probability of 0.02
  sh0 <- if (fam$nshape) fam$shape_link(fam$shape_starts[[1L]])
  st_zip <- if (link != "logit") NULL else tryCatch({
    zs <- .zi_zip_seeds(Xc, Zz, as.integer(Y), rep_len(off, n), w1, support = "ceiling")
    if (is.null(zs)) NULL else c(zs$b, zs$g, sh0) }, error = function(e) NULL)
  icz <- match("(Intercept)", colnames(Zz))
  st_low <- if (is.na(icz)) NULL else c(b0, replace(rep(0, pz), icz, stats::make.link(link)$linkfun(0.02)), sh0)
  ## the best optimum over a list of starts (original scale), mapped back to the original scale
  run <- function(starts) {
    o <- NULL
    for (st in starts) {
      oi <- tryCatch(stats::optim(st * sc, nll_s, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10)),
                     error = function(e) NULL)
      if (!is.null(oi)) oi <- .ud_nm_polish(oi, nll_s, maxit = 5000L, reltol = 1e-10)
      if (!is.null(oi) && (is.null(o) || oi$value < o$value)) o <- oi
    }
    if (is.null(o)) return(NULL)
    if (fam$nshape) {                                        # the shape's coordinates are unscaled (sc = 1)
      o$par <- .shape_clamp_search(o$par, fam, idx)
      o <- .shape_limit_candidate(o, function(v) nll0(v / sc), fam, idx)
    }
    o$par <- o$par / sc
    o
  }
  o <- run(Filter(Negate(is.null), list(st_tr, c(b0, g0, sh0), st_zip, st_low)))
  if (is.null(o)) stop("zi_count: optimization failed from all starting values.")
  ## the separation check, and a restart from any better point it finds (see ?zi_cpb, Separation)
  chk <- .zi_restart(o,
    check = function(o) {
      th <- if (fam$nshape) unname(fam$shape_inv(o$par[pc + pz + seq_len(fam$nshape)])) else NA_real_
      .zi_separation(Zz, o$par[pc + seq_len(pz)], Y, fam$p0(as.numeric(exp(rep_len(off, n) + Xc %*% o$par[seq_len(pc)])), th),
                     w, .ud_fe_group(data[rows, , drop = FALSE], zero_fe), .zi_fe_cols(Zz, zero_fe), link = link)
    },
    refit = function(o, g) run(list(c(o$par[seq_len(pc)], g, o$par[pc + pz + seq_len(fam$nshape)]))),
    loglik = function(o) -o$value)
  o <- chk$fit; sep <- chk$separation
  .zi_warn_separation(sep, "zi_count")
  beta <- o$par[seq_len(pc)]; names(beta) <- colnames(Xc)
  gamma <- o$par[pc + seq_len(pz)]; names(gamma) <- colnames(Zz)
  theta <- if (fam$nshape) unname(fam$shape_inv(o$par[pc + pz + seq_len(fam$nshape)])) else NA_real_
  ## the count effects and zero-equation dummies of units with no positive count are at their limits; not counted
  ## (.ud_zero_levels)
  zl_fe <- unlist(lapply(fe, function(f) .ud_zero_levels(Y, data[[f]][rows])))
  zl_zero <- .ud_fe_zero_levels(Y, data[rows, , drop = FALSE], zero_fe)
  ol_zero <- .ud_fe_one_levels(Y, data[rows, , drop = FALSE], zero_fe)   # the dummy at its other limit, in new data
  df <- pc + pz + fam$nshape - length(zl_fe) - length(zl_zero)
  se.beta <- setNames(rep(NA_real_, pc), colnames(Xc)); se.zero <- setNames(rep(NA_real_, pz), colnames(Zz))
  vcv <- NULL; vcv_full <- NULL
  if (se != "none") {
    nm <- c(paste0("count:", colnames(Xc)), paste0("zero:", colnames(Zz)), if (fam$nshape) fam$shape_name)
    info <- if (isTRUE(fam$kinked)) {
      offn <- rep_len(off, n)
      .count_expected_info(o$par, function(p, idx, y) .zi_llik_i(p, Xc[idx, , drop = FALSE], Zz[idx, , drop = FALSE], y, fam,
                                                                 linkinv, offn[idx]),
                           .count_ends(fam, as.numeric(exp(offn + Xc %*% beta)), theta), 0L, w, sc)
    }
    V <- tryCatch(.count_vcov(o$par, function(th) .zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off), se, clid, nm, w = w,
                              scale = sc, info = info),
                  error = function(e) .count_vcov_failed(e))
    if (!is.null(V)) {
      vcv_full <- V
      vcv <- V[seq_len(pc), seq_len(pc), drop = FALSE]; se.beta <- .count_se(vcv)
      se.zero <- .count_se(V[pc + seq_len(pz), pc + seq_len(pz), drop = FALSE])
    }
  }
  structure(list(coefficients = beta, zero.coefficients = gamma, zero_coef = gamma, theta = theta, family = family,
                 se.beta = se.beta, se.zero = se.zero, vcov = vcv, vcov_full = vcv_full,
                 loglik = -o$value, df = df, fe_zero_units = zl_fe, zero_fe_zero_units = zl_zero,
                 zero_fe_one_units = ol_zero, n = n,
                 nobs_weighted = if (is.null(w)) n else sum(w),
                 weights = w, link = link,
                 pi_full = as.numeric(linkinv(Zz %*% gamma)), separation = sep,
                 lambda_full = as.numeric(exp(rep_len(off, n) + Xc %*% beta)), y = as.integer(Y),
                 xref = .fe_dummy_ref(colMeans(Xc), Xc, Y, fe, .ident_count),
                 zref = .fe_dummy_ref(colMeans(Zz), Zz, Y, zero_fe, .ident_binary),
                 offset = off, offset_name = .ud_offset_name(offset),
                 mu = as.numeric(exp(rep_len(off, n) + Xc %*% beta)),
                 fitted.values = as.numeric((1 - linkinv(Zz %*% gamma)) * fam$meanfun(exp(rep_len(off, n) + Xc %*% beta), theta)),
                 Xc = Xc, Zz = Zz, fe = fe, zero_fe = zero_fe,
                 se.type = se, converged = o$convergence == 0, formula = cform,
                 zero.formula = zrhs, zero_terms = stats::terms(zrhs), levels = .cpb_xlevels(mf),
                 contrasts = attr(Xc, "contrasts"), zero_contrasts = attr(Zz, "contrasts"),
                 call = match.call()), class = "zi_count")
}

## ---------------------------------------------------------------------------
## S3 methods for count_reg / hurdle_count / zi_count
## ---------------------------------------------------------------------------

#' @export
print.count_reg <- function(x, ...) {
  cat(.count_fam(x$family)$label, "count regression",
      if (x$truncated) "(zero-truncated)" else "", if (!is.null(x$fe)) paste0("[FE: ", paste(x$fe, collapse = ", "), "]"), "\n")
  cf <- x$coefficients
  if (!is.null(x$fe)) {                                 # the absorbed unit dummies are not shown
    keep <- !.is_fe_col(names(cf), x$fe); nd <- sum(!keep); cf <- cf[keep]
  }
  cat("Call:  ", deparse(x$call), "\n\nCoefficients:\n", sep = ""); print(round(cf, 4))
  if (!is.null(x$fe) && nd > 0) cat("(", nd, " fixed-effect dummies not shown; see coef())\n", sep = "")
  sc <- .coef_scale(x)$scale
  if (sc != "log mean")
    cat("\nCoefficients are on the ", sc, " scale, not the log mean; predict(type = 'response') and ame() are on the ",
        "mean scale.\n", sep = "")
  cat(.shape_line(x$family, x$theta))
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}
#' @method coef count_reg
#' @export
coef.count_reg <- function(object, ...) object$coefficients
#' @method vcov count_reg
#' @export
vcov.count_reg <- function(object, ...) {
  if (is.null(object$vcov)) return(NULL)
  object$vcov
}
#' @method logLik count_reg
#' @export
logLik.count_reg <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")
#' @method nobs count_reg
#' @export
nobs.count_reg <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted
#' @method fitted count_reg
#' @export
fitted.count_reg <- function(object, ...) object$fitted.values          # the mean scale (E(Y | Y >= 1) when truncated)
#' @rdname residuals.underdisp
#' @method residuals count_reg
#' @export
residuals.count_reg <- function(object, type = c("response", "pearson"), ...) .ud_residuals(object, match.arg(type))
#' Predictions from a matched count-family fit
#'
#' @param object A `"count_reg"` object.
#' @param newdata Optional data frame of covariate profiles.
#' @param type `"response"` (the mean E(Y); for a zero-truncated fit the
#'   conditional mean E(Y | Y > 0)), `"link"` (the linear predictor), or
#'   `"prob"` (P(Y = `at`); for a zero-truncated fit this is P(Y = `at` | Y > 0)).
#' @param at Count value for `type = "prob"`.
#' @param offset Optional offset (log scale) for `newdata`: a numeric vector or a
#'   column name in `newdata`.
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict count_reg
#' @export
predict.count_reg <- function(object, newdata = NULL, type = c("response", "link", "prob"),
                              at = NULL, offset = NULL, se.fit = FALSE, interval = c("none", "confidence"),
                              level = 0.95, ...) {
  type <- match.arg(type); fam <- .count_fam(object$family)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset, at = at))
  eta <- if (is.null(newdata)) object$linear.predictors else .count_reg_eta(object, newdata, offset)
  bad <- if (is.null(newdata)) FALSE else is.na(eta)          # a missing covariate or offset predicts NA
  if (any(bad)) eta[bad] <- 0
  mu <- exp(eta)
  if (type == "link") return(.ud_mask(eta, bad))
  if (type == "prob") {
    if (is.null(at)) stop("type = \"prob\" needs 'at' (a non-negative integer count).")
    return(.ud_mask(vapply(mu, function(m) {
      pk <- fam$pvec(m, object$theta, at)[at + 1L]           # P(Y = at)
      if (isTRUE(object$truncated))                          # condition on Y > 0: P(Y=at)/(1-P(0))
        pk <- if (at == 0L) 0 else pk / (1 - fam$p0(m, object$theta))
      pk
    }, numeric(1)), bad))
  }
  m <- fam$meanfun(mu, object$theta)
  if (isTRUE(object$truncated)) m <- m / pmax(1 - fam$p0(mu, object$theta), 1e-12)   # E(Y | Y > 0)
  .ud_mask(m, bad)
}
## The linear predictor of a count_reg() fit in new data, offset included. A unit of an `fe` column the fit has not
## seen (or every row, when `newdata` lacks the column) enters at the reference unit `xref` (by default the fit's own,
## the units its data identify; a hurdle passes its intensity's), and with `all_factors` (held-out scoring) an unseen
## level of any factor; the rows the fit can evaluate are computed as before.
.count_reg_eta <- function(object, newdata, offset = NULL, xref = NULL, all_factors = FALSE) {
  Terms <- stats::delete.response(object$terms)
  miss <- setdiff(all.vars(Terms), c(names(newdata), object$fe))
  if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
  if (!.ud_has_unseen(object$levels, newdata, object$fe, all_factors)) {
    X <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    return(.ud_newdata_offset(object, newdata, offset, nrow(X)) + as.numeric(X %*% object$coefficients))
  }
  if (is.null(xref) && !is.null(object$X))
    xref <- .fe_dummy_ref(colMeans(object$X), object$X, object$Y, object$fe,
                          if (isTRUE(object$truncated)) .ident_trunc else .ident_count)
  .ud_newdata_offset(object, newdata, offset, nrow(newdata)) +
    .ud_dummy_eta(Terms, object$levels, object$coefficients, newdata, object$contrasts, xref, object$fe, all_factors)
}

#' @method summary count_reg
#' @export
summary.count_reg <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(coefficients = ctab, family = object$family, label = .count_fam(object$family)$label,
                 se_ratio = if (!is.null(object$se_model)) object$se.beta / object$se_model,
                 truncated = object$truncated, n = object$n, n_dropped = object$n_dropped, nobs = nobs(object),
                 se.type = object$se.type,
                 theta = object$theta, loglik = object$loglik, df = object$df, converged = object$converged, fe = object$fe),
            class = "summary.count_reg")
}
#' @export
print.summary.count_reg <- function(x, ...) {
  cat("\n", x$label, " count regression", if (x$truncated) " (zero-truncated)" else "", "\n", sep = "")
  cat("N =", x$n, if (!is.null(x$nobs) && x$nobs != x$n) paste0(" (weight total ", format(x$nobs), ")"),
      if (isTRUE(x$n_dropped > 0)) paste0(" (", x$n_dropped, " rows with a missing model variable dropped)"),
      "  inference:", x$se.type, "\n\n")
  ct <- x$coefficients
  if (!is.null(x$fe)) { keep <- !.is_fe_col(rownames(ct), x$fe); nd <- sum(!keep); ct <- ct[keep, , drop = FALSE] }
  stats::printCoefmat(ct, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  if (!is.null(x$fe) && nd > 0) cat("(", nd, " fixed-effect dummies not shown; see coef())\n", sep = "")
  if (!is.null(x$se_ratio)) {
    r <- x$se_ratio[rownames(ct)]
    cat(sprintf("%s / model-based standard errors: %s\n  (%s)\n",
                if (identical(x$se.type, "cluster")) "Cluster-robust" else "Robust",
                paste(sprintf("%s %.2f", names(r), r), collapse = ", "),
                if (identical(x$family, "poisson"))
                  "Poisson: a gap shows the variance is not the Poisson's; the robust errors stay valid for the mean parameters"
                else "King and Roberts 2015: a large gap signals a misspecified model"))
  }
  cat(.shape_line(x$family, x$theta))
  cat("logLik =", round(x$loglik, 2), "  AIC =", round(-2 * x$loglik + 2 * x$df, 2), "\n")
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  invisible(x)
}
#' @method coef zi_count
#' @export
coef.zi_count <- function(object, ...) .zi_coef(object)
#' @method vcov zi_count
#' @export
vcov.zi_count <- function(object, ...) {
  if (is.null(object$vcov_full)) return(NULL)
  k <- length(object$coefficients) + length(object$zero.coefficients)   # count + inflation blocks
  object$vcov_full[seq_len(k), seq_len(k), drop = FALSE]
}
#' @method logLik zi_count
#' @export
logLik.zi_count <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")
#' @method nobs zi_count
#' @export
nobs.zi_count <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted
#' @export
print.zi_count <- function(x, ...) {
  cat("Zero-inflated", .count_fam(x$family)$label, "regression\n")
  cat("Call:  ", deparse(x$call), "\n\nCount coefficients:\n", sep = ""); print(round(x$coefficients, 4))
  cat(sprintf("\nZero-inflation (%s link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n",
              x$link)); print(round(x$zero.coefficients, 4))
  cat(.shape_line(x$family, x$theta))
  if (length(x$fe_zero_units) && identical(sort(x$fe_zero_units), sort(x$zero_fe_zero_units))) {
    .ud_zero_df_note(x$fe_zero_units, "count effect and zero-equation dummy", "count effects and zero-equation dummies")
  } else {
    .ud_zero_df_note(x$fe_zero_units, "count effect", "count effects")
    .ud_zero_df_note(x$zero_fe_zero_units, "zero-equation dummy", "zero-equation dummies")
  }
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
  .zi_separation_note(x)
  invisible(x)
}
#' @method logLik hurdle_count
#' @export
logLik.hurdle_count <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = nobs(object), class = "logLik")
#' @method nobs hurdle_count
#' @export
nobs.hurdle_count <- function(object, ...) if (is.null(object$nobs_weighted)) object$n else object$nobs_weighted
#' @export
print.hurdle_count <- function(x, ...) {
  cat("Hurdle", .count_fam(x$family)$label, "regression")
  if (!is.null(x$fe)) cat(" (intensity FE: ", x$fe, ")", sep = ""); cat("\n")
  cat(sprintf("Units: %d (%d participate, %.0f%%)\n", x$n, x$n_participate,
              100 * x$n_participate / x$n))
  pc <- stats::coef(x$participation); pc <- pc[!.is_fe_col(names(pc), x$part_fe)]
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n",
              x$link)); print(round(pc, 4))
  cat("\nIntensity (zero-truncated ", .count_fam(x$family)$label, "):\n", sep = "")
  print(round(x$intensity$coefficients, 4))
  cat(.shape_line(x$family, x$theta))
  for (nt in .sep_note(x$part_separation, "participation")) cat(nt, "\n")
  .ud_zero_df_note(x$part_zero_units, "participation dummy", "participation dummies")
  invisible(x)
}

#' @rdname irr.count
#' @method irr count_reg
#' @export
irr.count_reg <- function(object, level = 0.95, ...) {
  .irr_drop_fe(.count_irr_rows(object, object$coefficients,
                               if (!is.null(object$vcov)) sqrt(pmax(diag(object$vcov), 0)) else NULL, level), object$fe)
}
