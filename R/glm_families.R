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
  theta <- if (fam$nshape) fam$shape_inv(par[p + 1L]) else NA_real_
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
.count_vcov <- function(par, ll_i, type, cluster = NULL, nm, w = NULL, scale = NULL) {
  n <- length(ll_i(par)); w <- .ud_w1(w, n)
  k <- length(par); sc <- if (is.null(scale)) rep(1, k) else scale
  ps <- par * sc; to_par <- function(v) v / sc
  ## differenced at an absolute step of 1e-3 in the scaled parameterization
  ## (numDeriv's default is a step of a tenth of the parameter, which at large
  ## fitted means leaves the family's support and returns NA)
  ma <- list(eps = 1e-3, d = 0.1, zero.tol = 1e-8, r = 4, v = 2, show.details = FALSE)
  H <- numDeriv::hessian(function(u) sum(w * ll_i(to_par(ps + u))), rep(0, k), method.args = ma)
  bread <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  V <- if (type == "analytic") {
    bread
  } else {
    S <- numDeriv::jacobian(function(u) ll_i(to_par(ps + u)), rep(0, k), method.args = ma)   # n x length(par) scores
    meat <- if (type == "cluster") {
      if (is.null(cluster)) stop("cluster-robust SE requires a 'cluster'.")
      G <- length(unique(cluster))
      (G / (G - 1)) * crossprod(as.matrix(rowsum(S * w, group = cluster)))   # sum weighted scores within cluster
    } else crossprod(S * sqrt(w))                            # HC (per-observation)
    bread %*% meat %*% bread
  }
  V <- V * tcrossprod(1 / sc)                                # back to the original parameterization
  dimnames(V) <- list(nm, nm)
  V
}

## core estimator shared by count_reg() and the intensity of hurdle/zi.
.count_fit <- function(X, Y, fam, truncated, offset = 0, start = NULL, w = NULL) {
  p <- ncol(X); off <- rep_len(offset, nrow(X)); w <- .ud_w1(w, nrow(X))
  s <- .ud_colscale(X, w); Xs <- sweep(X, 2, s, "/")           # optimize on unit-SD (weighted) columns
  sc <- c(s, if (fam$nshape) 1)
  if (!is.null(start)) start <- start * sc
  b0 <- tryCatch({ v <- stats::glm.fit(Xs, Y, weights = w, offset = off, family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  nll <- function(th) { v <- sum(w * .count_llik_i(th, Xs, Y, fam, truncated, off)); if (!is.finite(v)) 1e10 else -v }
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
  ## a Nelder-Mead polish from the BFGS solution (kept only if it improves)
  if (length(o$par) > 1L) {
    o2 <- tryCatch(stats::optim(o$par, nll, method = "Nelder-Mead", control = list(maxit = 2000, reltol = 1e-12)),
                   error = function(e) NULL)
    if (!is.null(o2) && o2$value < o$value) { o2$convergence <- 0L; o <- o2 }
  }
  if (!is.finite(o$value) || o$value >= 1e9) {            # never left the infeasible plateau
    warning("the optimizer found no finite likelihood: a fitted mean lies outside [1e-12, 1e8] at every ",
            "point it visited. The fit did not converge; the coefficients are its last iterate.", call. = FALSE)
    o$value <- Inf; o$convergence <- 1L
  }
  par <- o$par / sc                                          # back to the original scale
  beta <- par[seq_len(p)]; names(beta) <- colnames(X)
  theta <- if (fam$nshape) unname(fam$shape_inv(par[p + 1L])) else NA_real_
  list(par = par, beta = beta, theta = theta, loglik = -o$value, offset = off, scale = sc,
       mu = as.numeric(exp(off + X %*% beta)), converged = o$convergence == 0)
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
    warning(sprintf("only %d clusters; the cluster-robust SE has no small-sample correction and may be unreliable.", ng))
  clid
}

#' Matched count regressions: Poisson, negative binomial, COM-Poisson, generalized Poisson, gamma-count, double Poisson
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
#'   \eqn{\lambda \in (-1, 1)} and mean \eqn{\mu = \exp(x'\beta)}, so
#'   \eqn{\mathrm{Var}/\mathrm{Mean} = 1/(1-\lambda)^2}; \eqn{\lambda < 0} is
#'   underdispersion, on a finite support that the pmf is renormalized over
#'   (see [genpois-distribution]).
#' * `"gammacount"`: Winkelmann's (1995) gamma-count renewal-process model,
#'   underdispersed when the waiting times between events are more regular than
#'   exponential (\eqn{\alpha > 1}); \eqn{\exp(x'\beta)} is the long-run event
#'   rate and the exact mean is reported (see [gammacount-distribution]).
#' * `"doublepois"`: Efron's (1986) double Poisson with the exact normalizing
#'   constant, \eqn{\theta > 1} underdispersed (see [doublepois-distribution]).
#'
#' The Poisson is the equidispersed member of every family, so
#' [dispersion_test()] tests each family's dispersion parameter against it, and
#' [dispersion_profile()] compares the families' implied variance-to-mean curves
#' against the data.
#'
#' Runtime: the Poisson, negative binomial, gamma-count, and double Poisson and
#' generalized Poisson families fit in under a second at a few hundred rows; the
#' COM-Poisson families evaluate a normalizing sum per observation (and the
#' mean parameterization solves a root per observation), so they take seconds
#' at a few hundred rows and minutes at several thousand; `se = "robust"` and
#' `"cluster"` add numerical scores at the same cost per parameter.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param family One of `"poisson"`, `"negbin"`, `"compois"`, `"mpcmp"`,
#'   `"genpois"`, `"gammacount"`, `"doublepois"` (see Families).
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
#'   (cluster-robust; needs `cluster`), or `"none"`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector aligned to its rows) for `se = "cluster"`; supplying it selects
#'   `se = "cluster"` unless `se` is given explicitly.
#' @return An object of class `"count_reg"`. The component `$theta` holds the
#'   dispersion parameter on its natural scale (the negative-binomial size, the
#'   COM-Poisson \eqn{\nu}, the generalized-Poisson \eqn{\lambda}, the
#'   gamma-count \eqn{\alpha}, the double-Poisson \eqn{\theta}; `NA` for the
#'   Poisson); `$fitted.values` and `$residuals` are on the mean scale (the
#'   conditional mean E(Y | Y >= 1) for a zero-truncated fit), as `fitted()` and
#'   `residuals()` return them, and `$mu` is the family's natural parameter
#'   `exp(offset + x'b)`.
#' @references Huang, A. (2017). Mean-parametrized Conway--Maxwell--Poisson
#'   regression models for dispersed counts. \emph{Statistical Modelling},
#'   17(6), 359-380. Winkelmann, R. (1995). Duration dependence and dispersion
#'   in count-data models. \emph{Journal of Business & Economic Statistics},
#'   13(4), 467-474. Efron, B. (1986). Double exponential families and their
#'   use in generalized linear regression. \emph{Journal of the American
#'   Statistical Association}, 81(395), 709-721. Consul, P. C. and Famoye, F.
#'   (1992). Generalized Poisson regression model. \emph{Communications in
#'   Statistics -- Theory and Methods}, 21(1), 89-109.
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
                                                "gammacount", "doublepois"),
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
    formula <- stats::reformulate(c(labels(stats::terms(formula)), paste0("factor(", fe, ")")),  # fe may name >1 column (two-way FE)
                                  response = all.vars(formula)[1L])
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
  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL; vcv_full <- NULL
  if (se != "none") {
    nm <- c(colnames(X), if (fam$nshape) fam$shape_name)
    V <- tryCatch(.count_vcov(ft$par, function(th) .count_llik_i(th, X, Y, fam, truncated, off),
                              se, clid, nm, w = w, scale = ft$scale), error = function(e) NULL)
    if (!is.null(V)) {
      vcv_full <- V
      vcv <- V[seq_len(p), seq_len(p), drop = FALSE]
      se.beta <- .count_se(vcv)
    }
  }
  ## $mu is the family's natural parameter (the rate of the recursion); the
  ## stored fitted values and residuals are on the mean scale, as fitted() returns
  fv <- fam$meanfun(ft$mu, ft$theta)
  if (truncated) fv <- fv / pmax(1 - fam$p0(ft$mu, ft$theta), 1e-12)     # E(Y | Y >= 1)
  structure(list(
    coefficients = ft$beta, theta = ft$theta, family = family, truncated = truncated,
    se.beta = se.beta, vcov = vcv, vcov_full = vcv_full, loglik = ft$loglik, mu = ft$mu, fitted.values = fv,
    linear.predictors = as.numeric(rep_len(off, n) + X %*% ft$beta), residuals = Y - fv,
    offset = off, weights = w, n = n, nobs_weighted = if (is.null(w)) n else sum(w),
    df = df, p = p, fe = fe, se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = ft$converged, formula = formula, terms = attr(mf, "terms"),
    levels = .cpb_xlevels(mf), contrasts = attr(X, "contrasts"), X = X, Y = Y, call = cl),
    class = "count_reg")
}

## standard errors from a numerical covariance: a non-positive diagonal means
## the Hessian is not positive definite (the fit may sit on a feasibility
## boundary, with counts at their ceiling), and such an SE is NA, not 0
.count_se <- function(vcv) {
  d <- diag(vcv); bad <- !is.finite(d) | d <= 0
  if (any(bad)) {
    warning("the numerical Hessian is not positive definite for ", sum(bad), " coefficient(s): their standard ",
            "errors are NA. The fit may sit on the family's feasibility boundary (counts at their ceiling); ",
            "compare a bootstrap, or se = \"none\".", call. = FALSE)
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
#'   participation models (entered as factor dummies).
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`. Positive participation coefficients raise the probability of
#'   a *positive count* (participation).
#' @param offset Optional offset for the intensity, on the linear-predictor (log)
#'   scale: a numeric vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name),
#'   applied to both margins; see [count_reg()].
#' @param se,cluster Standard-error type and optional cluster for the intensity;
#'   see [count_reg()].
#' @return An object of class `"hurdle_count"`.
#' @examples
#' set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(0.4 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
#' hurdle_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
#'              participation = ~ z, se = "none")
#' @seealso [hurdle_cpb()], [zi_count()], [compare_models()]
#' @export
hurdle_count <- function(formula, data, family = c("poisson", "negbin", "compois", "mpcmp", "genpois",
                                                   "gammacount", "doublepois"),
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
  dfh <- length(stats::coef(pfit)) + ifit$df
  structure(list(participation = pfit, intensity = ifit,
                 converged = isTRUE(pfit$converged) && isTRUE(ifit$converged), y = as.integer(y),
                 p_full = as.numeric(stats::fitted(pfit)), lambda_full = lambda_full,
                 theta = ifit$theta, family = family, loglik = llh, df = dfh, link = link,
                 int_beta = beta, int_xref = colMeans(mm[, common, drop = FALSE]),
                 part_beta = stats::coef(pfit), part_xref = colMeans(stats::model.matrix(pfit)),
                 n = length(y), nobs_weighted = if (is.null(w_all)) length(y) else sum(w_all),
                 weights = w_all, n_participate = sum(d), fe = fe, part_fe = part_fe,
                 formula = formula, part_formula = participation, call = match.call()),
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
  theta <- if (fam$nshape) fam$shape_inv(par[pc + pz + 1L]) else NA_real_
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
#'   (factor dummies); zero-equation fixed effects are opt-in.
#' @param link Link for the inflation probability: `"logit"` (default),
#'   `"probit"`, or `"cloglog"`. Positive inflation coefficients raise the
#'   probability of a *structural zero*.
#' @param offset Optional offset for the count component, on the linear-predictor
#'   (log) scale: a numeric vector or the name of a column in `data`.
#' @param weights Optional frequency weights (a numeric vector or a column name);
#'   see [count_reg()].
#' @param se,cluster Standard-error type and optional cluster; see [count_reg()].
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
                                               "gammacount", "doublepois"),
                     zero = NULL, fe = NULL, zero_fe = NULL, link = c("logit", "probit", "cloglog"),
                     offset = NULL, weights = NULL, se = c("analytic", "robust", "cluster", "none"),
                     cluster = NULL) {
  .ud_no_formula_offset(formula, zero)
  se_missing <- missing(se)
  fam <- .count_fam(family); family <- fam$tag; se <- match.arg(se); link <- match.arg(link)
  if (!is.null(cluster) && se_missing) se <- "cluster"   # uniform cluster semantics
  linkinv <- stats::make.link(link)$linkinv
  cform <- formula
  if (!is.null(fe)) cform <- stats::reformulate(c(labels(stats::terms(formula)), paste0("factor(", fe, ")")),
                                                response = all.vars(formula)[1L])
  zrhs <- if (is.null(zero)) formula[-2L] else zero
  if (!is.null(zero_fe))                                        # zero-equation FE is opt-in
    zrhs <- stats::reformulate(c(labels(stats::terms(zrhs)), paste0("factor(", zero_fe, ")")))
  allform <- stats::reformulate(unique(c(labels(stats::terms(cform)), labels(stats::terms(zrhs)))),
                                response = all.vars(formula)[1L])
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
  nll <- function(th) { v <- sum(w1 * .zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off)); if (!is.finite(v)) 1e10 else -v }
  ## the optimizer works on the covariates of both equations scaled to unit
  ## (weighted) standard deviation, and each start is polished by chained
  ## Nelder-Mead restarts; the coefficients are mapped back at the end
  sc <- c(.ud_colscale(Xc, w1), .ud_colscale(Zz, w1), if (fam$nshape) 1)
  nll_s <- function(v) nll(v / sc)
  o <- NULL
  for (st in Filter(Negate(is.null), list(st_tr, c(b0, g0, if (fam$nshape) fam$shape_link(fam$shape_starts[1L]))))) {
    oi <- tryCatch(stats::optim(st * sc, nll_s, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10)),
                   error = function(e) NULL)
    if (!is.null(oi)) oi <- .ud_nm_polish(oi, nll_s, maxit = 5000L, reltol = 1e-10)
    if (!is.null(oi) && (is.null(o) || oi$value < o$value)) o <- oi
  }
  if (is.null(o)) stop("zi_count: optimization failed from all starting values.")
  o$par <- o$par / sc
  beta <- o$par[seq_len(pc)]; names(beta) <- colnames(Xc)
  gamma <- o$par[pc + seq_len(pz)]; names(gamma) <- colnames(Zz)
  theta <- if (fam$nshape) unname(fam$shape_inv(o$par[pc + pz + 1L])) else NA_real_
  df <- pc + pz + fam$nshape
  se.beta <- setNames(rep(NA_real_, pc), colnames(Xc)); se.zero <- setNames(rep(NA_real_, pz), colnames(Zz))
  vcv <- NULL; vcv_full <- NULL
  if (se != "none") {
    nm <- c(paste0("count:", colnames(Xc)), paste0("zero:", colnames(Zz)), if (fam$nshape) fam$shape_name)
    V <- tryCatch(.count_vcov(o$par, function(th) .zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off), se, clid, nm, w = w,
                              scale = sc),
                  error = function(e) NULL)
    if (!is.null(V)) {
      vcv_full <- V
      vcv <- V[seq_len(pc), seq_len(pc), drop = FALSE]; se.beta <- .count_se(vcv)
      se.zero <- .count_se(V[pc + seq_len(pz), pc + seq_len(pz), drop = FALSE])
    }
  }
  structure(list(coefficients = beta, zero.coefficients = gamma, zero_coef = gamma, theta = theta, family = family,
                 se.beta = se.beta, se.zero = se.zero, vcov = vcv, vcov_full = vcv_full,
                 loglik = -o$value, df = df, n = n, nobs_weighted = if (is.null(w)) n else sum(w),
                 weights = w, link = link,
                 pi_full = as.numeric(linkinv(Zz %*% gamma)),
                 lambda_full = as.numeric(exp(rep_len(off, n) + Xc %*% beta)), y = as.integer(Y),
                 xref = colMeans(Xc), zref = colMeans(Zz), offset = off,
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
    keep <- !grepl("^factor\\(", names(cf)); nd <- sum(!keep); cf <- cf[keep]
  }
  cat("Call:  ", deparse(x$call), "\n\nCoefficients:\n", sep = ""); print(round(cf, 4))
  if (!is.null(x$fe) && nd > 0) cat("(", nd, " fixed-effect dummies not shown; see coef())\n", sep = "")
  if (x$family == "compois")
    cat("\nCoefficients are on the log-rate scale; predict(type='response') is on the mean scale.\n")
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
#' @method residuals count_reg
#' @export
residuals.count_reg <- function(object, ...) object$Y - fitted(object)
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
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict count_reg
#' @export
predict.count_reg <- function(object, newdata = NULL, type = c("response", "link", "prob"),
                              at = NULL, offset = NULL, ...) {
  type <- match.arg(type); fam <- .count_fam(object$family)
  if (is.null(newdata)) { eta <- object$linear.predictors }
  else {
    Terms <- stats::delete.response(object$terms)
    miss <- setdiff(all.vars(Terms), names(newdata))
    if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    X <- .ud_newdata_matrix(Terms, newdata, object$levels, object$contrasts, names(object$coefficients))
    eta <- as.numeric(X %*% object$coefficients)
    if (!is.null(offset)) {                                   # offset for newdata (log scale)
      ov <- if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset
      eta <- eta + as.numeric(ov)
    }
  }
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
#' @method summary count_reg
#' @export
summary.count_reg <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(coefficients = ctab, family = object$family, label = .count_fam(object$family)$label,
                 truncated = object$truncated, n = object$n, nobs = nobs(object), se.type = object$se.type,
                 theta = object$theta, loglik = object$loglik, df = object$df, converged = object$converged, fe = object$fe),
            class = "summary.count_reg")
}
#' @export
print.summary.count_reg <- function(x, ...) {
  cat("\n", x$label, " count regression", if (x$truncated) " (zero-truncated)" else "", "\n", sep = "")
  cat("N =", x$n, if (!is.null(x$nobs) && x$nobs != x$n) paste0(" (weight total ", format(x$nobs), ")"),
      "  inference:", x$se.type, "\n\n")
  ct <- x$coefficients
  if (!is.null(x$fe)) { keep <- !grepl("^factor\\(", rownames(ct)); nd <- sum(!keep); ct <- ct[keep, , drop = FALSE] }
  stats::printCoefmat(ct, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  if (!is.null(x$fe) && nd > 0) cat("(", nd, " fixed-effect dummies not shown; see coef())\n", sep = "")
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
  if (isFALSE(x$converged)) cat("Note: the optimizer did not report convergence.\n")
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
  pc <- stats::coef(x$participation); pc <- pc[!grepl("^factor\\(", names(pc))]
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n",
              x$link)); print(round(pc, 4))
  cat("\nIntensity (zero-truncated ", .count_fam(x$family)$label, "):\n", sep = "")
  print(round(x$intensity$coefficients, 4))
  cat(.shape_line(x$family, x$theta))
  invisible(x)
}

#' @rdname irr.count
#' @method irr count_reg
#' @export
irr.count_reg <- function(object, level = 0.95, ...) {
  .ud_irr_wald(object$coefficients,
               if (!is.null(object$vcov)) sqrt(diag(object$vcov)) else NULL,
               level, "count", "IRR")
}
