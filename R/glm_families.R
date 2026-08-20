## ---------------------------------------------------------------------------
## Native Poisson and negative-binomial count families.
##
## Matched baselines for the CPB/GEC family: single-equation, fixed-effects,
## zero-truncated, hurdle, and zero-inflated Poisson and negative binomial, all
## fit and scored through one interface with identical df, log-likelihood, proper
## score, and robust/clustered-SE accounting. The value is comparability, not
## novelty: compare_models() can place a cpb_fe next to an nb_fe or a zi_pois on
## exactly the same footing, so model selection happens in one framework.
##
## Design: everything reduces to a per-observation log-likelihood in the mean mu
## (log link) plus an optional shape (negative-binomial theta). A single robust
## sandwich (.count_vcov) turns per-observation scores into analytic, HC-robust,
## or cluster-robust covariance, uniformly across every variant.
## ---------------------------------------------------------------------------

## family objects: pmf on the log scale, pmf vector over 0..kmax, and P(Y=0)
## Each family carries logpmf/pvec/p0 in its natural parameter `mu` (the mean for
## Poisson/NB, the rate lambda for COM-Poisson) plus `meanfun`, which maps that
## parameter to E(Y) so predict()/first_difference() always report the mean.
.count_families <- list(
  poisson = list(tag = "poisson", label = "Poisson", nshape = 0L,
    logpmf = function(y, mu, theta) stats::dpois(y, mu, log = TRUE),
    pvec   = function(mu, theta, kmax) stats::dpois(0:kmax, mu),
    p0     = function(mu, theta) exp(-mu),
    meanfun = function(mu, theta) mu),
  negbin = list(tag = "negbin", label = "NegBinomial", nshape = 1L,
    logpmf = function(y, mu, theta) stats::dnbinom(y, mu = mu, size = theta, log = TRUE),
    pvec   = function(mu, theta, kmax) stats::dnbinom(0:kmax, mu = mu, size = theta),
    p0     = function(mu, theta) stats::dnbinom(0, mu = mu, size = theta),
    meanfun = function(mu, theta) mu),
  compois = list(tag = "compois", label = "COM-Poisson", nshape = 1L,
    logpmf = function(y, mu, theta) y * log(mu) - theta * lgamma(y + 1) - .compois_logZ(mu, theta),
    pvec   = function(mu, theta, kmax) .compois_pmf1(mu, theta, kmax),
    p0     = function(mu, theta) exp(-.compois_logZ(mu, theta)),
    meanfun = function(mu, theta) .compois_mean(mu, theta))
)
.count_fam <- function(family) {
  family <- match.arg(family, c("poisson", "negbin", "compois"))
  .count_families[[family]]
}

## dispersion/shape display line for print/summary: theta for NB, nu for COM-Poisson,
## nothing for Poisson (which has no free shape).
.shape_line <- function(family, theta) {
  if (family == "negbin") sprintf("\ntheta (NB size): %.4f\n", theta)
  else if (family == "compois") sprintf("\nnu (COM-Poisson dispersion; > 1 underdispersed): %.4f\n", theta)
  else ""
}

## per-observation log-likelihood in (beta, log theta) for a count family, with
## optional zero-truncation. par = c(beta, log theta) when the family has a shape.
.count_llik_i <- function(par, X, Y, fam, truncated, offset = 0) {
  p <- ncol(X)
  eta <- offset + as.numeric(X %*% par[seq_len(p)])       # offset enters the log-mean, exp(offset + x'b)
  mu <- exp(eta)
  theta <- if (fam$nshape) exp(par[p + 1L]) else NA_real_
  li <- fam$logpmf(Y, mu, theta)
  if (truncated) li <- li - log1p(-fam$p0(mu, theta))   # condition on Y >= 1
  li
}

## robust/clustered sandwich from a per-observation log-likelihood `ll_i(par)`.
##   type "analytic" -> inverse observed information (model-based)
##   type "robust"   -> HC sandwich  bread %*% meat %*% bread
##   type "cluster"  -> clustered meat (scores summed within cluster)
## Used uniformly by count_reg() and zi_count(), so every variant's SEs share one
## estimator. The hurdle factorizes, so its margins carry their own SEs.
.count_vcov <- function(par, ll_i, type, cluster = NULL, nm) {
  H <- numDeriv::hessian(function(th) sum(ll_i(th)), par)
  bread <- tryCatch(solve(-H), error = function(e) MASS::ginv(-H))
  V <- if (type == "analytic") {
    bread
  } else {
    S <- numDeriv::jacobian(ll_i, par)                       # n x length(par) score matrix
    meat <- if (type == "cluster") {
      if (is.null(cluster)) stop("cluster-robust SE requires a 'cluster'.")
      crossprod(as.matrix(rowsum(S, group = cluster)))       # sum scores within cluster
    } else crossprod(S)                                      # HC (per-observation)
    bread %*% meat %*% bread
  }
  dimnames(V) <- list(nm, nm)
  V
}

## core estimator shared by count_reg() and the intensity of hurdle/zi.
.count_fit <- function(X, Y, fam, truncated, offset = 0, start = NULL) {
  p <- ncol(X); off <- rep_len(offset, nrow(X))
  b0 <- tryCatch({ v <- stats::glm.fit(X, Y, offset = off, family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, p))
  nll <- function(th) { v <- sum(.count_llik_i(th, X, Y, fam, truncated, off)); if (!is.finite(v)) 1e10 else -v }
  ## Shape families get a second dispersion start and keep the better optimum.
  ## For the negative binomial the likelihood can be monotone in theta toward
  ## the Poisson limit (underdispersed data), where a single theta = 1 start can
  ## leave BFGS stranded far from the plateau; start also near that limit. The
  ## COM-Poisson's Poisson point is nu = 1, its default start, but BFGS can
  ## stall there on real data whose optimum lies at nu > 1 (observed on the
  ## zero-truncated peacekeeping fit: profile optimum nu ~ 1.5, +78 log-lik
  ## points above the nu = 1 stall); start also from mild underdispersion.
  starts <- if (!is.null(start)) list(start)
            else if (fam$nshape && identical(fam$tag, "negbin")) list(c(b0, log(1)), c(b0, log(1e4)))
            else if (fam$nshape) list(c(b0, log(1)), c(b0, log(2))) else list(b0)
  o <- NULL
  for (st in starts) {
    oi <- tryCatch(stats::optim(st, nll, method = "BFGS", control = list(maxit = 500, reltol = 1e-10)),
                   error = function(e) NULL)
    if (!is.null(oi) && (is.null(o) || oi$value < o$value)) o <- oi
  }
  if (is.null(o)) stop("count model optimization failed from all starting values.")
  if (o$convergence != 0)
    o <- stats::optim(o$par, nll, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-12))
  beta <- o$par[seq_len(p)]; names(beta) <- colnames(X)
  theta <- if (fam$nshape) unname(exp(o$par[p + 1L])) else NA_real_
  list(par = o$par, beta = beta, theta = theta, loglik = -o$value, offset = off,
       mu = as.numeric(exp(off + X %*% beta)), converged = o$convergence == 0)
}

#' Poisson and negative-binomial count regression (matched CPB baselines)
#'
#' Fits a Poisson or negative-binomial regression with a log link, optionally
#' zero-truncated and/or with unit fixed effects, returning an object that
#' [compare_models()], [score()], and the \pkg{broom} methods treat on the same
#' footing as a [cpb()] or [gec()] fit. This is the matched baseline for the
#' underdispersed models: it exists so a Poisson/NB can be compared to a CPB with
#' identical degrees-of-freedom, log-likelihood, proper-score, and
#' robust-standard-error accounting, rather than reconciled across packages.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param family `"poisson"`, `"negbin"` (negative binomial), or `"compois"`
#'   (Conway--Maxwell--Poisson). For `"compois"` the coefficients are on the
#'   log-rate scale (\eqn{\log\lambda = x'\beta}; the rate \eqn{\lambda} is
#'   \emph{not} the mean, though `predict(type = "response")` and `fitted()`
#'   return the mean), and the shape in `$theta` is the dispersion \eqn{\nu}
#'   (\eqn{\nu > 1} underdispersed, \eqn{\nu = 1} Poisson, \eqn{\nu < 1}
#'   overdispersed).
#' @param truncated Logical; if `TRUE`, fit the zero-truncated form (requires all
#'   `Y >= 1`), the count analogue of the CPB's zero-truncated default.
#' @param fe Optional column name(s) for fixed effects, entered as factor dummies
#'   (so the degrees of freedom count each absorbed intercept, matching
#'   [cpb_fe()]). Pass a vector of two columns for two-way (e.g. unit and time)
#'   fixed effects.
#' @param offset Optional offset entered on the linear-predictor (log) scale --
#'   the log-mean for `"poisson"`/`"negbin"`, the log-rate (\eqn{\log\lambda}) for
#'   `"compois"`. A numeric vector or the name of a column in `data`, e.g.
#'   \eqn{\log(\text{exposure})}
#'   so the model becomes a rate model.
#' @param se Standard errors: `"analytic"` (inverse information, default),
#'   `"robust"` (heteroskedasticity-consistent sandwich), `"cluster"`
#'   (cluster-robust; needs `cluster`), or `"none"`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector aligned to its rows) for `se = "cluster"`.
#' @param ... Unused.
#' @return An object of class `"count_reg"`. The component `$theta` holds the
#'   shape parameter: the negative-binomial size for `"negbin"`, the COM-Poisson
#'   dispersion \eqn{\nu} for `"compois"`, and `NA` for `"poisson"`.
#' @seealso [cpb()], [compare_models()], [hurdle_count()], [zi_count()]
#' @examples
#' set.seed(1); n <- 400; x <- rnorm(n)
#' y <- rpois(n, exp(1 + 0.5 * x))
#' m <- count_reg(y ~ x, data = data.frame(y = y, x = x), family = "poisson")
#' compare_models(pois = m,
#'                nb = count_reg(y ~ x, data = data.frame(y = y, x = x), family = "negbin"))
#' @export
count_reg <- function(formula, data, family = c("poisson", "negbin", "compois"), truncated = FALSE,
                      fe = NULL, offset = NULL, se = c("analytic", "robust", "cluster", "none"),
                      cluster = NULL, ...) {
  se_missing <- missing(se)
  family <- match.arg(family); fam <- .count_fam(family); se <- match.arg(se)
  ## uniform cluster semantics across the family: supplying `cluster` implies
  ## clustered inference unless the user explicitly chose another se type
  if (!is.null(cluster) && se_missing) se <- "cluster"
  cl <- match.call()
  if (!is.null(fe)) formula <- stats::reformulate(c(labels(stats::terms(formula)),
                                                    paste0("factor(", fe, ")")),  # fe may name >1 column (two-way FE)
                                                  response = all.vars(formula)[1L])
  mf <- stats::model.frame(formula, data, na.action = stats::na.omit)
  Y <- stats::model.response(mf); X <- stats::model.matrix(formula, mf)
  n <- length(Y); p <- ncol(X)
  off <- 0                                                   # log-scale exposure offset
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    off <- as.numeric(ov[match(rownames(mf), rownames(data))])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  if (truncated && any(Y < 1)) stop("truncated = TRUE requires all Y >= 1.")

  clid <- NULL
  if (!is.null(cluster) && se != "cluster")
    warning("'cluster' is ignored unless se = \"cluster\".")
  if (se == "cluster") {
    cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    if (is.null(cv)) stop("se = \"cluster\" needs 'cluster' (a column name or vector).")
    clid <- cv[match(rownames(mf), rownames(data))]
    if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
    ng <- length(unique(clid))
    if (ng < 30L)
      warning(sprintf("only %d clusters; the cluster-robust SE has no small-sample correction and may be unreliable.", ng))
  }

  ft <- .count_fit(X, Y, fam, truncated, offset = off)
  df <- p + fam$nshape
  se.beta <- setNames(rep(NA_real_, p), colnames(X)); vcv <- NULL
  if (se != "none") {
    nm <- c(colnames(X), if (fam$nshape) "log(theta)")
    V <- tryCatch(.count_vcov(ft$par, function(th) .count_llik_i(th, X, Y, fam, truncated, off),
                              se, clid, nm), error = function(e) NULL)
    if (!is.null(V)) {
      vcv <- V[seq_len(p), seq_len(p), drop = FALSE]
      se.beta <- sqrt(pmax(diag(vcv), 0))
    }
  }
  structure(list(
    coefficients = ft$beta, theta = ft$theta, family = family, truncated = truncated,
    se.beta = se.beta, vcov = vcv, loglik = ft$loglik, fitted.values = ft$mu,
    linear.predictors = as.numeric(rep_len(off, n) + X %*% ft$beta), residuals = Y - ft$mu,
    offset = off, n = n, df = df, p = p, fe = fe, se.type = se, clustered = !is.null(clid),
    n_clusters = if (!is.null(clid)) length(unique(clid)) else NA_integer_,
    converged = ft$converged, formula = formula, terms = attr(mf, "terms"),
    levels = .cpb_xlevels(mf), contrasts = attr(X, "contrasts"), X = X, Y = Y, call = cl),
    class = "count_reg")
}

## ---------------------------------------------------------------------------
## Hurdle Poisson / negative binomial (participation logit + zero-truncated
## intensity). The likelihood factorizes, so the two margins are fit separately
## and each carries its own standard errors, exactly as hurdle_cpb().
## ---------------------------------------------------------------------------

#' Hurdle Poisson / negative-binomial regression (matched CPB baseline)
#'
#' Fits a participation logit joined to a zero-truncated Poisson or
#' negative-binomial intensity, the count analogue of [hurdle_cpb()]. Returned as
#' a `"hurdle_count"` object that [compare_models()] and [score()] accept, so a
#' hurdle-CPB and a hurdle-NB can be compared on one footing.
#'
#' @param formula Intensity formula (`y ~ x`).
#' @param data A data frame.
#' @param family `"poisson"`, `"negbin"`, or `"compois"` for the intensity; see
#'   [count_reg()] for the COM-Poisson parameterization.
#' @param participation Optional one-sided formula for the participation model;
#'   defaults to the intensity right-hand side.
#' @param fe,part_fe Optional fixed-effects column names for the intensity and
#'   participation models (entered as factor dummies).
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`. Positive participation coefficients raise the probability of
#'   a *positive count* (participation).
#' @param offset Optional offset for the intensity, on the linear-predictor (log)
#'   scale (log-mean for Poisson/NB, log-rate for COM-Poisson): a numeric vector
#'   or the name of a column in `data`.
#' @param se,cluster Standard-error type and optional cluster for the intensity;
#'   see [count_reg()].
#' @param ... Unused.
#' @return An object of class `"hurdle_count"`.
#' @examples
#' set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(0.4 + 0.8 * z)) == 1, rpois(n, exp(1 + 0.3 * x)) + 1L, 0L)
#' hurdle_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
#'              participation = ~ z, se = "none")
#' @seealso [hurdle_cpb()], [zi_count()], [compare_models()]
#' @export
hurdle_count <- function(formula, data, family = c("poisson", "negbin", "compois"),
                         participation = NULL, fe = NULL, part_fe = NULL,
                         link = c("logit", "probit", "cloglog"), offset = NULL,
                         se = c("analytic", "robust", "cluster", "none"), cluster = NULL, ...) {
  family <- match.arg(family); se <- match.arg(se); link <- match.arg(link)
  ## reduce to complete cases on all model variables so the two margins stay aligned
  mv <- unique(c(all.vars(formula), all.vars(if (is.null(participation)) formula[-2L] else participation),
                 fe, part_fe))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  if (!is.null(offset) && !is.character(offset)) offset <- offset[keep]
  data <- data[keep, , drop = FALSE]
  y <- stats::model.response(stats::model.frame(formula, data))
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  d <- as.integer(y > 0)
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, family = stats::binomial(link = link))
  pos <- data[y > 0, , drop = FALSE]
  int_offset <- if (is.null(offset)) NULL                       # offset applies to the intensity
                else if (is.character(offset) && length(offset) == 1L) offset
                else offset[y > 0]
  ifit <- count_reg(formula, data = pos, family = family, truncated = TRUE,
                    fe = fe, offset = int_offset, se = se, cluster = cluster)
  ## per-observation intensity mean for all units (missing FE levels -> reference)
  beta <- ifit$coefficients
  mm <- stats::model.matrix(stats::delete.response(stats::terms(ifit$formula)), data)
  common <- intersect(names(beta), colnames(mm))
  lambda_full <- exp(as.numeric(mm[, common, drop = FALSE] %*% beta[common]))
  llh <- as.numeric(stats::logLik(pfit)) + ifit$loglik
  dfh <- length(stats::coef(pfit)) + ifit$df
  structure(list(participation = pfit, intensity = ifit, y = as.integer(y),
                 p_full = as.numeric(stats::fitted(pfit)), lambda_full = lambda_full,
                 theta = ifit$theta, family = family, loglik = llh, df = dfh, link = link,
                 int_beta = beta, int_xref = colMeans(mm[, common, drop = FALSE]),
                 part_beta = stats::coef(pfit), part_xref = colMeans(stats::model.matrix(pfit)),
                 n = length(y), n_participate = sum(d), fe = fe, part_fe = part_fe,
                 formula = formula, part_formula = participation, call = match.call()),
            class = "hurdle_count")
}

## ---------------------------------------------------------------------------
## Zero-inflated Poisson / negative binomial (structural-zero mixture). Unlike
## the hurdle this does not factorize, so it is a single joint maximization.
## ---------------------------------------------------------------------------

## per-observation log-likelihood in (beta_count, gamma_zero, log theta). The
## structural-zero probability uses an arbitrary link (logit/probit/cloglog).
.zi_llik_i <- function(par, Xc, Zz, Y, fam, linkinv = stats::plogis, offset = 0) {
  pc <- ncol(Xc); pz <- ncol(Zz)
  mu <- exp(offset + as.numeric(Xc %*% par[seq_len(pc)]))
  pistar <- linkinv(as.numeric(Zz %*% par[pc + seq_len(pz)]))
  theta <- if (fam$nshape) exp(par[pc + pz + 1L]) else NA_real_
  f0 <- fam$p0(mu, theta)
  ifelse(Y == 0, log(pistar + (1 - pistar) * f0), log1p(-pistar) + fam$logpmf(Y, mu, theta))
}

#' Zero-inflated Poisson / negative-binomial regression (matched CPB baseline)
#'
#' Fits a structural-zero mixture with a Poisson or negative-binomial count
#' component, the count analogue of [zi_cpb()]. The count component uses a log
#' link (on the mean for Poisson/NB, on the rate \eqn{\lambda} for COM-Poisson)
#' and the inflation probability a link set by `link`. Returned as a `"zi_count"`
#' object that [compare_models()] and [score()] accept.
#'
#' @param formula Count formula (`y ~ x`).
#' @param data A data frame.
#' @param family `"poisson"`, `"negbin"`, or `"compois"` for the count component;
#'   see [count_reg()] for the COM-Poisson parameterization.
#' @param zero Optional one-sided formula for the inflation (structural-zero)
#'   model; defaults to the count right-hand side.
#' @param fe Optional fixed-effects column for the count equation (factor dummies).
#' @param zero_fe Optional fixed-effects column for the inflation equation
#'   (factor dummies); zero-equation fixed effects are opt-in.
#' @param link Link for the inflation probability: `"logit"` (default),
#'   `"probit"`, or `"cloglog"`. Positive inflation coefficients raise the
#'   probability of a *structural zero*.
#' @param offset Optional offset for the count component, on the linear-predictor
#'   (log) scale (log-mean for Poisson/NB, log-rate for COM-Poisson): a numeric
#'   vector or the name of a column in `data`.
#' @param se,cluster Standard-error type and optional cluster; see [count_reg()].
#' @param ... Unused.
#' @return An object of class `"zi_count"`.
#' @examples
#' set.seed(2); n <- 300; x <- rnorm(n); z <- rnorm(n)
#' y <- ifelse(rbinom(n, 1, plogis(-0.5 + 0.8 * z)) == 1, 0L, rpois(n, exp(1 + 0.3 * x)))
#' zi_count(y ~ x, data.frame(y = y, x = x, z = z), family = "poisson",
#'          zero = ~ z, se = "none")
#' @seealso [zi_cpb()], [hurdle_count()], [compare_models()]
#' @export
zi_count <- function(formula, data, family = c("poisson", "negbin", "compois"), zero = NULL,
                     fe = NULL, zero_fe = NULL, link = c("logit", "probit", "cloglog"),
                     offset = NULL, se = c("analytic", "robust", "cluster", "none"), cluster = NULL, ...) {
  se_missing <- missing(se)
  family <- match.arg(family); fam <- .count_fam(family); se <- match.arg(se); link <- match.arg(link)
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
  mf <- stats::model.frame(allform, data, na.action = stats::na.omit)
  Y <- stats::model.response(mf); Xc <- stats::model.matrix(cform, mf); Zz <- stats::model.matrix(zrhs, mf)
  n <- length(Y); pc <- ncol(Xc); pz <- ncol(Zz)
  if (any(Y < 0) || any(Y != floor(Y))) stop("Response must be non-negative integer counts.")
  clid <- NULL
  if (!is.null(cluster) && se != "cluster")
    warning("'cluster' is ignored unless se = \"cluster\".")
  if (se == "cluster") {
    cv <- if (is.character(cluster) && length(cluster) == 1L) data[[cluster]] else cluster
    if (is.null(cv)) stop("se = \"cluster\" needs 'cluster'.")
    clid <- cv[match(rownames(mf), rownames(data))]
    if (anyNA(clid)) stop("'cluster' has missing values on the estimation rows.")
    ng <- length(unique(clid))
    if (ng < 30L) warning(sprintf("only %d clusters; the cluster-robust SE may be unreliable.", ng))
  }
  off <- 0                                                     # offset enters the count component
  if (!is.null(offset)) {
    ov <- if (is.character(offset) && length(offset) == 1L) data[[offset]] else offset
    off <- as.numeric(ov[match(rownames(mf), rownames(data))])
    if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  }
  b0 <- tryCatch({ v <- stats::glm.fit(Xc, Y, offset = rep_len(off, n), family = stats::poisson())$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, pc))
  g0 <- tryCatch({ v <- stats::glm.fit(Zz, as.integer(Y == 0), family = stats::binomial(link = link))$coefficients
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, pz))
  ## pscl-style seeding: also size the count component to the POSITIVE counts via
  ## a zero-truncated fit, so heavy inflation cannot drag the starting rate toward
  ## a poor basin; keep the all-data Poisson start as a fallback and take the
  ## better optimum.
  st_tr <- tryCatch({
    pos <- Y > 0
    ft  <- .count_fit(Xc[pos, , drop = FALSE], Y[pos], fam, truncated = TRUE,
                      offset = rep_len(off, n)[pos])
    c(ft$beta, g0, if (fam$nshape) log(max(ft$theta, 1e-3)))
  }, error = function(e) NULL)
  nll <- function(th) { v <- sum(.zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off)); if (!is.finite(v)) 1e10 else -v }
  o <- NULL
  for (st in Filter(Negate(is.null), list(st_tr, c(b0, g0, if (fam$nshape) log(1))))) {
    oi <- tryCatch(stats::optim(st, nll, method = "BFGS", control = list(maxit = 1000, reltol = 1e-10)),
                   error = function(e) NULL)
    if (!is.null(oi) && (is.null(o) || oi$value < o$value)) o <- oi
  }
  if (is.null(o)) stop("zi_count: optimization failed from all starting values.")
  beta <- o$par[seq_len(pc)]; names(beta) <- colnames(Xc)
  gamma <- o$par[pc + seq_len(pz)]; names(gamma) <- colnames(Zz)
  theta <- if (fam$nshape) unname(exp(o$par[pc + pz + 1L])) else NA_real_
  df <- pc + pz + fam$nshape
  se.beta <- setNames(rep(NA_real_, pc), colnames(Xc)); vcv <- NULL
  if (se != "none") {
    nm <- c(colnames(Xc), paste0("zero_", colnames(Zz)), if (fam$nshape) "log(theta)")
    V <- tryCatch(.count_vcov(o$par, function(th) .zi_llik_i(th, Xc, Zz, Y, fam, linkinv, off), se, clid, nm),
                  error = function(e) NULL)
    if (!is.null(V)) { vcv <- V[seq_len(pc), seq_len(pc), drop = FALSE]; se.beta <- sqrt(pmax(diag(vcv), 0)) }
  }
  structure(list(coefficients = beta, zero.coefficients = gamma, zero_coef = gamma, theta = theta, family = family,
                 se.beta = se.beta, vcov = vcv, loglik = -o$value, df = df, n = n, link = link,
                 pi_full = as.numeric(linkinv(Zz %*% gamma)),
                 lambda_full = as.numeric(exp(rep_len(off, n) + Xc %*% beta)), y = as.integer(Y),
                 xref = colMeans(Xc), zref = colMeans(Zz), offset = off,
                 fitted.values = as.numeric(exp(rep_len(off, n) + Xc %*% beta)), Xc = Xc, Zz = Zz, fe = fe, zero_fe = zero_fe,
                 se.type = se, converged = o$convergence == 0, formula = cform,
                 zero.formula = zrhs, call = match.call()), class = "zi_count")
}

## ---------------------------------------------------------------------------
## S3 methods for count_reg / hurdle_count / zi_count
## ---------------------------------------------------------------------------

#' @export
print.count_reg <- function(x, ...) {
  cat(.count_fam(x$family)$label, "count regression",
      if (x$truncated) "(zero-truncated)" else "", if (!is.null(x$fe)) paste0("[FE: ", x$fe, "]"), "\n")
  cat("Call:  ", deparse(x$call), "\n\nCoefficients:\n", sep = ""); print(round(x$coefficients, 4))
  if (x$family == "compois")
    cat("\nCoefficients are on the log-rate scale; predict(type='response') is on the mean scale.\n")
  cat(.shape_line(x$family, x$theta))
  invisible(x)
}
#' @method coef count_reg
#' @export
coef.count_reg <- function(object, ...) object$coefficients
#' @method vcov count_reg
#' @export
vcov.count_reg <- function(object, ...) {
  if (is.null(object$vcov)) stop("No covariance; refit with se != \"none\".")
  object$vcov
}
#' @method logLik count_reg
#' @export
logLik.count_reg <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
#' @method nobs count_reg
#' @export
nobs.count_reg <- function(object, ...) object$n
#' @method fitted count_reg
#' @export
fitted.count_reg <- function(object, ...)
  .count_fam(object$family)$meanfun(object$fitted.values, object$theta)  # mean scale (identity for pois/nb)
#' @method residuals count_reg
#' @export
residuals.count_reg <- function(object, ...) object$Y - fitted(object)
#' Predictions from a matched Poisson/NB/COM-Poisson fit
#'
#' @param object A `"count_reg"` object.
#' @param newdata Optional data frame of covariate profiles.
#' @param type `"response"` (the mean E(Y)), `"link"` (the linear predictor), or
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
    X <- stats::model.matrix(stats::delete.response(object$terms), newdata, xlev = object$levels)
    eta <- as.numeric(X[, names(object$coefficients), drop = FALSE] %*% object$coefficients)
    if (!is.null(offset)) {                                   # offset for newdata (log scale)
      ov <- if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset
      eta <- eta + as.numeric(ov)
    }
  }
  mu <- exp(eta)
  if (type == "link") return(eta)
  if (type == "prob") {
    if (is.null(at)) stop("type = \"prob\" needs 'at' (a non-negative integer count).")
    return(vapply(mu, function(m) {
      pk <- fam$pvec(m, object$theta, at)[at + 1L]           # P(Y = at)
      if (isTRUE(object$truncated))                          # condition on Y > 0: P(Y=at)/(1-P(0))
        pk <- if (at == 0L) 0 else pk / (1 - fam$p0(m, object$theta))
      pk
    }, numeric(1)))
  }
  fam$meanfun(mu, object$theta)
}
#' @method summary count_reg
#' @export
summary.count_reg <- function(object, ...) {
  z <- object$coefficients / object$se.beta
  ctab <- cbind(Estimate = object$coefficients, `Std. Error` = object$se.beta,
                `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  structure(list(coefficients = ctab, family = object$family, label = .count_fam(object$family)$label,
                 truncated = object$truncated, n = object$n, se.type = object$se.type,
                 theta = object$theta, loglik = object$loglik, df = object$df),
            class = "summary.count_reg")
}
#' @export
print.summary.count_reg <- function(x, ...) {
  cat("\n", x$label, " count regression", if (x$truncated) " (zero-truncated)" else "", "\n", sep = "")
  cat("N =", x$n, "  inference:", x$se.type, "\n\n")
  stats::printCoefmat(x$coefficients, P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  cat(.shape_line(x$family, x$theta))
  cat("logLik =", round(x$loglik, 2), "  AIC =", round(-2 * x$loglik + 2 * x$df, 2), "\n")
  invisible(x)
}
#' @method coef zi_count
#' @export
coef.zi_count <- function(object, ...) .zi_coef(object)
#' @method vcov zi_count
#' @export
vcov.zi_count <- function(object, ...) {
  if (is.null(object$vcov)) stop("No covariance; refit with se != \"none\".")
  object$vcov
}
#' @method logLik zi_count
#' @export
logLik.zi_count <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
#' @method nobs zi_count
#' @export
nobs.zi_count <- function(object, ...) object$n
#' @export
print.zi_count <- function(x, ...) {
  cat("Zero-inflated", .count_fam(x$family)$label, "regression\n")
  cat("Call:  ", deparse(x$call), "\n\nCount coefficients:\n", sep = ""); print(round(x$coefficients, 4))
  cat(sprintf("\nZero-inflation (%s link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n",
              x$link)); print(round(x$zero.coefficients, 4))
  cat(.shape_line(x$family, x$theta))
  invisible(x)
}
#' @method logLik hurdle_count
#' @export
logLik.hurdle_count <- function(object, ...)
  structure(object$loglik, df = object$df, nobs = object$n, class = "logLik")
#' @method nobs hurdle_count
#' @export
nobs.hurdle_count <- function(object, ...) object$n
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
