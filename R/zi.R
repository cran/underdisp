## ---------------------------------------------------------------------------
## Zero-inflated continuous parameter binomial (mixture variant), fit by EM.
## Use when a zero is ambiguous between structural nonparticipation and a
## sampling zero from a low-rate count; for a pure participation process where
## participants have Y >= 1 by construction, prefer hurdle_cpb().
## ---------------------------------------------------------------------------

#' Simulate from the zero-inflated continuous parameter binomial
#'
#' @param n Number of units.
#' @param lambda Intensity mean; scalar or length-`n`.
#' @param alpha CPB shape parameter in (0, 1).
#' @param pi Structural-zero probability; scalar or length-`n`.
#' @return An integer vector: a structural zero with probability `pi`, otherwise
#'   an (untruncated) CPB draw, which may itself be zero.
#' @examples
#' set.seed(1)
#' y <- rzicpb(1000, lambda = 3, alpha = 0.5, pi = 0.3)
#' mean(y == 0)
#' @export
rzicpb <- function(n, lambda, alpha, pi) {
  s <- stats::rbinom(n, 1, pi)
  y <- integer(n); idx <- which(s == 0)
  if (length(idx)) y[idx] <- rcpb(length(idx), rep_len(lambda, n)[idx], alpha, truncated = FALSE)
  y
}

#' Zero-inflated continuous parameter binomial regression
#'
#' Fits the zero-inflated CPB. The model mixes a structural-zero process
#' (a logistic model for the probability `pi` that a unit is a structural zero)
#' with an untruncated CPB for the count, so a zero can arise either structurally
#' or as a sampling zero from the CPB. Unlike [hurdle_cpb()], the likelihood does
#' not factorize, so the parameters are estimated by direct joint maximum
#' likelihood, seeded from separate fits (a zero-truncated CPB on the positives
#' and a logit of the zero indicator) and cross-checked against an
#' expectation-maximization climb of the same observed-data likelihood, keeping
#' the better optimum. An EM-only path (`method = "em"`) is retained for
#' comparison.
#'
#' The inflation equation uses a logit link. For a probit or cloglog inflation
#' link, use [zi_count()] (Poisson/NB/COM-Poisson count), whose `link` argument
#' covers the binary stage.
#'
#' @param formula Intensity (count) model formula (`y ~ x`).
#' @param data A data frame.
#' @param zero Optional one-sided formula (`~ z`) for the zero-inflation model;
#'   defaults to the intensity model's right-hand side.
#' @param fe Optional column name for unit fixed effects in the intensity (count)
#'   component. Under `method = "ml"` the unit effects enter as plug-in offsets
#'   from a first-stage fit; under `method = "em"` a classification
#'   (hard-assignment) step handles the structural zeros.
#' @param zero_fe Optional column name for fixed effects in the zero-inflation
#'   equation, entered as factor dummies (opt-in).
#' @param method Estimation method: `"ml"` (default) is robust direct joint
#'   maximum likelihood seeded from separate fits; `"em"` is expectation-
#'   maximization, retained for comparison but prone to a degenerate
#'   structural-zero-probability collapse under heavy zero-inflation.
#' @param se Coefficient inference: `"none"` (default) or `"bootstrap"`.
#' @param B Bootstrap resamples when `se = "bootstrap"`.
#' @param cluster Optional cluster for the bootstrap (a column name or vector);
#'   under fixed effects the units are resampled by default. The mixture EM makes
#'   the bootstrap costly, so keep `B` modest.
#' @param max.support Maximum support for the CPB pmf.
#' @param maxit Optimizer iteration budget (scaled internally for the joint
#'   maximization; also the EM iteration cap under `method = "em"`).
#' @param offset Optional exposure offset (log scale) for the count component:
#'   a numeric vector or the name of a column in `data`, making the count a
#'   rate model. Supported for `method = "ml"` without `fe` (the fixed-effects
#'   path's plug-in unit effects have no offset handling yet, and errors
#'   loudly). The bootstrap carries the offset through every replicate.
#' @param tol Relative convergence tolerance on the observed-data log-likelihood.
#' @return An object of class `"zi_cpb"`.
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 800; x <- rnorm(n); z <- rnorm(n)
#' y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
#' zi_cpb(y ~ x, data = data.frame(y = y, x = x, z = z), zero = ~ z)
#' }
#' @export
zi_cpb <- function(formula, data, zero = NULL, fe = NULL, zero_fe = NULL, method = c("ml", "em"),
                   se = c("none", "bootstrap"), B = 500, cluster = NULL, max.support = 500,
                   maxit = 200, tol = 1e-6, offset = NULL) {
  se <- match.arg(se); method <- match.arg(method)
  if (!is.null(offset)) {
    if (!is.null(fe))
      stop("'offset' with 'fe' is not yet supported for zi_cpb (the fixed effects enter as ",
           "plug-in offsets of their own); fit without 'fe' or fold the exposure into the data.")
    if (identical(method, "em"))
      stop("'offset' is supported for method = \"ml\" only.")
    if (!(is.character(offset) && length(offset) == 1L)) {
      data[[".zi_off"]] <- as.numeric(offset); offset <- ".zi_off"   # vector -> column, so it subsets with data
    }
  }
  res <- .zi_cpb_fit(formula, data, zero, fe, zero_fe, max.support, maxit, tol, method, offset)
  if (se == "bootstrap")
    res <- .zi_cpb_boot(res, formula, data, zero, fe, zero_fe, cluster, B, max.support, maxit, tol,
                        method, offset)
  res
}

## cluster/pairs bootstrap for the zero-inflated CPB (resample clusters, or the
## fixed-effects units by default; each replicate refits the EM).
.zi_cpb_boot <- function(res, formula, data, zero, fe, zero_fe, cluster, B, max.support, maxit, tol,
                         method = "ml", offset = NULL) {
  use_fe <- !is.null(fe)
  cvar <- if (!is.null(cluster)) cluster else if (use_fe) fe else NULL
  if (!is.null(cvar)) {
    cval <- if (is.character(cvar) && length(cvar) == 1L) data[[cvar]] else cvar
    grp  <- split(seq_len(nrow(data)), cval)
    fev  <- if (use_fe) as.character(data[[fe]]) else NULL
  }
  nb <- length(res$coefficients); nz <- length(res$zero_coef)
  Bc <- matrix(NA_real_, B, nb); Zc <- matrix(NA_real_, B, nz)
  for (b in seq_len(B)) {
    if (is.null(cvar)) { bd <- data[sample.int(nrow(data), nrow(data), replace = TRUE), , drop = FALSE]; ffe <- fe }
    else {
      gs <- sample.int(length(grp), length(grp), replace = TRUE)
      rows <- unlist(grp[gs], use.names = FALSE); bd <- data[rows, , drop = FALSE]
      if (use_fe) {
        bd[[".bootunit"]] <- unlist(lapply(seq_along(gs),
          function(k) paste0(k, "_", fev[grp[[gs[k]]]])), use.names = FALSE); ffe <- ".bootunit"
      } else ffe <- NULL
    }
    fb <- tryCatch(.zi_cpb_fit(formula, bd, zero, ffe, zero_fe, max.support, maxit, tol, method, offset),
                   error = function(e) NULL)
    if (!is.null(fb)) { Bc[b, ] <- fb$coefficients; Zc[b, ] <- fb$zero_coef }
  }
  okB <- Bc[stats::complete.cases(Bc), , drop = FALSE]
  okZ <- Zc[stats::complete.cases(Zc), , drop = FALSE]
  if (nrow(okB) >= 2) {
    z <- stats::qnorm(0.975)
    res$se.beta <- setNames(apply(okB, 2, stats::sd), names(res$coefficients))
    res$se.zero <- setNames(apply(okZ, 2, stats::sd), names(res$zero_coef))
    res$ci.beta <- cbind(lower = res$coefficients - z * res$se.beta,
                         upper = res$coefficients + z * res$se.beta)
  } else warning("Too few bootstrap resamples converged for stable inference.")
  res$se.type <- "bootstrap"
  res$cluster <- if (!is.null(cvar)) (if (is.character(cvar) && length(cvar) == 1L) cvar else "custom") else "obs"
  res
}

.zi_cpb_fit <- function(formula, data, zero = NULL, fe = NULL, zero_fe = NULL, max.support = 500, maxit = 200,
                        tol = 1e-6, method = "ml", offset = NULL) {
  mf <- stats::model.frame(formula, data)
  y  <- as.integer(stats::model.response(mf))
  if (any(y < 0)) stop("The response must be nonnegative integer counts.")
  X  <- stats::model.matrix(formula, data)
  zrhs <- if (is.null(zero)) formula[-2L] else zero
  if (!is.null(zero_fe))                                          # fixed effects in the zero equation
    zrhs <- stats::reformulate(c(labels(stats::terms(zrhs)), paste0("factor(", zero_fe, ")")))
  Zt <- stats::terms(stats::update(zrhs, ~ .))
  Z  <- stats::model.matrix(Zt, data)
  int_terms0 <- stats::delete.response(stats::terms(formula))            # stored factor levels so
  int_xlev  <- stats::.getXlevels(int_terms0, mf)                        # predict(newdata=) does not
  zero_xlev <- stats::.getXlevels(Zt, stats::model.frame(Zt, data))      # drop levels absent from newdata
  n  <- length(y); kmax <- max(y); is0 <- y == 0

  ## fixed-effects path. The default ("ml") is a fast TWO-STEP estimator: the unit
  ## count intercepts are taken from a first-stage zero-truncated FE-CPB fit to the
  ## positive counts and held fixed as an offset, after which the observed-data
  ## mixture likelihood is maximized over (beta, alpha, gamma). Because the offset
  ## is estimated from the positives alone it does not use the zeros, and the
  ## inflation equation carries only the terms in `zero` (a global pi unless
  ## `zero_fe` adds unit dummies, which the outer optimizer then fits directly).
  ## method = "em" is a seeded classification-EM that instead re-estimates the unit
  ## count effects from the FULL data (positives and sampling zeros) and fits the
  ## inflation fixed effects (`zero_fe`) cheaply via glm at each step; seeding it
  ## from the separate-stage fits avoids the degenerate pi -> 0 basin that an
  ## unseeded EM can fall into on heavily-zero panels. Use "em" when the inflation
  ## equation needs its own unit fixed effects (a comparison symmetric to a hurdle
  ## with participation fixed effects).
  if (!is.null(fe)) {
    ord  <- order(data[[fe]]); data <- data[ord, , drop = FALSE]
    y    <- y[ord]; is0 <- is0[ord]; Z <- Z[ord, , drop = FALSE]
    Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
    Xcov <- Xall[, setdiff(colnames(Xall), "(Intercept)"), drop = FALSE]
    uf   <- factor(data[[fe]]); nu <- nlevels(uf)
    ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
    pcov <- ncol(Xcov); pg <- ncol(Z); ms <- as.integer(max.support); iit <- 20L
    cf <- tryCatch(cpb_fe(formula, data = data[y > 0, , drop = FALSE], fe = fe, truncated = TRUE,
                          max.support = max.support), error = function(e) NULL)
    b0 <- if (!is.null(cf)) cf$coefficients[colnames(Xcov)] else stats::setNames(rep(0, pcov), colnames(Xcov))
    b0[!is.finite(b0)] <- 0
    a0 <- if (!is.null(cf)) cf$alpha else 0.5
    g0 <- suppressWarnings(stats::glm.fit(Z, as.numeric(is0), family = stats::binomial())$coefficients)
    g0[!is.finite(g0)] <- 0

    if (method == "em") {
      b <- b0; a <- a0; g <- g0; fe_hat <- if (!is.null(cf)) cf$fe else stats::setNames(rep(0, nu), levels(uf))
      lam_fun <- function() { fi <- fe_hat[as.character(data[[fe]])]; fi[is.na(fi)] <- mean(fe_hat)
        as.numeric(exp(fi + Xcov %*% b)) }
      p0_fun <- function(lam) vapply(lam, function(l) .cpb_pmf1(l, a, kmax, FALSE)[1L], numeric(1))
      ll_old <- -Inf; s_old <- rep(NA, length(y)); best <- NULL; iters <- 0L
      for (it in seq_len(maxit)) {
        iters <- it; lam <- lam_fun(); pit <- as.numeric(stats::plogis(Z %*% g))
        w <- ifelse(is0, pit / (pit + (1 - pit) * p0_fun(lam)), 0)
        g <- suppressWarnings(stats::glm.fit(Z, w, family = stats::quasibinomial()))$coefficients
        g[!is.finite(g)] <- 0
        s <- is0 & (w > 0.5); if (it > 1L && all(s == s_old)) break; s_old <- s
        cf2 <- tryCatch(cpb_fe(formula, data = data[!s, , drop = FALSE], fe = fe, truncated = FALSE,
                               max.support = max.support, maxit = 800L, inner_it = 20L), error = function(e) NULL)
        if (!is.null(cf2)) { b <- cf2$coefficients[colnames(Xcov)]; b[!is.finite(b)] <- 0; a <- cf2$alpha; fe_hat <- cf2$fe }
        lam <- lam_fun(); pit <- as.numeric(stats::plogis(Z %*% g))
        py <- vapply(seq_along(y), function(i) .cpb_pmf1(lam[i], a, kmax, FALSE)[y[i] + 1L], numeric(1))
        ## Both mixture branches floored at 1e-12 (the score() convention): an early
        ## iterate can put pit ~ 0 and P(0) ~ 0 on a zero observation (feasibility
        ## leaves the count component no zero mass), and an unfloored log() then
        ## sends ll to -Inf, which turns the convergence test into if (NaN).
        ll <- sum(ifelse(is0, log(pmax(pit + (1 - pit) * p0_fun(lam), 1e-12)),
                         log(pmax(1 - pit, 1e-12)) + log(pmax(py, 1e-12))))
        if (is.finite(ll) && (is.null(best) || ll > best$ll)) best <- list(ll = ll, b = b, a = a, fe_hat = fe_hat, g = g)
        if (is.finite(ll) && is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) break
        ll_old <- ll
      }
      if (is.null(best)) stop("zi_cpb (EM): the log-likelihood was non-finite at every iteration; ",
                              "the count component may not cover the data (try method = \"ml\").")
      b <- best$b; a <- best$a; fe_hat <- best$fe_hat; g <- best$g; ll <- best$ll
      ## Units whose rows are all classified structural (e.g. units with no
      ## positive counts) never enter the count fit, so fe_hat covers only the
      ## units it was estimated on. Expand to all levels with the mean-effect
      ## fallback lam_fun used throughout the iterations, and count only the
      ## genuinely estimated unit effects toward the degrees of freedom.
      fe_df <- length(fe_hat)
      fe_full <- stats::setNames(rep(NA_real_, nu), levels(uf))
      fe_full[names(fe_hat)] <- fe_hat
      fe_full[is.na(fe_full)] <- mean(fe_hat)
      fe_hat <- fe_full
    } else {
      ## Two-step robust ML. PLUG IN the unit count intercepts from a first-stage
      ## zero-truncated FE-CPB fit to the positives (held fixed as an offset -- this
      ## is a plug-in, not a within-mixture concentration, so the offset does not use
      ## the zeros), then jointly maximize the observed-data mixture likelihood over
      ## (beta, alpha, gamma) with that offset fixed. The likelihood is the vectorized
      ## non-FE mixture nll (cpb_lp0_cpp with the offset folded in as a unit-coefficient
      ## column), so this is fast. For an estimator whose unit effects use the zeros,
      ## or for per-unit inflation fixed effects, use method = "em".
      off_lvl <- if (!is.null(cf)) cf$fe[levels(uf)] else stats::setNames(rep(0, nu), levels(uf))
      off_lvl[!is.finite(off_lvl)] <- NA
      off_lvl[is.na(off_lvl)] <- min(off_lvl, na.rm = TRUE) - 2   # all-zero / unfit units -> ~zero rate
      Xoff <- cbind(offset = off_lvl[as.integer(uf)], Xcov)       # offset column, coefficient fixed at 1
      nll <- function(par) {
        m   <- cpb_lp0_cpp(c(1, par[1:pcov], par[pcov + 1]), Xoff, y, rep(0, nrow(Xoff)), ms, FALSE)
        pit <- as.numeric(stats::plogis(Z %*% par[(pcov + 2):(pcov + 1 + pg)]))
        v <- -sum(ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
        if (is.finite(v)) v else 1e10
      }
      best <- NULL
      for (as0 in unique(c(a0, 0.3, 0.6))) {
        op <- tryCatch(stats::optim(c(b0, stats::qlogis(min(max(as0, 0.05), 0.95)), g0), nll,
                                    method = "Nelder-Mead", control = list(maxit = 20L * maxit, reltol = tol)),
                       error = function(e) NULL)
        if (!is.null(op) && op$value < 1e9 && (is.null(best) || op$value < best$value)) best <- op
      }
      if (is.null(best)) stop("zi_cpb (fe): all optimization starts failed.")
      b <- best$par[1:pcov]; a <- stats::plogis(best$par[pcov + 1]); g <- best$par[(pcov + 2):(pcov + 1 + pg)]
      ll <- -best$value; iters <- as.integer(best$counts[1]); fe_hat <- off_lvl
      fe_df <- length(fe_hat)                       # every offset level enters the plug-in likelihood
    }
    names(b) <- colnames(Xcov); names(g) <- colnames(Z); names(fe_hat) <- levels(uf)
    fi <- fe_hat[as.character(data[[fe]])]; fi[is.na(fi)] <- mean(fe_hat)
    lambda_full <- as.numeric(exp(fi + Xcov %*% b))
    return(structure(list(coefficients = b, alpha = a, zero_coef = g, zero.coefficients = g, loglik = ll, iterations = iters,
        df = length(b) + fe_df + 1L + length(g), n = length(y), y = y,
        pi_full = as.numeric(stats::plogis(Z %*% g)), lambda_full = lambda_full,
        ceiling = lambda_full / (1 - a), fe = fe, fe_hat = fe_hat, Z = Z,
        int_beta = b, int_xref = colMeans(Xcov), int_fe_ref = mean(fe_hat),
        zero_beta = g, zero_xref = colMeans(Z),
        int_terms = int_terms0, zero_terms = Zt, int_xlev = int_xlev, zero_xlev = zero_xlev,
        formula = formula, zero_formula = zero, call = match.call()), class = "zi_cpb"))
  }

  ## per-observation log P(Y = y) and P(Y = 0) for the untruncated CPB (C++)
  ms <- as.integer(max.support)
  cpb_lp <- function(beta, alpha) {
    ## `off` resolves from the enclosing frame at call time (defined below,
    ## before any call site)
    m <- cpb_lp0_cpp(c(beta, stats::qlogis(alpha)), X, y, off, ms, FALSE)
    list(lp = pmax(m[, 1], -700), p0 = m[, 2])
  }

  ## Robust estimation: good starting values from separate fits (a zero-truncated
  ## CPB on the positive counts for (beta, alpha), and a logit of the zero
  ## indicator for gamma), then DIRECT joint maximum likelihood over all three,
  ## from several dispersion starts. Starting the count model at its positives-only
  ## fit keeps lambda sized to the positives and the structural-zero probability
  ## away from the degenerate pi -> 0 optimum that coordinate-wise EM falls into
  ## under heavy zero-inflation. (This mirrors how pscl::zeroinfl seeds ZIP/ZINB.)
  pb <- ncol(X); pg <- ncol(Z)
  ## exposure offset (log scale), column-name form only by the time we get here
  off <- if (is.null(offset)) rep(0, nrow(X))
         else as.numeric(data[[offset]][match(rownames(stats::model.frame(formula, data)),
                                              rownames(data))])
  if (anyNA(off)) stop("'offset' has missing values on the estimation rows.")
  pf <- tryCatch(cpb(formula, data[y > 0, , drop = FALSE], truncated = TRUE, se = "none",
                     max.support = max.support,
                     offset = if (is.null(offset)) NULL else offset),
                 error = function(e) NULL)
  b0 <- if (!is.null(pf)) as.numeric(pf$coefficients)
        else { v <- stats::glm.fit(X, y, offset = off, family = stats::poisson())$coefficients
               v[!is.finite(v)] <- 0; v }
  a0 <- if (!is.null(pf)) pf$alpha else 0.5
  g0 <- suppressWarnings(stats::glm.fit(Z, as.numeric(is0), family = stats::binomial())$coefficients)
  g0[!is.finite(g0)] <- 0
  nll <- function(par) {
    m   <- cpb_lp0_cpp(par[1:(pb + 1)], X, y, off, ms, FALSE)
    pit <- as.numeric(stats::plogis(Z %*% par[(pb + 2):(pb + 1 + pg)]))
    -sum(ifelse(is0, log(pit + (1 - pit) * m[, 2]), log(1 - pit) + pmax(m[, 1], -700)))
  }
  run_em <- function() {
    ## Expectation-maximization from the same seeds. The default "ml" branch also
    ## calls this as a cross-check: both climb the same observed-data likelihood
    ## by different routes, and on some data each escapes a basin that strands
    ## the other.
    b <- b0; a <- a0; g <- g0; ll_old <- -Inf; bestem <- NULL
    for (it in seq_len(maxit)) {
      pit <- as.numeric(stats::plogis(Z %*% g)); cl <- cpb_lp(b, a)
      w <- ifelse(is0, pit / (pit + (1 - pit) * cl$p0), 0)                # E-step
      g <- suppressWarnings(stats::glm.fit(Z, w, family = stats::quasibinomial()))$coefficients
      g[!is.finite(g)] <- 0                                              # M-step: zero-inflation
      op <- stats::optim(c(b, stats::qlogis(a)),                         # M-step: weighted CPB
                         function(par) cpb_wnll_cpp(par, X, y, 1 - w, rep(0, nrow(X)), ms, FALSE),
                         method = "Nelder-Mead", control = list(maxit = 400, reltol = 1e-7))
      b <- op$par[-length(op$par)]; a <- stats::plogis(op$par[length(op$par)])
      pit <- as.numeric(stats::plogis(Z %*% g)); cl <- cpb_lp(b, a)
      ## floored as in the FE branch: an unfloored zero-branch log() can hit -Inf
      ## and turn the convergence test into if (NaN)
      ll <- sum(ifelse(is0, log(pmax(pit + (1 - pit) * cl$p0, 1e-12)),
                       log(pmax(1 - pit, 1e-12)) + cl$lp))
      if (is.finite(ll) && (is.null(bestem) || ll > bestem$ll)) bestem <- list(ll = ll, b = b, a = a, g = g, it = it)
      if (is.finite(ll) && is.finite(ll_old) && abs(ll - ll_old) < tol * (abs(ll_old) + tol)) break
      ll_old <- ll
    }
    bestem
  }
  if (method == "em") {
    bestem <- run_em()
    if (is.null(bestem)) stop("zi_cpb (EM): the log-likelihood was non-finite at every iteration; ",
                              "the count component may not cover the data (try method = \"ml\").")
    b <- bestem$b; a <- bestem$a; g <- bestem$g; ll <- bestem$ll; it <- bestem$it
  } else {
    best <- NULL                                                          # direct joint ML (default)
    for (as0 in unique(c(a0, 0.3, 0.7))) {
      op <- tryCatch(stats::optim(c(b0, stats::qlogis(as0), g0), nll, method = "Nelder-Mead",
                                  control = list(maxit = 20L * maxit, reltol = tol)),
                     error = function(e) NULL)
      ## Restart each start from its own solution until the polish stops paying:
      ## a single Nelder-Mead run can stall on a collapsed simplex far from the
      ## optimum (observed: a start stranded 26 log-lik units short recovered
      ## most of the gap in one restart), and the stalled value can then win the
      ## multi-start by accident.
      for (rs in 1:3) {
        if (is.null(op)) break
        op2 <- tryCatch(stats::optim(op$par, nll, method = "Nelder-Mead",
                                     control = list(maxit = 20L * maxit, reltol = tol)),
                        error = function(e) NULL)
        if (is.null(op2) || op2$value > op$value - 1e-4) { if (!is.null(op2) && op2$value < op$value) op <- op2; break }
        op <- op2
      }
      if (!is.null(op) && (is.null(best) || op$value < best$value)) best <- op
    }
    if (is.null(best)) stop("zi_cpb: all optimization starts failed.")
    ## EM cross-check: the EM climbs the same observed-data likelihood by a
    ## different route and cheaply escapes the degenerate pi -> 0 basin that can
    ## strand every Nelder-Mead start (observed: all three ML starts returned a
    ## basin 17 log-lik units below the EM's optimum). Evaluate its solution
    ## under the ML objective, polish, and keep the better optimum.
    bestem <- if (is.null(offset)) tryCatch(run_em(), error = function(e) NULL) else NULL
    if (!is.null(bestem)) {
      st <- c(bestem$b, stats::qlogis(min(max(bestem$a, 1e-6), 1 - 1e-6)), bestem$g)
      v0 <- tryCatch(nll(st), error = function(e) Inf)
      ope <- tryCatch(stats::optim(st, nll, method = "Nelder-Mead",
                                   control = list(maxit = 20L * maxit, reltol = tol)),
                      error = function(e) NULL)
      cand <- if (!is.null(ope) && ope$value <= v0) ope
              else if (is.finite(v0)) list(par = st, value = v0, counts = c(0L, 0L))
              else NULL
      if (!is.null(cand) && cand$value < best$value) best <- cand
    }
    b <- best$par[1:pb]; a <- stats::plogis(best$par[pb + 1]); g <- best$par[(pb + 2):(pb + 1 + pg)]
    ll <- -best$value; it <- as.integer(best$counts[1])
  }
  ## degeneracy guard: the untruncated CPB cannot always satisfy its feasibility
  ## constraint under heavy zero-inflation with a wide count range
  if (mean(cpb_lp(b, a)$lp <= -699) > 0.01)
    warning("zi_cpb: the count component cannot satisfy the CPB feasibility constraint for a ",
            "nontrivial share of observations; the mixture fit may be unreliable and a hurdle ",
            "model is likely more appropriate.")

  names(b) <- colnames(X); names(g) <- colnames(Z)
  lambda_full <- as.numeric(exp(off + X %*% b))
  structure(list(coefficients = b, alpha = a, zero_coef = g, zero.coefficients = g, loglik = ll,
                 iterations = it, offset = off, df = length(b) + 1L + length(g), n = n, y = y,
                 pi_full = as.numeric(stats::plogis(Z %*% g)), lambda_full = lambda_full,
                 ceiling = lambda_full / (1 - a), X = X, Z = Z, fe = NULL,
                 int_beta = b, int_xref = colMeans(X), int_fe_ref = 0,
                 zero_beta = g, zero_xref = colMeans(Z),
                 int_terms = int_terms0, zero_terms = Zt, int_xlev = int_xlev, zero_xlev = zero_xlev,
                 formula = formula, zero_formula = zero,
                 call = match.call()),
            class = "zi_cpb")
}

#' @method print zi_cpb
#' @export
print.zi_cpb <- function(x, ...) {
  cat("Zero-Inflated Continuous Parameter Binomial")
  if (!is.null(x$fe)) cat(" (intensity fixed effects on '", x$fe, "')", sep = "")
  cat("\nCall:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("N: %d   EM iterations: %d   logLik: %.1f\n", x$n, x$iterations, x$loglik))
  cat("\nIntensity (CPB) coefficients:\n"); print(round(x$coefficients, 4))
  cat("Intensity alpha (shape):", round(x$alpha, 4), "\n")
  cat("\nZero-inflation (logit link) -- positive coefficients raise P(structural zero), i.e. lower the\nchance of a positive count (the opposite direction from a hurdle participation model):\n")
  print(round(x$zero_coef, 4))
  cat("Mean structural-zero probability:", round(mean(x$pi_full), 3), "\n")
  invisible(x)
}

#' Predict from a zi_cpb fit
#'
#' @param object A `"zi_cpb"` object.
#' @param newdata Data frame of covariate profiles.
#' @param type `"response"` for the marginal expected count E(Y), `"zero"` for
#'   the structural-zero probability, or `"intensity"` for the CPB mean.
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict zi_cpb
#' @export
predict.zi_cpb <- function(object, newdata = NULL, type = c("response", "zero", "intensity"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    pit <- object$pi_full; lam <- object$lambda_full
  } else {
    miss <- setdiff(unique(c(all.vars(object$int_terms), all.vars(object$zero_terms))), names(newdata))
    if (length(miss))
      stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".")
    Xn <- stats::model.matrix(object$int_terms, newdata, xlev = object$int_xlev)
    Zn <- stats::model.matrix(object$zero_terms, newdata, xlev = object$zero_xlev)
    lam <- if (is.null(object$fe)) as.numeric(exp(Xn %*% object$coefficients))
           else as.numeric(exp(mean(object$fe_hat) +
                               Xn[, names(object$coefficients), drop = FALSE] %*% object$coefficients))
    pit <- as.numeric(stats::plogis(Zn %*% object$zero_coef))
  }
  if (type == "zero") return(pit)
  if (type == "intensity") return(lam)
  (1 - pit) * lam
}
