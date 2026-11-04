## ---------------------------------------------------------------------------
## Separation in the zero-inflation equation.
##
## A zero contributes log(pi + (1 - pi) f0) to a zero-inflated log-likelihood,
## which rises with pi, and a positive count log(1 - pi) + log f(y), which falls
## with it. A zero-equation coefficient can therefore run to infinity without
## loss when the zeros it pushes towards pi = 1 and the positives it pushes
## towards pi = 0 go there freely, and the zeros it pushes towards pi = 0 are
## covered by the count component's f0. Unlike separation in a logistic
## regression, this needs no linear separation of the zeros from the positives,
## because f0 catches the zeros on the wrong side: a hurdle's participation logit
## on the same design can have a finite maximum while the mixture's zero equation
## has none. The optimizer then stops somewhere along the divergent direction,
## and the reported log-likelihood stands for a supremum.
##
## The check is the profile likelihood of each covariate coefficient of the zero
## equation (not the intercept, whose limits are the no-inflation boundary that
## zi_test() examines, and not the zero_fe dummies, whose limits are boundaries of
## single units, as an all-zero unit is in a fixed-effects logit). The
## coefficient is fixed at its estimate, then 10 and 100 standard deviations of
## its covariate further from zero, and the other zero-equation coefficients are
## re-maximized each time, with the count component held at its estimate. The
## coefficient is flagged when its profile 100 standard deviations out is less
## than qchisq(0.95, 1) / 2 below its value at the estimate: its 95% profile
## interval then reaches past an effect of 100 on the log-odds per standard
## deviation, and the data put no useful bound on it. That covers a coefficient
## with no finite estimate (the profile flat to infinity) and one whose finite
## maximum lies on a near-flat ridge where the fitted pi are numerically 0 or 1
## (on real panels, a profile can lose well under one log-likelihood unit over
## 100 standard deviations), which a fit cannot tell apart. Holding the count
## component fixed can only lower the profile, so the full model's interval is at
## least as wide. The mixture likelihood is not concave, and a profile can rise
## above its value at the estimate before it falls; such a point is a better
## point of the likelihood. So is the re-maximization at the estimate itself
## when the fit is not a stationary point of the zero equation (an optimizer that
## stopped early). The best point found, against the fit, is returned as
## `short_of_max`, and the fits restart their final polish from it (.zi_restart).
##
## Every profile can reach the no-inflation boundary, where each pi is zero and
## the zero equation drops out: the intercept goes to minus infinity whatever
## the fixed coefficient. A coefficient whose profile at the estimate lies within
## qchisq(0.95, 1) / 2 of that limit is unbounded because the data show little
## inflation, not by separation (a zero-inflated fit to data with no structural
## zeros), so it is not flagged; the fit is then reported as near the boundary
## (`no_inflation`), which zi_test() tests.
##
## Alongside, the observations whose pi is numerically 0 or 1 are counted (glm's
## rule, 10 * .Machine$double.eps), leaving out the groups (zero_fe units, else
## the whole sample) whose every pi lies within 1e-8 of the same edge.
## ---------------------------------------------------------------------------

## Z, g: the zero equation's design and coefficients (original units); y: the
## counts; f0: the count component's P(Y = 0) at the fit; w: frequency weights;
## group: the zero_fe values (NULL without zero-equation fixed effects); fe_cols:
## the names of Z's zero_fe dummy columns.
.zi_separation <- function(Z, g, y, f0, w = NULL, group = NULL, fe_cols = character(0), link = "logit",
                           steps = c(10, 100), bound = 100) {
  n <- length(y); is0 <- y == 0; w <- if (is.null(w)) rep(1, n) else as.numeric(w)
  g <- as.numeric(g); f0 <- pmin(pmax(as.numeric(f0), 0), 1); lf0 <- log(f0); l1f0 <- log1p(-f0)
  lk <- stats::make.link(link); logit <- identical(link, "logit")
  ## per-observation log-likelihood of the zero/positive split and its derivative in the index
  lld <- function(e) {
    if (logit) {
      lp <- stats::plogis(e, log.p = TRUE); lq <- stats::plogis(-e, log.p = TRUE)
      a <- pmax(lp, lq + lf0); lden <- a + log(exp(lp - a) + exp(lq + lf0 - a))
      list(l = ifelse(is0, lden, lq), d = ifelse(is0, exp(l1f0 + lp + lq - lden), -exp(lp)))
    } else {
      p <- lk$linkinv(e); me <- lk$mu.eta(e); den <- pmax(p + (1 - p) * f0, 1e-300); q <- pmax(1 - p, 1e-300)
      list(l = ifelse(is0, log(den), log(q)), d = ifelse(is0, (1 - f0) * me / den, -me / q))
    }
  }
  eta <- as.numeric(Z %*% g); pit <- lk$linkinv(eta); l_fit <- sum(w * lld(eta)$l)
  ## the no-inflation limit (every pi at zero: each zero contributes log f0, each positive nothing)
  cut <- stats::qchisq(0.95, 1) / 2; l_none <- sum(w[is0] * lf0[is0])
  ## numerically 0 or 1, outside groups that sit wholly at one edge
  grp <- if (is.null(group)) rep(1L, n) else as.integer(factor(group))
  edge <- as.logical(stats::ave(pit, grp, FUN = function(p) all(p < 1e-8) || all(p > 1 - 1e-8)))
  ext <- (pit < 10 * .Machine$double.eps | pit > 1 - 10 * .Machine$double.eps) & !edge
  ## the profile of each covariate coefficient, the rest of the zero equation re-maximized
  cols <- setdiff(seq_len(ncol(Z)), c(which(colnames(Z) == "(Intercept)"), which(colnames(Z) %in% fe_cols)))
  ## two starts for the re-maximization: the other coefficients at the fit, and all of them scaled with the fixed
  ## one (the ray through the fit), since on a saturated ridge the gradient vanishes and a search from the fit
  ## alone cannot find the compensating move (a threshold in the covariate keeps its place only along the ray)
  prof <- function(j, gj) {
    free <- setdiff(seq_along(g), j); zj <- Z[, j] * gj; Zf <- Z[, free, drop = FALSE]
    f <- function(v) -sum(w * lld(zj + as.numeric(Zf %*% v))$l)
    gr <- function(v) -as.numeric(crossprod(Zf, w * lld(zj + as.numeric(Zf %*% v))$d))
    best <- -Inf; bv <- NULL
    for (st in list(g[free], g[free] * (gj / g[j]))) {
      v0 <- -f(st); if (is.finite(v0) && v0 > best) { best <- v0; bv <- st }
      o <- tryCatch(stats::optim(st, f, gr, method = "BFGS", control = list(maxit = 500L, reltol = 1e-12)),
                    error = function(e) NULL)
      if (!is.null(o) && is.finite(o$value) && -o$value > best) { best <- -o$value; bv <- o$par }
    }
    gg <- g; gg[j] <- gj; if (!is.null(bv)) gg[free] <- bv
    list(value = best, g = gg)
  }
  tab <- NULL; top <- list(value = l_fit, g = NULL)
  for (j in cols) {
    sdj <- stats::sd(Z[, j])
    if (!is.finite(g[j]) || g[j] == 0 || !is.finite(sdj) || sdj <= 0) next
    s <- sign(g[j]) / sdj; st <- sort(unique(c(steps, bound)))
    pr <- lapply(c(0, st), function(k) prof(j, g[j] + k * s))
    for (q in pr) if (is.finite(q$value) && q$value > top$value) top <- q
    pv <- vapply(pr, `[[`, numeric(1), "value"); p0 <- pv[1L]; ps <- pv[-1L] - p0
    row <- data.frame(term = colnames(Z)[j], estimate = g[j], profile = p0)
    row[paste0("step", st)] <- as.list(ps)
    row$unbounded <- ps[st == bound] > -cut && p0 - l_none >= cut
    tab <- rbind(tab, row)
  }
  ## a profile point above the fit, at the estimate or beyond it, is a better point of the likelihood (the count
  ## component unchanged): the optimizer stopped at least that far short of the maximum, and `best` holds the
  ## zero-equation coefficients there, from which the fit restarts (.zi_restart)
  sc <- top$value - l_fit; short <- is.finite(sc) && sc > 1e-6 * (1 + abs(l_fit))
  list(diverging = if (is.null(tab)) character(0) else tab$term[tab$unbounded], profile = tab,
       short_of_max = if (short) sc else 0, best = if (short) stats::setNames(top$g, colnames(Z)),
       no_inflation = isTRUE(top$value - l_none < cut), inflation_gain = top$value - l_none,
       n_extreme = sum(ext), groups_extreme = length(unique(grp[ext])))
}

## The fit restarted from the best point the check found, while that pays: `check(fit)` runs .zi_separation on a
## fit, `refit(fit, g)` re-polishes the whole model from the fit with its zero-equation coefficients set to g, and
## `loglik(fit)` reads the log-likelihood. A restart is kept only when it gains; at most `rounds` of them, while the
## check finds at least `min_gain`. Returns the final fit and its check, which records the restarts and their gain.
.zi_restart <- function(fit, check, refit, loglik, rounds = 3L, min_gain = 0.01) {
  sep <- check(fit); k <- 0L; gain <- 0
  while (k < rounds && sep$short_of_max >= min_gain && !is.null(sep$best)) {
    nf <- tryCatch(suppressWarnings(refit(fit, sep$best)), error = function(e) NULL)
    if (is.null(nf) || !is.finite(loglik(nf)) || loglik(nf) <= loglik(fit) + 1e-8 * (1 + abs(loglik(fit)))) break
    gain <- gain + loglik(nf) - loglik(fit); fit <- nf; k <- k + 1L
    sep <- check(fit)
  }
  sep$restarts <- k; sep$restart_gain <- gain
  list(fit = fit, separation = sep)
}

## The zero design's zero_fe dummy columns
.zi_fe_cols <- function(Z, zero_fe)
  if (is.null(zero_fe)) character(0) else colnames(Z)[.is_fe_col(colnames(Z), zero_fe)]

## The fit's warning (the fit the user asked for, not a start or a bootstrap refit)
.zi_warn_separation <- function(sep, fun) {
  if (is.null(sep) || !length(sep$diverging)) return(invisible(NULL))
  pl <- length(sep$diverging) > 1L
  warning(fun, ": the zero-inflation equation is separated or nearly so. The profile likelihood of the coefficient",
          if (pl) "s" else "", " on ", paste(sep$diverging, collapse = ", "), " falls by less than ",
          "qchisq(0.95, 1) / 2 over 100 standard deviations of ", if (pl) "their covariates" else "its covariate",
          " further from zero, so the data put no useful bound on ", if (pl) "them" else "it", ". The reported ",
          "coefficients and log-likelihood are where the optimizer stopped on a flat or divergent ridge (",
          sep$n_extreme, " observations have a structural-zero probability numerically 0 or 1). Compare the ",
          "log-likelihood with care; see the Separation section of ?zi_cpb.", call. = FALSE)
}

## The print methods' note
.zi_separation_note <- function(x) {
  sep <- x$separation
  if (is.null(sep)) return(invisible(NULL))
  if (length(sep$diverging))
    cat("Note: the zero-inflation equation is separated or nearly so (no useful bound on: ",
        paste(sep$diverging, collapse = ", "), "); see ?zi_cpb, Separation.\n", sep = "")
  if (isTRUE(sep$no_inflation))
    cat(sprintf(paste0("Note: the fit is near the no-inflation boundary (structural-zero probabilities of zero lower the ",
                       "log-likelihood by %.2f, less than qchisq(0.95, 1) / 2), so the zero-inflation coefficients are not ",
                       "identified; zi_test() tests for inflation.\n"), max(sep$inflation_gain, 0)))
  if (isTRUE(sep$restarts > 0L))
    cat(sprintf("Note: the fit was restarted from a better point the separation check found (%d restart%s, %+.2f log-likelihood).\n",
                sep$restarts, if (sep$restarts > 1L) "s" else "", sep$restart_gain))
  if (isTRUE(sep$short_of_max >= 0.01))
    cat(sprintf("Note: the separation check found a point %.2f log-likelihood units higher; the optimizer stopped short.\n",
                sep$short_of_max))
}
