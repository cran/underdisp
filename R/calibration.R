## ---------------------------------------------------------------------------
## Calibration tooling: predicted pmf, rootogram, PIT, proper scores
## ---------------------------------------------------------------------------

## internal SINGLE SOURCE for the CPB pmf: the normalized untruncated pmf over
## 0..min(floor(nu), kcap), computed in log space. Every CPB pmf consumer
## (.cpb_pmf1 here, .cpb_pmf_row in qoi.R, rcpb in hurdle.R) goes through this
## core, so the numerics cannot drift apart. kcap bounds the allocation against
## a pathological ceiling (alpha -> 1); at 1e5 it never binds in real use.
.cpb_pmf_core <- function(lambda, alpha, kcap = 100000L) {
  nn <- lambda / (1 - alpha); K <- min(floor(nn + 1e-9), kcap); k <- 0:K   # tolerance matches the C++ core
  lw <- lgamma(nn + 1) - lgamma(k + 1) - lgamma(nn - k + 1) +
        k * log(1 - alpha) + (nn - k) * log(alpha)
  pr <- exp(lw - max(lw))
  pr / sum(pr)
}

## zero-truncation of a pmf matrix over its FULL support: the zero column is
## removed and each row divided by 1 - P(Y = 0) (the predictive distribution of
## every zero-truncated fit and of every hurdle intensity)
.zt_rows <- function(P) { p0 <- P[, 1]; P[, 1] <- 0; P / pmax(1 - p0, 1e-12) }

## internal: pmf of the (possibly truncated) CPB for one lambda, over 0..kmax.
## Normalized over the FULL support then cut to kmax (rows may sum to < 1 when
## the ceiling exceeds kmax) -- calibration consumers rely on this semantics.
.cpb_pmf1 <- function(lambda, alpha, kmax, truncated = FALSE) {
  pr <- .cpb_pmf_core(lambda, alpha)
  K  <- length(pr) - 1L
  if (truncated) {
    if (K < 1) { out <- numeric(kmax + 1); if (kmax >= 1) out[2L] <- 1; return(out) }   # ceiling below 1: the point mass at 1
    pr[1] <- 0; pr <- pr / sum(pr)
  }
  out <- numeric(kmax + 1); m <- min(K + 1, kmax + 1); out[seq_len(m)] <- pr[seq_len(m)]
  out
}

## Exact moments of the (possibly zero-truncated) CPB at rate lambda from its
## normalized pmf. The CPB is a binomial with non-integer N = lambda/(1-alpha)
## renormalized over 0..floor(N), so its mean equals lambda and its variance
## alpha*lambda only when N is an integer; the fitted mean, the residuals, the
## first differences, and the dispersion profile all use the exact values.
.cpb_moments <- function(lambda, alpha, truncated = FALSE) {
  out <- vapply(lambda, function(l) {
    pr <- .cpb_pmf_core(l, alpha); k <- seq_along(pr) - 1
    if (truncated) { if (length(pr) < 2L) return(c(1, 0)); pr[1] <- 0; pr <- pr / sum(pr) }   # ceiling below 1: the point mass at 1
    m <- sum(k * pr); c(m, max(sum(k * k * pr) - m * m, 0))
  }, numeric(2))
  list(mean = out[1, ], var = out[2, ])
}
.cpb_mean <- function(lambda, alpha, truncated = FALSE) .cpb_moments(lambda, alpha, truncated)$mean
## conditional mean E(Y | Y >= 1) at rate lambda
.cpb_ztmean <- function(lambda, alpha) .cpb_mean(lambda, alpha, TRUE)

## the CPB rate of a fit: stored as $rate; older objects carry it in fitted.values
.cpb_rate <- function(fit) if (!is.null(fit$rate)) fit$rate else fit$fitted.values

## internal: observed y and an n x (kmax+1) predicted-pmf matrix for a fit
.pmf_and_y <- function(fit, kmax) {
  if (inherits(fit, "hurdle_cpb")) {
    alpha <- fit$intensity$alpha
    P <- t(mapply(function(pi, li) {
      row <- pi * .cpb_pmf1(li, alpha, kmax, truncated = TRUE)
      row[1] <- 1 - pi
      row
    }, fit$p_full, fit$lambda_full))
    y <- fit$y
  } else if (inherits(fit, "zi_cpb")) {
    alpha <- fit$alpha
    P <- t(mapply(function(pi, li) {
      cpb0 <- .cpb_pmf1(li, alpha, kmax, truncated = FALSE)
      row <- (1 - pi) * cpb0
      row[1] <- pi + (1 - pi) * cpb0[1]
      row
    }, fit$pi_full, fit$lambda_full))
    y <- fit$y
  } else if (inherits(fit, "hurdle_gec")) {
    U <- .zt_rows(gec_pmf_cpp(fit$lambda_full, fit$delta, kmax, fit$max.support))
    P <- U * fit$p_full; P[, 1] <- 1 - fit$p_full
    y <- fit$y
  } else if (inherits(fit, "zi_gec")) {
    U <- gec_pmf_cpp(fit$lambda_full, fit$delta, kmax, fit$max.support)
    P <- U * (1 - fit$pi_full); P[, 1] <- fit$pi_full + (1 - fit$pi_full) * U[, 1]
    y <- fit$y
  } else if (inherits(fit, "cpb")) {
    y <- fit$Y                                   # the stored rate carries the offset
    P <- t(vapply(.cpb_rate(fit), .cpb_pmf1, numeric(kmax + 1),
                  alpha = fit$alpha, kmax = kmax, truncated = isTRUE(fit$truncated)))
  } else if (inherits(fit, "cpb_fe")) {
    y <- fit$Y                                   # the rate already carries the fixed effects
    P <- t(vapply(.cpb_rate(fit), .cpb_pmf1, numeric(kmax + 1),
                  alpha = fit$alpha, kmax = kmax, truncated = isTRUE(fit$truncated)))
  } else if (inherits(fit, "gec_fe")) {
    y <- fit$Y                                   # rate already carries the fixed effects
    P <- gec_pmf_cpp(fit$rate, fit$delta, kmax, fit$max.support)
  } else if (inherits(fit, "gec")) {
    y <- fit$Y                                   # the linear predictor carries the offset
    P <- gec_pmf_cpp(as.numeric(exp(fit$linear.predictors)), fit$delta, kmax, fit$max.support)
    if (isTRUE(fit$truncated)) P <- .zt_rows(P)
  } else if (inherits(fit, "hurdle_count")) {
    fam <- .count_fam(fit$family)
    P <- t(mapply(function(pi, li) {
      zt <- fam$pvec(li, fit$theta, kmax); zt[1] <- 0; zt <- zt / max(1 - fam$p0(li, fit$theta), 1e-12)
      row <- pi * zt; row[1] <- 1 - pi; row
    }, fit$p_full, fit$lambda_full))
    y <- fit$y
  } else if (inherits(fit, "zi_count")) {
    fam <- .count_fam(fit$family)
    P <- t(mapply(function(pi, li) {
      f0 <- fam$pvec(li, fit$theta, kmax)
      row <- (1 - pi) * f0; row[1] <- pi + (1 - pi) * f0[1]; row
    }, fit$pi_full, fit$lambda_full))
    y <- fit$y
  } else if (inherits(fit, "count_reg")) {
    fam <- .count_fam(fit$family); y <- fit$Y
    P <- t(vapply(fit$mu, function(mu) {
      v <- fam$pvec(mu, fit$theta, kmax)
      if (isTRUE(fit$truncated)) { v[1] <- 0; v <- v / max(1 - fam$p0(mu, fit$theta), 1e-12) }   # full support
      v
    }, numeric(kmax + 1)))
  } else {
    stop("Unsupported model class: ", paste(class(fit), collapse = "/"))
  }
  list(y = as.integer(y), P = P)
}

.obs_counts <- function(fit)
  if (inherits(fit, c("hurdle_cpb", "zi_cpb", "hurdle_gec", "zi_gec", "hurdle_count", "zi_count"))) fit$y else fit$Y

## internal: linear predictor for a dummy-coded binary equation that tolerates
## factor levels unseen in training (e.g. fixed-effect units present only in a
## held-out fold): unseen levels are scored at the AVERAGE estimated effect of
## that factor's levels (reference level included at zero) -- the binary-equation
## analogue of the mean-fixed-effect fallback the intensity side already uses.
.binary_eta_newdata <- function(terms_obj, xlev, coefs, newdata, contrasts = NULL) {
  coefs[is.na(coefs)] <- 0
  extra <- rep(0, nrow(newdata))
  for (fv in names(xlev)) {
    v <- if (grepl("^factor\\(.*\\)$", fv)) sub("^factor\\((.*)\\)$", "\\1", fv) else fv
    if (!v %in% names(newdata)) next
    vals <- as.character(newdata[[v]])
    seen <- vals %in% xlev[[fv]]
    if (all(seen)) next
    dcoef <- coefs[startsWith(names(coefs), fv)]
    extra[!seen] <- extra[!seen] + (if (length(dcoef)) mean(c(0, unname(dcoef))) else 0)
    vals[!seen] <- xlev[[fv]][1L]                    # reference level for the design
    newdata[[v]] <- vals
  }
  MM <- .ud_newdata_matrix(terms_obj, newdata, xlev, contrasts, names(coefs))
  as.numeric(MM %*% coefs) + extra
}

## internal: observed y and an n x (kmax+1) predicted-pmf matrix for a fit evaluated
## on NEW data (held-out scoring). The dispersion (alpha / theta / delta) is taken
## from the fit; the per-observation rate lambda and, for two-part models, the
## participation/inflation probability, are recomputed from the fitted coefficients
## applied to newdata's design matrices. Fixed-effect units unseen in training fall
## back to the mean fixed effect (intensity) and to the average estimated unit
## effect (participation/inflation dummies). Offsets are assumed absent on newdata.
.pmf_and_y_newdata <- function(fit, newdata, kmax) {
  cpbmat <- function(lam, alpha, trunc)
    t(vapply(lam, .cpb_pmf1, numeric(kmax + 1), alpha = alpha, kmax = kmax, truncated = trunc))
  yof <- function(tf) as.integer(stats::model.response(stats::model.frame(tf, newdata)))

  if (inherits(fit, "cpb")) {
    Tm <- stats::delete.response(fit$terms)
    X  <- .ud_newdata_matrix(Tm, newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    P  <- cpbmat(as.numeric(exp(X %*% fit$coefficients)), fit$alpha, isTRUE(fit$truncated))
    y  <- yof(fit$terms)
  } else if (inherits(fit, "cpb_fe")) {
    X  <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    fe <- fit$fe[as.character(newdata[[fit$fe_var]])]; fe[is.na(fe)] <- mean(fit$fe)
    P  <- cpbmat(as.numeric(exp(fe + X %*% fit$coefficients)), fit$alpha, isTRUE(fit$truncated))
    y  <- yof(fit$formula)
  } else if (inherits(fit, "hurdle_cpb")) {
    pg   <- fit$participation
    peta <- .binary_eta_newdata(stats::delete.response(stats::terms(pg)), pg$xlevels,
                                stats::coef(pg), newdata, contrasts = pg$contrasts)
    p    <- as.numeric(stats::family(pg)$linkinv(peta))
    b    <- fit$int_beta
    Xb   <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$intensity$levels, fit$intensity$contrasts, names(b)) %*% b
    lam  <- if (is.null(fit$fe)) as.numeric(exp(Xb))
            else { fe <- fit$intensity$fe[as.character(newdata[[fit$fe]])]; fe[is.na(fe)] <- fit$int_fe_ref
                   as.numeric(exp(fe + Xb)) }
    a <- fit$intensity$alpha
    P <- t(mapply(function(pi, li) { row <- pi * .cpb_pmf1(li, a, kmax, truncated = TRUE); row[1] <- 1 - pi; row },
                  p, lam))
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else if (inherits(fit, "zi_cpb")) {
    b  <- fit$coefficients
    Xb <- .ud_newdata_matrix(fit$int_terms, newdata, fit$int_xlev, fit$int_contrasts, names(b)) %*% b
    lam <- if (is.null(fit$fe)) as.numeric(exp(Xb))
           else { fe <- fit$fe_hat[as.character(newdata[[fit$fe]])]; fe[is.na(fe)] <- mean(fit$fe_hat)
                  as.numeric(exp(fe + Xb)) }
    pit <- as.numeric(stats::plogis(.binary_eta_newdata(fit$zero_terms, fit$zero_xlev,
                                                        fit$zero_coef, newdata, contrasts = fit$zero_contrasts)))
    a <- fit$alpha
    P <- t(mapply(function(pi, li) { c0 <- .cpb_pmf1(li, a, kmax, truncated = FALSE)
                                     row <- (1 - pi) * c0; row[1] <- pi + (1 - pi) * c0[1]; row }, pit, lam))
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else if (inherits(fit, "count_reg")) {
    fam <- .count_fam(fit$family)
    X   <- .ud_newdata_matrix(fit$terms, newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    lam <- as.numeric(exp(X %*% fit$coefficients))
    P   <- t(vapply(lam, function(mu) { v <- fam$pvec(mu, fit$theta, kmax)
                                        if (isTRUE(fit$truncated)) { v[1] <- 0; v <- v / max(1 - fam$p0(mu, fit$theta), 1e-12) }; v },
                    numeric(kmax + 1)))
    y <- yof(fit$terms)
  } else if (inherits(fit, "hurdle_count")) {
    fam  <- .count_fam(fit$family)
    p    <- .ud_glm_predict(fit$participation, newdata)
    lam  <- exp(as.numeric(predict(fit$intensity, newdata = newdata, type = "link")))   # natural parameter
    P <- t(mapply(function(pi, li) { zt <- fam$pvec(li, fit$theta, kmax); zt[1] <- 0
                                     zt <- zt / max(1 - fam$p0(li, fit$theta), 1e-12)
                                     row <- pi * zt; row[1] <- 1 - pi; row }, p, lam))
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else if (inherits(fit, "zi_count")) {
    fam <- .count_fam(fit$family); linkinv <- stats::make.link(fit$link)$linkinv
    Xc  <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$levels, fit$contrasts)
    Zz  <- .ud_newdata_matrix(if (!is.null(fit$zero_terms)) fit$zero_terms else stats::terms(fit$zero.formula), newdata, fit$levels, fit$zero_contrasts)
    lam <- as.numeric(exp(Xc[, names(fit$coefficients), drop = FALSE] %*% fit$coefficients))
    pit <- as.numeric(linkinv(Zz[, names(fit$zero.coefficients), drop = FALSE] %*% fit$zero.coefficients))
    P <- t(mapply(function(pi, li) { f0 <- fam$pvec(li, fit$theta, kmax)
                                     row <- (1 - pi) * f0; row[1] <- pi + (1 - pi) * f0[1]; row }, pit, lam))
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else if (inherits(fit, "gec_fe")) {                  # check before "gec": gec_fe inherits gec
    Tm  <- stats::delete.response(stats::terms(fit$formula))
    X   <- .ud_newdata_matrix(Tm, newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    fe  <- fit$fe[as.character(newdata[[fit$fe_var]])]; fe[is.na(fe)] <- mean(fit$fe)
    P   <- gec_pmf_cpp(as.numeric(exp(fe + X %*% fit$coefficients)), fit$delta, kmax, fit$max.support)
    y   <- yof(fit$formula)
  } else if (inherits(fit, "gec")) {
    X   <- .ud_newdata_matrix(fit$terms, newdata, fit$levels, fit$contrasts, names(fit$coefficients))
    lam <- as.numeric(exp(X %*% fit$coefficients))
    P   <- gec_pmf_cpp(lam, fit$delta, kmax, fit$max.support)
    if (isTRUE(fit$truncated)) P <- .zt_rows(P)
    y   <- yof(fit$terms)
  } else if (inherits(fit, "hurdle_gec")) {
    p    <- .ud_glm_predict(fit$participation, newdata)
    Xall <- .ud_newdata_matrix(stats::terms(fit$formula), newdata, fit$intensity$levels, fit$intensity$contrasts)
    lam  <- as.numeric(exp(Xall[, names(fit$int_beta), drop = FALSE] %*% fit$int_beta))
    U <- gec_pmf_cpp(lam, fit$delta, kmax, fit$max.support); p0 <- U[, 1]; U[, 1] <- 0; U <- U / pmax(1 - p0, 1e-12)
    P <- U * p; P[, 1] <- 1 - p
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else if (inherits(fit, "zi_gec")) {
    Ti  <- stats::delete.response(fit$int_terms)
    Xn  <- .ud_newdata_matrix(Ti, newdata, fit$int_xlev, fit$int_contrasts, names(fit$coefficients))
    Zn  <- .ud_newdata_matrix(fit$zero_terms, newdata, fit$zero_xlev, fit$zero_contrasts, names(fit$zero_coef))
    lam <- as.numeric(exp(Xn[, names(fit$coefficients), drop = FALSE] %*% fit$coefficients))
    pit <- as.numeric(stats::plogis(Zn[, names(fit$zero_coef), drop = FALSE] %*% fit$zero_coef))
    U <- gec_pmf_cpp(lam, fit$delta, kmax, fit$max.support)
    P <- U * (1 - pit); P[, 1] <- pit + (1 - pit) * U[, 1]
    y <- as.integer(stats::model.response(stats::model.frame(fit$formula, newdata)))
  } else {
    stop("Held-out scoring not supported for model class: ", paste(class(fit), collapse = "/"))
  }
  list(y = as.integer(y), P = P)
}

#' Proper scoring rules for a fitted count model
#'
#' Computes the mean logarithmic score and the ranked probability score (RPS) of
#' a fitted model's predicted distribution against observed counts. Lower is
#' better for both, and both are proper, so they are the natural way to compare
#' the calibration of competing count models.
#'
#' By default the scores are computed **in sample** (against the data the model
#' was fit to, as means weighted by its frequency weights). Supply `newdata` to
#' score a fitted model on **held-out**
#' observations, or use [cv_score()] for a cross-validated score; an in-sample
#' log score equals `-logLik/n` and does not penalize model complexity, so for
#' comparing models of different size the held-out or cross-validated score is
#' the honest criterion.
#'
#' @param fit A fitted `underdisp` count model (`cpb`, `cpb_fe`, `hurdle_cpb`,
#'   `zi_cpb`, `count_reg`, `hurdle_count`, `zi_count`, `gec`, `gec_fe`,
#'   `hurdle_gec`, or `zi_gec`).
#' @param newdata Optional data frame of held-out observations (including the
#'   response). When supplied, the fitted model's predicted distribution is
#'   evaluated on these rows; fixed-effect units unseen in training fall back to
#'   the mean fixed effect. Offsets are assumed absent on `newdata`.
#' @param kmax Highest count to evaluate; defaults to the maximum observed count
#'   in the fitting data, raised to the maximum held-out count when `newdata`
#'   contains larger values. Predicted probabilities are floored at 1e-12 in
#'   the log score, so an observation outside a hard-ceiling model's support
#'   contributes 27.6 to it.
#' @param weights Frequency weights for the rows of `newdata`: a column name of
#'   `newdata` or a numeric vector with one value per row. The held-out scores
#'   are then weighted means, as the in-sample scores of a weighted fit are.
#' @return A named numeric vector `c(logscore, rps)`.
#' @seealso [cv_score()]
#' @examples
#' set.seed(1)
#' y <- rcpb(400, lambda = 3, alpha = 0.5)
#' fit <- cpb(y ~ 1, data = data.frame(y = y), truncated = FALSE, se = "none")
#' score(fit)
#' @export
score <- function(fit, newdata = NULL, kmax = NULL, weights = NULL) {
  if (!is.null(weights) && is.null(newdata))
    stop("'weights' weights the rows of 'newdata'; in sample, the fit's own frequency weights apply.")
  if (is.null(kmax)) kmax <- max(.obs_counts(fit))
  pu <- if (is.null(newdata)) .pmf_and_y(fit, kmax) else .pmf_and_y_newdata(fit, newdata, kmax)
  if (max(pu$y) > kmax) {                       # held-out counts above the in-sample maximum
    kmax <- max(pu$y)
    pu <- if (is.null(newdata)) .pmf_and_y(fit, kmax) else .pmf_and_y_newdata(fit, newdata, kmax)
  }
  y <- pu$y; P <- pu$P
  ## frequency-weighted means: in sample by the fit's own weights, as the
  ## likelihood is; held out by the weights given for the rows of newdata
  w <- if (is.null(newdata)) .ud_fit_w(fit, length(y)) else .ud_weights(weights, newdata, seq_len(nrow(newdata)))
  if (!is.null(w) && length(w) != length(y))
    stop("'weights' must have one value per scored row of 'newdata' (", length(y), ").")
  if (is.null(w)) w <- rep(1, length(y))
  yk <- pmin(y, kmax)
  py <- P[cbind(seq_along(y), yk + 1L)]
  ls <- -sum(w * log(pmax(py, 1e-12))) / sum(w)
  cdf <- t(apply(P, 1, cumsum))
  ind <- outer(yk, 0:kmax, function(a, b) as.numeric(a <= b))
  rps <- sum(w * rowSums((cdf - ind)^2)) / sum(w)
  c(logscore = ls, rps = rps)
}

#' Cross-validated proper scores for a count model
#'
#' K-fold cross-validated log score and RPS: the data are split into `k` folds,
#' the model is refit on each training set via `fitfun`, and its predicted
#' distribution is scored on the held-out fold. Because the score is genuinely
#' out of sample it penalizes over-parameterization, so it is the honest
#' criterion for comparing models that differ in the number of parameters (for
#' example a hurdle with per-unit participation fixed effects against a
#' zero-inflated model). Per-observation held-out scores are returned so two
#' models fit on the *same* folds can be compared with a paired, cluster-robust
#' test.
#'
#' @param fitfun A function of one argument (a training data frame) returning a
#'   fitted `underdisp` model, e.g.
#'   `function(d) hurdle_cpb(y ~ x, d, participation = ~ z, fe = "unit", se = "none")`.
#' @param data The full data frame.
#' @param k Number of folds (default 5).
#' @param kmax Highest count to evaluate; if `NULL`, taken from a fit on the full
#'   data (one extra fit).
#' @param folds Optional integer vector of length `nrow(data)` giving a fixed
#'   fold assignment (so two models can be scored on identical folds). If `NULL`,
#'   folds are drawn at random.
#' @param cores Worker processes for the fold refits (default 1); see [cpb()].
#' @param weights Optional frequency weights for the held-out scores: a column
#'   name of `data` or a numeric vector with one value per row. The mean scores
#'   are then weighted means, so a row of weight w counts as w held-out
#'   observations (the folds assign whole rows). The fold fits see the weights
#'   only through `fitfun`, so pass them there as well.
#' @return An object of class `"cv_score"`: a list with the mean held-out
#'   `logscore` and `rps`, the per-observation vectors `logscore_i`/`rps_i`, and
#'   the `folds` used.
#' @seealso [score()]
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(y = rcpb(500, lambda = 3, alpha = 0.5))
#' cv <- cv_score(function(tr) cpb(y ~ 1, tr, truncated = FALSE, se = "none"), d, k = 5)
#' cv$logscore
#' }
#' @export
cv_score <- function(fitfun, data, k = 5, kmax = NULL, folds = NULL, cores = 1L, weights = NULL) {
  n <- nrow(data)
  wt <- .ud_weights(weights, data, seq_len(n))
  if (is.null(folds)) folds <- sample(rep(seq_len(k), length.out = n))
  if (length(folds) != n) stop("'folds' must have length nrow(data).")
  if (is.null(kmax)) {
    f0 <- tryCatch(fitfun(data), error = function(e) stop("fitfun failed on the full data: ", conditionMessage(e)))
    kmax <- max(.obs_counts(f0))
  }
  ls_i <- rps_i <- rep(NA_real_, n)
  one <- function(kk) {
    idx <- which(folds == kk)
    fit <- tryCatch(fitfun(data[-idx, , drop = FALSE]), error = function(e) NULL)
    if (is.null(fit)) return(list(idx = idx, msg = paste0("fold ", kk, " fit failed; its rows are dropped.")))
    pu <- tryCatch(.pmf_and_y_newdata(fit, data[idx, , drop = FALSE], kmax), error = function(e) NULL)
    if (is.null(pu)) return(list(idx = idx, msg = paste0("fold ", kk, " scoring failed; its rows are dropped.")))
    y <- pu$y; P <- pu$P; yk <- pmin(y, kmax)
    cdf <- t(apply(P, 1, cumsum)); ind <- outer(yk, 0:kmax, function(a, b) as.numeric(a <= b))
    list(idx = idx, ls = -log(pmax(P[cbind(seq_along(y), yk + 1L)], 1e-12)), rps = rowSums((cdf - ind)^2))
  }
  for (r in .ud_lapply(sort(unique(folds)), one, cores)) {
    if (!is.null(r$msg)) { warning(r$msg); next }
    ls_i[r$idx] <- r$ls; rps_i[r$idx] <- r$rps
  }
  ok <- is.finite(ls_i)
  avg <- function(v) if (is.null(wt)) mean(v[ok]) else sum(wt[ok] * v[ok]) / sum(wt[ok])
  structure(list(logscore = avg(ls_i), rps = avg(rps_i),
                 logscore_i = ls_i, rps_i = rps_i, folds = folds, k = k, n_scored = sum(ok)),
            class = "cv_score")
}
#' @export
print.cv_score <- function(x, ...) {
  cat(sprintf("%d-fold cross-validated scores (%d of %d obs scored)\n", x$k, x$n_scored, length(x$folds)))
  cat(sprintf("  held-out log score: %.4f\n  held-out RPS:       %.4f\n", x$logscore, x$rps))
  invisible(x)
}

#' Hanging rootogram for a fitted count model
#'
#' Draws a Tukey hanging rootogram: bars for the observed frequencies hang from
#' the curve of expected frequencies, both on the square-root scale. Bars that
#' hang below the zero line mark counts the model under-predicts; bars that stop
#' short mark counts it over-predicts. A fit with frequency weights counts each
#' row as many times as its weight in both the observed and the expected
#' frequencies.
#'
#' @param fit A fitted `underdisp` count model (any class [score()] takes).
#' @param kmax Highest count to display; defaults to the maximum observed count.
#' @param main,xlab,ylab Plot labels.
#' @param ... Passed to [graphics::plot()].
#' @return Invisibly, a data frame of `count`, `observed`, and `expected`
#'   frequencies.
#' @examples
#' set.seed(1)
#' y <- rcpb(400, lambda = 3, alpha = 0.5)
#' fit <- cpb(y ~ 1, data = data.frame(y = y), truncated = FALSE, se = "none")
#' rootogram(fit)
#' @export
rootogram <- function(fit, kmax = NULL, main = "Rootogram", xlab = "Count",
                      ylab = "sqrt(frequency)", ...) {
  y0 <- .obs_counts(fit)
  if (is.null(kmax)) kmax <- max(y0)
  pu <- .pmf_and_y(fit, kmax); y <- pu$y; P <- pu$P
  ks  <- 0:kmax; wt <- .ud_fit_w(fit, length(y))
  if (is.null(wt)) {
    obs  <- as.numeric(tabulate(factor(pmin(y, kmax), levels = ks), nbins = length(ks)))
    expf <- colSums(P)
  } else {                                      # frequency weights: w copies of a row
    obs  <- vapply(ks, function(k) sum(wt[pmin(y, kmax) == k]), numeric(1))
    expf <- colSums(P * wt)
  }
  ro <- sqrt(obs); re <- sqrt(expf)
  graphics::plot(ks, re, type = "n", ylim = range(0, re, re - ro),
                 xlab = xlab, ylab = ylab, main = main, ...)
  graphics::rect(ks - 0.4, re - ro, ks + 0.4, re, col = "grey85", border = "grey55")
  graphics::lines(ks, re, col = "firebrick", lwd = 2)
  graphics::points(ks, re, col = "firebrick", pch = 19, cex = 0.6)
  graphics::abline(h = 0, col = "grey40")
  invisible(data.frame(count = ks, observed = obs, expected = expf))
}

#' Non-randomized PIT histogram for a fitted count model
#'
#' Draws the non-randomized probability integral transform histogram of
#' Czado, Gneiting, and Held (2009). A well-calibrated model yields a flat
#' histogram at height one (the reference line); a U shape indicates
#' under-dispersion in the predictive distribution and a hump indicates
#' over-dispersion.
#'
#' @param fit A fitted `underdisp` count model (any class [score()] takes); a
#'   fit with frequency weights counts each row as many times as its weight.
#' @param bins Number of histogram bins.
#' @param main,xlab,ylab Plot labels.
#' @param ... Passed to [graphics::barplot()].
#' @return Invisibly, the vector of bin heights (normalized so that a calibrated
#'   model gives heights near one).
#' @examples
#' set.seed(1)
#' y <- rcpb(400, lambda = 3, alpha = 0.5)
#' fit <- cpb(y ~ 1, data = data.frame(y = y), truncated = FALSE, se = "none")
#' pit_hist(fit)
#' @export
pit_hist <- function(fit, bins = 10, main = "PIT histogram", xlab = "PIT",
                     ylab = "Relative frequency", ...) {
  y0 <- .obs_counts(fit); kmax <- max(y0)
  pu <- .pmf_and_y(fit, kmax); y <- pu$y; P <- pu$P
  cdf <- t(apply(P, 1, cumsum)); n <- length(y)
  yk <- pmin(y, kmax)
  Fy   <- cdf[cbind(seq_len(n), yk + 1L)]
  ## F(y - 1), zero at y = 0; an index-matrix row with a zero is dropped from the
  ## result, so the lookup runs over the positive counts only
  Fym1 <- numeric(n); pos <- yk > 0
  Fym1[pos] <- cdf[cbind(which(pos), yk[pos])]
  denom <- pmax(Fy - Fym1, 1e-12)
  u <- seq(0, 1, length.out = bins + 1)
  wt <- .ud_fit_w(fit, n)
  Fbar <- vapply(u, function(uu) { v <- pmin(pmax((uu - Fym1) / denom, 0), 1)
                                   if (is.null(wt)) mean(v) else sum(wt * v) / sum(wt) }, numeric(1))
  h <- diff(Fbar) * bins
  mids <- (u[-1] + u[-(bins + 1)]) / 2
  graphics::barplot(h, names.arg = round(mids, 2), space = 0, main = main,
                    xlab = xlab, ylab = ylab, col = "grey85", border = "grey55", ...)
  graphics::abline(h = 1, col = "firebrick", lwd = 2, lty = 2)
  invisible(h)
}
