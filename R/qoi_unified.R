## ---------------------------------------------------------------------------
## The unified first-difference contract ("ud_fd"): every first_difference()
## method in the package returns the same object -- a data.frame of class
## c("ud_fd", "data.frame") with columns
##   component | from | to | diff | lower | upper | method
## Single-equation fits return one row (the requested quantity); two-part fits
## return one row per margin LEVEL: the binary-stage probability, the
## count-stage mean, and the marginal mean, each evaluated at `from` and `to`
## (so the rows are quantities, and `stage=` controls WHERE the covariate
## moves: in both equations, or in one holding the other at its reference).
## `method` records the uncertainty source per row: "bootstrap (stored)",
## "bootstrap (refit, B=...)", "delta", or "none" (bounds NA); the rule is
## stated at .qoi_ci() below.
## Unknown arguments error loudly (no silent ... swallowing).
## ---------------------------------------------------------------------------

## strict-dots guard: called by every method as .fd_dots(...)
.fd_dots <- function(...) {
  if (...length() > 0L) {
    nm <- names(list(...))
    nm <- if (is.null(nm)) "<positional>" else ifelse(nzchar(nm), nm, "<positional>")
    stop("unused argument(s) for first_difference(): ", paste(nm, collapse = ", "),
         ". See ?first_difference for this model's arguments.", call. = FALSE)
  }
  invisible(NULL)
}

## constructor: assemble + class the unified return
.ud_fd <- function(component, from, to, lower = NA_real_, upper = NA_real_,
                   method = "none") {
  out <- data.frame(component = component, from = from, to = to,
                    diff = to - from,
                    lower = rep_len(lower, length(component)),
                    upper = rep_len(upper, length(component)),
                    method = rep_len(method, length(component)),
                    row.names = NULL, stringsAsFactors = FALSE)
  class(out) <- c("ud_fd", "data.frame")
  out
}

#' @export
print.ud_fd <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  if (all(is.na(y$lower))) y$lower <- y$upper <- NULL
  if (all(y$method == "none")) y$method <- NULL
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

## ---------------------------------------------------------------------------
## The unified rate/odds-ratio contract ("ud_irr"): every irr() method returns
## a data.frame of class c("ud_irr","data.frame") with columns
##   term | equation ("count", "participation" or "inflation") | ratio ("IRR" or "OR") |
##   estimate | lower | upper | method
## Two-part fits stack both equations; single-equation fits return the count
## rows only. `method` records the interval source per row.
## ---------------------------------------------------------------------------

.ud_irr <- function(term, equation, ratio, estimate,
                    lower = NA_real_, upper = NA_real_, method = "none") {
  out <- data.frame(term = term, equation = rep_len(equation, length(term)),
                    ratio = rep_len(ratio, length(term)),
                    estimate = estimate,
                    lower = rep_len(lower, length(term)),
                    upper = rep_len(upper, length(term)),
                    method = rep_len(method, length(term)),
                    row.names = NULL, stringsAsFactors = FALSE)
  class(out) <- c("ud_irr", "data.frame")
  out
}

## The scale of a fit's count-equation coefficients and the ratio irr() names exp(beta): the log mean and "IRR" for
## the mean-parameterized families; the log rate for the CPB and GEC (whose mean is their rate renormalized over a
## finite support; "IRR", with ?irr.count's caveat); otherwise the parameter exp(beta) multiplies, whose ratio is
## not a ratio of means (ame() and first_difference() give those)
.coef_scale <- function(x) {
  if (inherits(x, c("cpb", "cpb_fe", "hurdle_cpb", "zi_cpb", "gec", "hurdle_gec", "zi_gec")))
    return(list(scale = "log rate", ratio = "IRR"))
  fam <- if (is.character(x$family)) x$family else ""
  switch(fam,
         compois = list(scale = "log lambda", ratio = "lambda ratio"),
         ewp3 = , ewp2 = list(scale = "log centrality", ratio = "centrality ratio"),
         gammacount = list(scale = "log long-run rate", ratio = "long-run rate ratio"),
         doublepois = list(scale = "log location", ratio = "location ratio"),
         genpois = if (isTRUE(x$theta < 0)) list(scale = "log mu", ratio = "mu ratio")
                   else list(scale = .mean_scale(x), ratio = "IRR"),
         list(scale = .mean_scale(x), ratio = "IRR"))
}
## a zero-truncated fit's or a hurdle intensity's coefficients are on the log of the untruncated mean
.mean_scale <- function(x) if (isTRUE(x$truncated) || inherits(x, "hurdle_count")) "log untruncated mean" else "log mean"
## a short name of a fit's model, for a table column
.model_short <- function(x) {
  fams <- c(poisson = "Poisson", negbin = "NB", compois = "COM-Poisson", mpcmp = "MP COM-Poisson", genpois = "GP",
            gammacount = "Gamma count", doublepois = "Double Poisson", hyperpois = "Hyper-Poisson", ewp3 = "EWP3",
            ewp2 = "EWP2")
  k <- intersect(class(x), c("cpb_fe", "gec_fe", "cpb", "gec", "count_reg", "hurdle_cpb", "hurdle_gec", "hurdle_count",
                             "zi_cpb", "zi_gec", "zi_count"))[1L]
  f <- if (is.character(x$family) && x$family %in% names(fams)) fams[[x$family]] else ""
  switch(k, cpb = if (isTRUE(x$truncated)) "ZT CPB" else "CPB", cpb_fe = "CPB FE", gec = "GEC", gec_fe = "GEC FE",
         hurdle_cpb = "Hurdle CPB", hurdle_gec = "Hurdle GEC", zi_cpb = "ZI CPB", zi_gec = "ZI GEC",
         count_reg = paste0(if (isTRUE(x$truncated)) "ZT ", f, if (!is.null(x$fe)) " FE"),
         hurdle_count = paste("Hurdle", f), zi_count = paste("ZI", f), "")
}
## irr()'s count rows of a count family, under the ratio .coef_scale() names, with a message where it is not a ratio
## of means
.count_irr_rows <- function(object, b, se, level) {
  r <- .coef_scale(object)$ratio
  if (r != "IRR")
    message("irr(): exp(coefficient) of the ", .count_fam(object$family)$label, " count equation is a ", r,
            ", not a ratio of means; ame() and first_difference() give the mean contrasts.")
  .ud_irr_wald(b, se, level, "count", r)
}

## Wald rows on the ratio scale from coefficients + SEs (NULL se -> no interval; a non-finite se -> none for its row)
.ud_irr_wald <- function(b, se, level, equation, ratio, method = "wald") {
  if (is.null(se)) return(.ud_irr(names(b), equation, ratio, exp(unname(b))))
  se <- unname(se); ok <- is.finite(se); z <- stats::qnorm(1 - (1 - level) / 2)
  .ud_irr(names(b), equation, ratio, exp(unname(b)),
          lower = ifelse(ok, exp(unname(b) - z * se), NA_real_), upper = ifelse(ok, exp(unname(b) + z * se), NA_real_),
          method = ifelse(ok, method, "none"))
}

## ratio rows from the stored bootstrap draws: the percentile interval of exp(coefficient), which is the
## exponential of the coefficient's percentile interval (confint())
.ud_irr_draws <- function(object, b, level, equation = "count") {
  D <- .ud_draws(object, names(b))
  if (is.null(D)) return(.ud_irr_wald(b, NULL, level, equation, "IRR"))
  a <- (1 - level) / 2
  .ud_irr(names(b), equation, "IRR", exp(unname(b)), lower = exp(apply(D, 2L, stats::quantile, a, names = FALSE)),
          upper = exp(apply(D, 2L, stats::quantile, 1 - a, names = FALSE)), method = "bootstrap (stored)")
}

## The rows of terms the data do not bound (separation; see ?zi_cpb) are reported as NA: a ratio such as 6e7 with
## the interval [0, Inf] is where the optimizer stopped on a divergent ridge, not an estimate
.irr_mask <- function(tab, terms, equation) {
  if (is.null(tab) || !length(terms)) return(tab)
  hit <- tab$term %in% terms & tab$equation == equation
  if (any(hit)) {
    tab$estimate[hit] <- NA_real_; tab$lower[hit] <- NA_real_; tab$upper[hit] <- NA_real_; tab$method[hit] <- "separated"
    message("irr(): the data put no useful bound on ", paste(tab$term[hit], collapse = ", "), " in the ", equation,
            " equation (separation), so ", if (sum(hit) > 1L) "their ratios are" else "its ratio is", " NA.")
  }
  tab
}
## a ratio table without the dummies of the fixed-effects variables `fe`
.irr_drop_fe <- function(tab, fe) if (is.null(tab) || is.null(fe)) tab else tab[!.is_fe_col(tab$term, fe), , drop = FALSE]

#' @export
print.ud_irr <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, digits)
  if (all(is.na(y$lower))) y$lower <- y$upper <- NULL
  if (all(y$method == "none")) y$method <- NULL
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

## ---------------------------------------------------------------------------
## The uncertainty of a quantity of interest: one rule for every class.
## A quantity q(theta) of the fitted parameters theta gets its interval
##   (1) from the fit's stored bootstrap draws (se = "bootstrap" at fit time),
##       jointly over the coefficients, the dispersion parameter and, for a
##       fixed-effects fit, the reference unit effect: q is computed at every
##       draw and the percentile interval taken. The draws of a fixed-effects fit
##       are first recentered on the reported estimate, since its replicates
##       share the incidental-parameters bias (the reason its coefficient
##       intervals are normal ones centred on the estimate) and a
##       jackknife-corrected fit reports the corrected estimate. A recentred
##       draw whose dispersion parameter leaves its range (`bounds`, an
##       interval (lower, upper] per parameter: the CPB's alpha above 1, say)
##       describes no distribution and is dropped;
##   (2) otherwise by the delta method over the fit's full covariance, the
##       dispersion parameter included: the count families' analytic, robust or
##       clustered covariance, and for a hurdle the covariances of its two
##       margins, which are independent because its likelihood factorizes (the
##       intensity's from its stored draws);
##   (3) otherwise not at all (method "none"); a row whose parameters lack a
##       covariance has no interval.
## theta is named as the draws' columns or the covariance's dimnames. `hold`
## names parameters the draws do not carry (a zero_fe fit's unit dummies, which
## are not the same units across replicates), held at their estimates. The
## standard error is the draws' standard deviation, or the delta method's.
## `each`: the quantities (predictions, one per row) take their draws one by one,
## each over the draws at which it is finite; otherwise a draw at which any
## quantity is not finite is dropped for all of them.
## ---------------------------------------------------------------------------
.qoi_ci <- function(qfun, theta, draws = NULL, vc = NULL, level = 0.95, recenter = FALSE, hold = character(0),
                    each = FALSE, bounds = NULL) {
  q0 <- qfun(theta); k <- length(q0); a <- (1 - level) / 2
  if (!is.null(draws)) {
    for (h in setdiff(intersect(hold, names(theta)), colnames(draws)))
      draws <- cbind(draws, matrix(theta[[h]], nrow(draws), 1L, dimnames = list(NULL, h)))
    if (all(names(theta) %in% colnames(draws))) {
      D <- draws[, names(theta), drop = FALSE]
      if (recenter) {
        D <- sweep(D, 2L, colMeans(D) - theta)
        for (nm in intersect(names(bounds), colnames(D)))
          D[D[, nm] <= bounds[[nm]][1L] | D[, nm] > bounds[[nm]][2L], ] <- NA
      }
      Q <- matrix(unlist(lapply(seq_len(nrow(D)), function(i)
                    if (anyNA(D[i, ])) rep(NA_real_, k) else qfun(stats::setNames(D[i, ], names(theta))))),
                  ncol = k, byrow = TRUE)
      if (each) {
        Q[!is.finite(Q)] <- NA
        nf <- colSums(!is.na(Q)); use <- nf >= 2L
        if (any(use)) {
          qq <- function(p) ifelse(use, apply(Q, 2L, stats::quantile, p, na.rm = TRUE, names = FALSE), NA_real_)
          return(list(lower = qq(a), upper = qq(1 - a),
                      se = ifelse(use, apply(Q, 2L, stats::sd, na.rm = TRUE), NA_real_),
                      method = ifelse(use, "bootstrap (stored)", "none")))
        }
      } else {
        Q <- Q[stats::complete.cases(Q), , drop = FALSE]
        if (nrow(Q) >= 2L)
          return(list(lower = apply(Q, 2L, stats::quantile, a, names = FALSE),
                      upper = apply(Q, 2L, stats::quantile, 1 - a, names = FALSE),
                      se = apply(Q, 2L, stats::sd), method = rep("bootstrap (stored)", k)))
      }
    }
  }
  lo <- hi <- sev <- rep(NA_real_, k); m <- rep("none", k)
  if (!is.null(vc) && all(names(theta) %in% rownames(vc))) {
    V <- vc[names(theta), names(theta), drop = FALSE]
    ## the derivative is taken along each parameter in steps of its standard error (a thousandth, refined by
    ## Richardson extrapolation), so it does not depend on the covariates' units: a step relative to the
    ## parameter's value moves the linear index by a hundred or more for a coefficient on a covariate in persons
    dv <- diag(V); bad <- !is.finite(dv) | dv <= 0               # a generalized inverse can leave a negative variance
    s <- numeric(length(dv)); s[!bad] <- sqrt(dv[!bad]); s[bad] <- pmax(abs(theta[bad]), 1) * 1e-2
    G <- tryCatch(sweep(matrix(numDeriv::jacobian(function(u) qfun(theta + u * s), rep(0, length(theta)),
                                                  method.args = list(eps = 1e-3, d = 1e-4, zero.tol = Inf, r = 4,
                                                                     v = 2)), nrow = k), 2L, s, "/"),
                  error = function(e) NULL)
    if (!is.null(G)) {
      z <- stats::qnorm(1 - a)
      for (r in seq_len(k)) {
        u <- which(G[r, ] != 0)
        if (anyNA(G[r, ]) || anyNA(V[u, u])) next
        se <- if (length(u)) sqrt(max(drop(G[r, u] %*% V[u, u, drop = FALSE] %*% G[r, u]), 0)) else 0
        lo[r] <- q0[r] - z * se; hi[r] <- q0[r] + z * se; sev[r] <- se; m[r] <- "delta"
      }
    }
  }
  list(lower = lo, upper = hi, se = sev, method = m)
}

## a covariance over the parameters `nm` from independent blocks (named matrices): the blocks on the diagonal,
## zero between them, NA where no block covers a parameter
.qoi_vc <- function(nm, ...) {
  V <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  blocks <- Filter(Negate(is.null), list(...))
  sets <- lapply(blocks, function(B) intersect(rownames(B), nm))
  for (k in seq_along(blocks)) V[sets[[k]], sets[[k]]] <- blocks[[k]][sets[[k]], sets[[k]]]
  for (k in seq_along(sets)) for (l in seq_along(sets)) if (k != l) V[sets[[k]], sets[[l]]] <- 0
  V
}

## a first difference from a class's specification (.fd_spec): `lev(theta)` gives the quantities at `from` and
## `to` (list(from, to), one value per row of the result)
.fd_build <- function(spec, level) {
  pt <- spec$lev(spec$theta)
  ci <- .qoi_ci(function(th) { L <- spec$lev(th); L$to - L$from }, spec$theta, spec$draws, spec$vc, level,
                isTRUE(spec$recenter), if (is.null(spec$hold)) character(0) else spec$hold, bounds = spec$bounds)
  .ud_fd(spec$component, pt$from, pt$to, ci$lower, ci$upper, ci$method)
}

## `variable` must be one model-matrix column that no other column contains; `at` fixes other columns
.fd_check_var <- function(variable, cols, at, msg = NULL) {
  if (!is.character(variable) || length(variable) != 1L || !variable %in% cols)
    stop(if (is.null(msg)) paste0("'variable' must name a model-matrix column: ", paste(cols, collapse = ", ")) else msg,
         call. = FALSE)
  .fd_check_terms(variable, cols)
  .fd_check_at(at, variable, cols)
}
.fd_check_at <- function(at, variable, cols) {
  if (is.null(at)) return(invisible(NULL))
  if (!is.numeric(at) || is.null(names(at)) || any(!nzchar(names(at))) || anyNA(at))
    stop("'at' must be a named numeric vector of model-matrix columns and their values, e.g. at = c(z = 1).",
         call. = FALSE)
  bad <- setdiff(names(at), cols)
  if (length(bad))
    stop("'at' names column(s) not in the model: ", paste(bad, collapse = ", "), ". The model-matrix columns are: ",
         paste(cols, collapse = ", "), ".", call. = FALSE)
  if (variable %in% names(at))
    stop("'at' fixes the other columns; 'variable' moves from 'from' to 'to'.", call. = FALSE)
  invisible(NULL)
}
.fd_at <- function(xref, at) {
  if (!is.null(at)) { k <- intersect(names(at), names(xref)); xref[k] <- at[k] }
  xref
}
## the linear index of a profile `x` under the parameters `th`, the coefficients named with `prefix`
.fd_index <- function(x, th, prefix = "") sum(x * th[paste0(prefix, names(x))])

## A two-part fit's reference profiles at `from` and `to`: the covariate moves in the count equation (`x`), the
## binary one (`z`), or both, by `stage`.
.fd_twopart_profiles <- function(xref, zref, variable, from, to, at, stage, count_name, binary_name) {
  inx <- variable %in% names(xref); inz <- variable %in% names(zref)
  if (!inx && !inz)
    stop("'variable' is in neither the ", count_name, " nor the ", binary_name, " equation.", call. = FALSE)
  cols <- unique(c(names(xref), names(zref)))
  .fd_check_terms(variable, cols); .fd_check_at(at, variable, cols)
  xref <- .fd_at(xref, at); zref <- .fd_at(zref, at)
  mx <- inx && stage %in% c("both", "intensity", "count"); mz <- inz && stage %in% c("both", "participation", "zero")
  list(x = list(from = if (mx) replace(xref, variable, from) else xref, to = if (mx) replace(xref, variable, to) else xref),
       z = list(from = if (mz) replace(zref, variable, from) else zref, to = if (mz) replace(zref, variable, to) else zref))
}
## the levels of a two-part fit's rows: the binary-stage probability, the count-stage mean, and the marginal mean
.fd_twopart_lev <- function(P, pfun, mfun, marginal) function(th) {
  pf <- pfun(P$z$from, th); pt <- pfun(P$z$to, th); mf <- mfun(P$x$from, th); mt <- mfun(P$x$to, th)
  list(from = c(pf, mf, marginal(pf, mf)), to = c(pt, mt, marginal(pt, mt)))
}

## The specification of a class's first difference: the levels `lev(theta)` and the components, with the
## parameters `theta` and the draws or covariance its interval comes from (.qoi_par(), qoi_predict.R; see .qoi_ci).
.fd_spec <- function(object, ...) UseMethod(".fd_spec")

## A single-equation fit's parameters: the coefficients by their names, then the dispersion as ".disp" and a
## fixed-effects fit's reference unit as ".fe_ref" (names no model-matrix column takes, so a covariate called
## alpha cannot collide with them). Its draws are in that layout (bootstrap.R), and are named the same.
.fd_theta1 <- function(object, disp, fe_ref = NULL) {
  theta <- c(object$coefficients, .disp = disp, if (!is.null(fe_ref)) c(.fe_ref = fe_ref))
  D <- .ud_draws(object)
  if (!is.null(D) && ncol(D) == length(theta)) colnames(D) <- names(theta) else D <- NULL
  list(theta = theta, draws = D)
}

.fd_spec.cpb <- function(object, variable, from, to, at = NULL, quantity = "mean", y = NULL, ...) {
  b <- object$coefficients; .fd_check_var(variable, names(b), at)
  if (quantity == "prob" && is.null(y)) stop("quantity = \"prob\" requires 'y'.", call. = FALSE)
  x0 <- .fd_at(colMeans(object$X), at); off <- .mean_offset(object); sup <- .cpb_support(object)
  tr <- isTRUE(object$truncated)
  lev <- function(th) {
    l <- exp(off + c(.fd_index(replace(x0, variable, from), th), .fd_index(replace(x0, variable, to), th)))
    a <- th[[".disp"]]
    q <- switch(quantity, mean = .cpb_mean(l, a, tr, sup), ceiling = l / (1 - a),
                prob = .cpb_prob_at(l, a, c(y, y), tr, object$max.support, support = sup))
    list(from = q[1L], to = q[2L])
  }
  c(list(component = quantity, lev = lev), .qoi_unc(.qoi_par(object)))
}

.fd_spec.cpb_fe <- function(object, variable, from, to, at = NULL, ...) {
  b <- object$coefficients; .fd_check_var(variable, names(b), at, "'variable' is not a covariate in the model.")
  x0 <- .fd_at(object$xref, at); off <- .mean_offset(object); sup <- .cpb_support(object)
  tr <- isTRUE(object$truncated)
  lev <- function(th) {
    m <- .cpb_mean(exp(th[[".fe_ref"]] + off + c(.fd_index(replace(x0, variable, from), th),
                                                .fd_index(replace(x0, variable, to), th))), th[[".disp"]], tr, sup)
    list(from = m[1L], to = m[2L])
  }
  c(list(component = "mean", lev = lev), .qoi_unc(.qoi_par(object)))
}

## gec() and gec_fe() (whose class is c("gec_fe", "gec"))
.fd_spec.gec <- function(object, variable, from, to, at = NULL, ...) {
  fe <- inherits(object, "gec_fe")
  b <- object$coefficients; .fd_check_var(variable, names(b), at, "'variable' is not a covariate in the model.")
  x0 <- .fd_at(if (fe) object$xref else colMeans(object$X), at); off <- .mean_offset(object)
  ms <- object$max.support; tr <- isTRUE(object$truncated)
  lev <- function(th) {
    d <- th[[".disp"]]; base <- off + if (fe) th[[".fe_ref"]] else 0
    mu <- exp(base + c(.fd_index(replace(x0, variable, from), th), .fd_index(replace(x0, variable, to), th)))
    m <- gec_mean_cpp(mu, d, ms)
    if (tr) m <- m / pmax(1 - gec_pmf_cpp(mu, d, 0L, ms)[, 1], 1e-12)          # E(Y | Y > 0)
    list(from = m[1L], to = m[2L])
  }
  c(list(component = "mean", lev = lev), .qoi_unc(.qoi_par(object)))
}

.fd_spec.count_reg <- function(object, variable, from, to, at = NULL, ...) {
  fam <- .count_fam(object$family); b <- object$coefficients
  .fd_check_var(variable, names(b), at, "'variable' is not a covariate in the model.")
  x0 <- .fd_at(.fe_dummy_ref(colMeans(object$X), object$X, object$Y, object$fe,
                             if (isTRUE(object$truncated)) .ident_trunc else .ident_count), at)
  off <- .mean_offset(object); sh <- if (fam$nshape) fam$shape_name else character(0)
  lev <- function(th) {
    tn <- if (fam$nshape) unname(fam$shape_inv(th[sh])) else NA_real_
    mu <- exp(off + c(.fd_index(replace(x0, variable, from), th), .fd_index(replace(x0, variable, to), th)))
    m <- fam$meanfun(mu, tn)
    if (isTRUE(object$truncated)) m <- m / pmax(1 - fam$p0(mu, tn), 1e-12)   # E(Y | Y > 0)
    list(from = m[1L], to = m[2L])
  }
  c(list(component = "mean", lev = lev), .qoi_unc(.qoi_par(object)))
}

## the participation equation of a hurdle: its coefficients (aliased ones at zero) and covariance, prefixed
.fd_part <- function(object) {
  g <- object$part_beta; g[is.na(g)] <- 0
  V <- .part_vcov(object)
  nm <- paste0("participation:", names(g))
  if (!is.null(V)) dimnames(V) <- list(paste0("participation:", rownames(V)), paste0("participation:", colnames(V)))
  list(theta = stats::setNames(g, nm), vc = V)
}
## the covariance of an intensity fit's stored draws over `cols`, the coefficients renamed with `prefix`
.fd_draws_vc <- function(fit, cols, coefs, prefix) {
  D <- .ud_draws(fit, cols)
  if (is.null(D)) return(NULL)
  V <- stats::cov(D); nm <- ifelse(colnames(D) %in% coefs, paste0(prefix, colnames(D)), colnames(D))
  dimnames(V) <- list(nm, nm); V
}

.fd_spec.hurdle_count <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  fam <- .count_fam(object$family); linkinv <- stats::make.link(object$link)$linkinv
  P <- .fd_twopart_profiles(object$int_xref, object$part_xref, variable, from, to, at, stage, "intensity", "participation")
  base <- if (is.null(object$int_fe_ref)) 0 else object$int_fe_ref                 # the offset at its mean
  sh <- if (fam$nshape) paste0("intensity:", fam$shape_name) else character(0)
  pfun <- function(z, th) linkinv(.fd_index(z, th, "participation:"))
  mfun <- function(x, th) {
    tn <- if (fam$nshape) unname(fam$shape_inv(th[sh])) else NA_real_; lam <- exp(base + .fd_index(x, th, "intensity:"))
    fam$meanfun(lam, tn) / pmax(1 - fam$p0(lam, tn), 1e-12)                                    # E(Y | Y > 0)
  }
  c(list(component = c("participation", "intensity", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, `*`)),
    .qoi_unc(.qoi_par(object)))
}

.fd_spec.zi_count <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  fam <- .count_fam(object$family); linkinv <- stats::make.link(object$link)$linkinv
  P <- .fd_twopart_profiles(object$xref, object$zref, variable, from, to, at, stage, "count", "inflation")
  off <- .mean_offset(object); sh <- if (fam$nshape) fam$shape_name else character(0)
  pfun <- function(z, th) linkinv(.fd_index(z, th, "zero:"))
  mfun <- function(x, th) fam$meanfun(exp(off + .fd_index(x, th, "count:")),
                                      if (fam$nshape) unname(fam$shape_inv(th[sh])) else NA_real_)
  c(list(component = c("zero", "count", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, function(p, m) (1 - p) * m)),
    .qoi_unc(.qoi_par(object)))
}

.fd_spec.hurdle_cpb <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  int <- object$intensity; fe <- inherits(int, "cpb_fe"); sup <- .cpb_support(object)
  linkinv <- stats::make.link(if (is.null(object$link)) "logit" else object$link)$linkinv
  P <- .fd_twopart_profiles(object$int_xref, object$part_xref, variable, from, to, at, stage, "intensity", "participation")
  offm <- if (fe) object$int_fe_ref - .fe_ref(int) else object$int_fe_ref         # the offset at its mean
  pfun <- function(z, th) linkinv(.fd_index(z, th, "participation:"))
  mfun <- function(x, th)                                                          # E(Y | Y > 0)
    .cpb_mean(exp(offm + (if (fe) th[["fe_ref"]] else 0) + .fd_index(x, th, "intensity:")), th[["alpha"]], TRUE, sup)
  c(list(component = c("participation", "intensity", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, `*`)),
    .qoi_unc(.qoi_par(object)))
}

.fd_spec.hurdle_gec <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  int <- object$intensity; ms <- object$max.support
  linkinv <- stats::make.link(if (is.null(object$link)) "logit" else object$link)$linkinv
  P <- .fd_twopart_profiles(object$int_xref, object$part_xref, variable, from, to, at, stage, "intensity", "participation")
  pfun <- function(z, th) linkinv(.fd_index(z, th, "participation:"))
  mfun <- function(x, th) {                                                        # E(Y | Y > 0)
    mu <- exp(object$int_fe_ref + .fd_index(x, th, "intensity:")); d <- th[["delta"]]
    gec_mean_cpp(mu, d, ms) / pmax(1 - gec_pmf_cpp(mu, d, 0L, ms)[, 1], 1e-12)
  }
  c(list(component = c("participation", "intensity", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, `*`)),
    .qoi_unc(.qoi_par(object)))
}

.fd_spec.zi_cpb <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  fe <- !is.null(object$fe); sup <- .cpb_support(object)
  P <- .fd_twopart_profiles(object$int_xref, object$zero_xref, variable, from, to, at, stage, "count", "inflation")
  pfun <- function(z, th) stats::plogis(.fd_index(z, th, "zero:"))
  mfun <- function(x, th)
    .cpb_mean(exp((if (fe) th[["fe_ref"]] else object$int_fe_ref) + .fd_index(x, th, "count:")), th[["alpha"]], FALSE, sup)
  c(list(component = c("zero", "count", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, function(p, m) (1 - p) * m)),
    .qoi_unc(.qoi_par(object)))
}

.fd_spec.zi_gec <- function(object, variable, from, to, at = NULL, stage = "both", ...) {
  zref <- if (!is.null(object$zref)) object$zref else colMeans(object$Z)
  P <- .fd_twopart_profiles(colMeans(object$X), zref, variable, from, to, at, stage, "count", "inflation")
  off <- .mean_offset(object); ms <- object$max.support
  pfun <- function(z, th) stats::plogis(.fd_index(z, th, "zero:"))
  mfun <- function(x, th) gec_mean_cpp(exp(off + .fd_index(x, th, "count:")), th[["delta"]], ms)
  c(list(component = c("zero", "count", "marginal"), lev = .fd_twopart_lev(P, pfun, mfun, function(p, m) (1 - p) * m)),
    .qoi_unc(.qoi_par(object)))
}
## the zero-equation parameters a zero-inflated fit's draws do not carry (its zero_fe unit dummies)
.fd_zero_hold <- function(nm, draws)
  if (is.null(draws)) character(0) else setdiff(nm[startsWith(nm, "zero:")], colnames(draws))

## The refit bootstrap of a two-part first difference (B > 0): the model refit on resampled rows, or whole
## clusters (a unit drawn twice relabeled as two units in the columns `relabel`), through .ud_boot(); each
## replicate's first difference comes from its own specification, and the interval is the percentile one.
.fd_refit_boot <- function(pt, refit, spec_of, data, B, level, cores, cl = NULL, relabel = NULL) {
  a <- (1 - level) / 2
  clv <- if (is.null(cl)) NULL else data[[cl]]
  if (is.factor(clv)) clv <- droplevels(clv)                 # an unused level would be drawn as an empty cluster
  D <- .ud_boot(nrow(data), clv, function(idx, blk) {
    d <- data[idx, , drop = FALSE]
    if (!is.null(blk)) for (v in relabel) d[[v]] <- paste0(blk, "_", d[[v]])
    f2 <- .ud_quiet_guard(.ud_replicate(refit(d)))
    if (!isTRUE(f2$converged)) return(NULL)
    s2 <- spec_of(f2); L <- s2$lev(s2$theta); L$to - L$from
  }, pt$component, B, cores, label = if (is.null(cl)) "obs" else cl)
  ok <- D[stats::complete.cases(D), , drop = FALSE]
  if (nrow(ok) >= 2L) {
    pt$lower <- apply(ok, 2L, stats::quantile, a, names = FALSE)
    pt$upper <- apply(ok, 2L, stats::quantile, 1 - a, names = FALSE)
    pt$method <- sprintf("bootstrap (refit, B=%d)", B)
  }
  pt
}
