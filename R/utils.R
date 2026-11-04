## ---------------------------------------------------------------------------
## Internal utilities shared by the estimators: estimation-row bookkeeping,
## frequency weights, the parallel replicate driver used by every bootstrap.
## ---------------------------------------------------------------------------

## indices (into the data frame passed to model.frame) of the rows that
## survived na.omit -- independent of row names
.ud_kept_rows <- function(mf, data) {
  na <- attr(mf, "na.action")
  idx <- seq_len(nrow(data))
  if (is.null(na)) idx else idx[-as.integer(na)]
}

## resolve a `weights` argument (column name or numeric vector, aligned to
## `data`) to the estimation rows `rows`; NULL stays NULL
.ud_weights <- function(weights, data, rows) {
  if (is.null(weights)) return(NULL)
  w <- if (is.character(weights) && length(weights) == 1L) {
    if (!weights %in% names(data)) stop("'weights' must name a column of 'data' or be a numeric vector.")
    data[[weights]]
  } else weights
  if (!is.numeric(w)) stop("'weights' must be numeric (frequency weights).")
  if (length(w) != nrow(data)) stop("'weights' must have one value per row of 'data' (", nrow(data), ").")
  w <- as.numeric(w[rows])
  if (anyNA(w)) stop("'weights' has missing values on the estimation rows.")
  if (any(!is.finite(w)) || any(w < 0)) stop("'weights' must be finite and non-negative.")
  w
}

## weight vector for the likelihood: ones when no weights were supplied
.ud_w1 <- function(w, n) if (is.null(w)) rep(1, n) else w

## the frequency weights a calibration tool applies to a fit's rows: NULL for an
## unweighted fit, so the unweighted computations run exactly as before
.ud_fit_w <- function(fit, n) {
  w <- fit$weights
  if (!is.null(w) && length(w) == n) as.numeric(w) else NULL
}

## A count-valued argument (max.support, B, inner_it) is one whole number from
## `lower` to the largest integer: as.integer() turns a larger value into NA,
## which the C++ cannot take, and a negative or fractional value would reach the
## optimizer as a meaningless support. Refused by name before any fitting.
.ud_check_whole <- function(x, name, lower = 1, null_ok = FALSE) {
  if (null_ok && is.null(x)) return(invisible(NULL))
  top <- .Machine$integer.max - 1L                 # the C++ loops run one past it (max.support + 1)
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < lower || x > top)
    stop("'", name, "' must be ", if (null_ok) "NULL or ", "a whole number from ", lower, " to ", top, ".", call. = FALSE)
  invisible(NULL)
}

## The support guard's warning ("... reaches max.support ...") and King's CPB's
## equidispersed-limit warning muffled and every other warning passed on, for a
## fit that is a step inside another computation: a start value, a bootstrap or
## jackknife refit, the nest inside zi_test(). The fit the user asked for
## reports them itself.
.ud_quiet_guard <- function(expr)
  withCallingHandlers(expr, warning = function(w)
    if (grepl("reaches max.support", conditionMessage(w), fixed = TRUE) ||
        grepl("equidispersed limit", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))

## Parallel replicate driver for the bootstraps and the parametric-bootstrap
## screen threshold. `cores = 1` runs lapply(); `cores > 1` runs a PSOCK cluster
## (portable across Windows, macOS, and Linux) that loads this package from the
## calling session's library path and draws its random numbers from L'Ecuyer
## streams seeded from the calling session's RNG, so `set.seed()` before the
## call makes a parallel run reproducible (the parallel and serial streams
## differ, so a parallel result equals a serial one only in distribution).
.ud_lapply <- function(X, FUN, cores = 1L) {
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L || length(X) < 2L) return(lapply(X, FUN))
  cores <- min(cores, length(X))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ## the workers attach the package (a user's fitfun names cpb(), count_reg(), ... unqualified)
  ## the library the session loaded underdisp from comes first, so that a development build loaded through
  ## lib.loc is the one the workers run
  lp <- unique(c(dirname(system.file(package = "underdisp")), .libPaths()))
  parallel::clusterCall(cl, function(lp) { .libPaths(lp); suppressPackageStartupMessages(library("underdisp", character.only = TRUE)); NULL },
                        lp)
  seed <- sample.int(.Machine$integer.max, 1L)     # drawn here: the master's RNG advances
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::parLapply(cl, X, FUN)
}

## Column scales for the optimizer: each non-constant design column is divided
## by its (frequency-weighted) standard deviation before optimization and the
## coefficients mapped back, so the fit (and the numerical Hessians differenced
## in the scaled coordinates) do not depend on the units the covariates are
## measured in, and a weight of w and w copies of the row give the same scaled
## problem. Constant columns (the intercept, a column of zeros) keep scale 1.
.ud_colscale <- function(X, w = NULL) {
  wv <- if (is.null(w)) rep(1, nrow(X)) else as.numeric(w)
  sw <- sum(wv)
  s <- apply(X, 2, function(v) {
    m <- sum(wv * v) / sw
    r <- if (sw > 1) sqrt(sum(wv * (v - m)^2) / (sw - 1)) else NA_real_
    if (!is.finite(r) || r <= 1e-10 * max(1, max(abs(v)))) 1 else r
  })
  s[colnames(X) %in% "(Intercept)"] <- 1
  unname(s)
}

## rank check for a design matrix: refuse an aliased column (an empty factor
## level, an exact linear combination) with a message naming it, rather than
## letting a flat likelihood direction reach the optimizer and the bootstrap
.ud_rank_check <- function(X, what = "design matrix") {
  if (ncol(X) == 0L) return(invisible(TRUE))
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    stop("the ", what, " is rank deficient: column(s) ", paste(bad, collapse = ", "),
         " are collinear with the others (an empty factor level or an exact linear combination); ",
         "drop them or use droplevels().", call. = FALSE)
  }
  invisible(TRUE)
}

## within-unit rank check for the concentrated fixed-effects fits: a covariate
## constant within units is absorbed by the unit intercepts
.ud_rank_check_fe <- function(X, unit) {
  if (ncol(X) == 0L) return(invisible(TRUE))
  Xd <- X - apply(X, 2, function(v) stats::ave(v, unit))
  wvar <- apply(Xd, 2, function(v) max(abs(v)))
  const <- colnames(X)[wvar <= 1e-10 * pmax(1, apply(abs(X), 2, max))]
  if (length(const))
    stop("covariate(s) ", paste(const, collapse = ", "), " do not vary within units and are absorbed by ",
         "the unit fixed effects; remove them from the formula.", call. = FALSE)
  qx <- qr(Xd)
  if (qx$rank < ncol(X)) {
    bad <- colnames(X)[qx$pivot[(qx$rank + 1L):ncol(X)]]
    stop("the within-unit design is rank deficient: column(s) ", paste(bad, collapse = ", "),
         " are collinear with the others after removing the unit effects.", call. = FALSE)
  }
  invisible(TRUE)
}

## A gradient-based pass (BFGS) from a multistart's best solution, kept only if it
## improves on it: Nelder-Mead stops on a flat or collapsed simplex short of the
## optimum, and BFGS's line search accepts only decreases, so the pass cannot make
## the fit worse. `fn` may return non-finite values at infeasible points.
.ud_bfgs_pass <- function(fit, fn, maxit = 500L, reltol = 1e-10) {
  f <- function(p) { v <- fn(p); if (is.finite(v)) v else 1e10 }
  op <- tryCatch(stats::optim(fit$par, f, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
                 error = function(e) NULL)
  if (!is.null(op) && op$value < fit$value) op else fit
}

## Nelder-Mead restart until the improvement stops paying (a single run can
## stall on a collapsed simplex short of the optimum)
.ud_nm_polish <- function(fit, fn, maxit, reltol, restarts = 3L, tol = 1e-6) {
  if (length(fit$par) == 1L) {                       # Nelder-Mead is unreliable in one dimension
    op <- tryCatch(stats::optimize(fn, c(fit$par - 1, fit$par + 1), tol = 1e-10), error = function(e) NULL)
    if (!is.null(op) && is.finite(op$objective) && op$objective < fit$value) { fit$par <- op$minimum; fit$value <- op$objective }
    return(fit)
  }
  for (rs in seq_len(restarts)) {
    op2 <- tryCatch(stats::optim(fit$par, fn, method = "Nelder-Mead",
                                 control = list(maxit = maxit, reltol = reltol)), error = function(e) NULL)
    if (is.null(op2)) break
    gain <- fit$value - op2$value                    # improvement over the incumbent
    if (gain > 0) fit <- op2
    if (gain < tol) break                            # the chain stops when a restart gains less than tol
  }
  fit
}

## `cluster` resolved to a character vector on the estimation rows: a column
## name in `data`, or a vector with one value per row of `data` (checked before
## the rows are reduced, so a misaligned vector is refused rather than recycled
## or padded); missing values on the estimation rows are refused.
.ud_cluster_values <- function(cluster, data, keep = NULL) {
  if (is.null(cluster)) return(NULL)
  cv <- if (is.character(cluster) && length(cluster) == 1L) {
    if (is.null(data[[cluster]])) stop("'cluster' must be a column name in 'data' or a vector.", call. = FALSE)
    data[[cluster]]
  } else cluster
  if (length(cv) != nrow(data)) stop("'cluster' must have one value per row of 'data'.", call. = FALSE)
  if (!is.null(keep)) cv <- cv[keep]
  if (anyNA(cv)) stop("'cluster' has missing values on the estimation rows.", call. = FALSE)
  as.character(cv)
}

## guard for a cluster bootstrap: no resampling variation with one cluster
.ud_cluster_guard <- function(clid) {
  ng <- length(unique(clid))
  if (ng < 2L) stop("the cluster bootstrap needs at least two clusters; 'cluster' is constant on the estimation rows.", call. = FALSE)
  if (ng < 30L) warning(sprintf("only %d clusters; the cluster bootstrap may be unreliable with so few.", ng), call. = FALSE)
  invisible(ng)
}

## first differences hold every other model-matrix column at its reference
## value; a variable that also enters another column (an interaction, I(x^2),
## log(x), poly(x, 2), a spline basis) cannot be moved through one column alone.
## The variable is matched as a whole name inside the other column names.
.fd_check_terms <- function(variable, cols) {
  esc <- gsub("([^[:alnum:]_])", "\\\\\\1", variable, perl = TRUE)      # every other character literal
  other <- setdiff(unique(cols), variable)
  hit <- other[grepl(paste0("(^|[^[:alnum:]._])", esc, "([^[:alnum:]._]|$)"), other, perl = TRUE)]
  if (length(hit))
    stop("'", variable, "' also enters the term(s) ", paste(hit, collapse = ", "),
         " (an interaction or a transformation of it); first_difference() moves one model-matrix column and ",
         "cannot recompute the others. ",
         "Use predict(newdata = ) with the two full covariate profiles instead.", call. = FALSE)
  invisible(TRUE)
}

## Outer optimizer for the concentrated fixed-effects likelihoods (cpb_fe, gec_fe).
## `nll_raw(par, awarm)` returns list(nll, unit intercepts, edge rows) and
## `grad_raw(par, a, edge)` the gradient of the concentrated objective at those
## intercepts (the envelope theorem with the breakpoint correction; see
## src/cpb_fe.cpp). The closure caches the last evaluation for the gradient.
## Each dispersion start runs BFGS from the Poisson slopes (Nelder-Mead if BFGS
## fails); the incumbent is polished by Nelder-Mead and challenged from two
## perturbed restarts. The returned `$a` are the intercepts that attain `$value`.
.fe_outer <- function(nll_raw, grad_raw, bs, dstarts, maxit, reltol) {
  ## every evaluation searches each unit from cold: a warm start carried from a
  ## distant trial point can leave a unit on the wrong tooth, and the resulting
  ## path dependence defeats the gradient methods
  env <- new.env(); env$last <- NULL
  fn <- function(par) {
    r <- nll_raw(par, numeric(0)); v <- r[[1L]]
    if (!is.finite(v)) v <- 1e10                       # an overflowed trial point is infeasible
    if (v < 1e9) env$last <- list(par = as.numeric(par), a = r[[2L]], edge = r[[3L]])
    v
  }
  gr <- function(par) {
    if (is.null(env$last) || !identical(as.numeric(par), env$last$par)) fn(par)
    if (is.null(env$last) || !identical(as.numeric(par), env$last$par)) return(rep(0, length(par)))
    g <- grad_raw(par, env$last$a, env$last$edge)
    if (all(is.finite(g))) return(g)
    ## a non-finite analytic gradient (an edge case of the breakpoint correction): central differences of the
    ## objective instead, and NA, which stops BFGS for the Nelder-Mead fallback, if those are not finite either
    h <- 1e-6 * pmax(1, abs(par))
    g <- vapply(seq_along(par), function(k) {
      e <- replace(numeric(length(par)), k, h[k]); fp <- fn(par + e); fm <- fn(par - e)
      if (fp >= 1e9 || fm >= 1e9) NA_real_ else (fp - fm) / (2 * h[k]) }, numeric(1))   # an infeasible side: NA
    fn(par)                                              # leave the cache at `par`
    if (all(is.finite(g))) g else rep(NA_real_, length(par))
  }
  run <- function(start) {
    o <- tryCatch(stats::optim(start, fn, gr, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$value) || o$value >= 1e9)
      o <- tryCatch(stats::optim(start, fn, method = "Nelder-Mead", control = list(maxit = maxit, reltol = reltol)),
                    error = function(e) NULL)
    o
  }
  cand <- lapply(dstarts, function(d0) run(c(bs, d0)))
  cand <- cand[!vapply(cand, is.null, logical(1))]
  if (!length(cand)) return(NULL)
  best <- cand[[which.min(vapply(cand, function(f) f$value, numeric(1)))]]
  if (best$value >= 1e9) return(best)
  best <- .ud_nm_polish(best, fn, maxit, reltol, restarts = 1L)
  ## the concentrated objective can jump in the slopes where a unit sits on its
  ## feasibility floor (see src/cpb_fe.cpp), so the incumbent may rest against a
  ## cliff: restart from two deterministic perturbations and keep the best
  p <- length(best$par)
  for (dlt in list(rep(0.1, p), -0.1 * (-1)^seq_len(p))) {
    o2 <- run(best$par + dlt)
    if (!is.null(o2) && is.finite(o2$value) && o2$value < best$value - 1e-8)
      best <- .ud_nm_polish(o2, fn, maxit, reltol, restarts = 1L)
  }
  rc <- nll_raw(best$par, numeric(0)); best$value <- rc[[1L]]; best$a <- rc[[2L]]
  best
}

## offset() terms in a formula are refused: every estimator takes the exposure
## through `offset =` (log scale), and a formula term that model.matrix() drops
## would otherwise be ignored silently
.ud_no_formula_offset <- function(...) {
  for (f in list(...)) {
    if (!inherits(f, "formula")) next
    if (!is.null(attr(stats::terms(f, allowDotAsName = TRUE), "offset")))
      stop("offset() terms in a formula are not supported; supply the exposure on the log scale ",
           "through offset = (a column name or a vector).", call. = FALSE)
    rhs <- f[[length(f)]]
    if (is.call(rhs) && identical(rhs[[1L]], as.name("|")))
      stop("a two-part formula 'y ~ x | z' is not supported; give the second equation through ",
           "participation = ~ z (hurdle) or zero = ~ z (zero-inflated).", call. = FALSE)
  }
  invisible(TRUE)
}

## The offset of new rows (predict(), score(), ceiling_exceedance()). An `offset`
## given to the method (a column of newdata or a vector) is used as given.
## Otherwise a fit made with an offset reads it from the same column of
## newdata, and refuses when it was given as a vector or newdata lacks the
## column, rather than predicting the new rows at exposure 1. A fit without an
## offset uses 0. `object` is the fit (or the two-part fit's count component)
## that stores `offset` (its vector) and `offset_name` (its column, if any).
.ud_newdata_offset <- function(object, newdata, offset = NULL, n = nrow(newdata)) {
  if (!is.null(offset))
    return(rep_len(as.numeric(if (is.character(offset) && length(offset) == 1L) newdata[[offset]] else offset), n))
  o <- object[["offset"]]
  if (is.null(o) || !length(o) || all(o == 0)) return(rep(0, n))
  nm <- object[["offset_name"]]
  if (!is.null(nm) && nm %in% names(newdata)) return(as.numeric(newdata[[nm]]))
  stop("The model was fitted with an offset; give the new rows' offset through offset = (a column of newdata ",
       "or a vector)", if (!is.null(nm)) paste0(", or include the column '", nm, "' in newdata"), ".", call. = FALSE)
}
## A count formula with its unit dummies added (fe = NULL: unchanged but for a
## '.' expanded). The response keeps its own expression (a transformed response
## stays transformed), the intercept and the environment are kept, and a '.'
## expands over the data's columns other than the unit identifiers.
.ud_add_fe_terms <- function(formula, fe, data) {
  tf <- stats::terms(formula, data = data[, setdiff(names(data), fe), drop = FALSE])
  if (is.null(fe) && !any(all.names(formula) == ".")) return(formula)
  stats::reformulate(c(labels(tf), if (!is.null(fe)) paste0("factor(", fe, ")")), response = formula[[2L]],
                     intercept = attr(tf, "intercept") == 1L, env = environment(formula))
}

## A fixed-effects fit estimates on rows sorted by unit (each unit's rows
## contiguous for the concentrated likelihood); its per-row results go back to
## the data's own row order, so that fitted(), residuals(), predict() and
## simulate() line up with the rows supplied. `perm[i]` is the data position of
## sorted row i; `fields` name the per-row vectors and matrices of the fit.
.ud_unsort <- function(fit, perm, fields) {
  inv <- order(perm)
  for (f in fields) {
    v <- fit[[f]]
    if (is.matrix(v) && nrow(v) == length(perm)) fit[[f]] <- v[inv, , drop = FALSE]
    else if (!is.null(v) && !is.matrix(v) && length(v) == length(perm)) fit[[f]] <- v[inv]
  }
  fit
}

## the column name an `offset` argument gives, kept with the fit for new data
.ud_offset_name <- function(offset) if (is.character(offset) && length(offset) == 1L) offset

## The units whose fixed effect the data identify. A unit with no positive count
## has its effect at minus infinity (the search stops at a sentinel), and under
## King's zero-truncated support a unit whose positive counts are all 1 fits
## every intercept below the one that puts its ceiling at 1 equally well (its
## truncated support is then {1}). Their values are arbitrary, so they stay out
## of the reference unit that new rows and first differences use (`fe_ref`).
## On the floor support the zero-truncated flat set is bounded and is kept.
.fe_identified <- function(Y, unit, nu, truncated = FALSE, support = "ceiling") {
  ymax <- vapply(split(Y, factor(unit, levels = seq_len(nu))), function(v) if (length(v)) max(v) else 0, numeric(1))
  if (!truncated) ymax > 0 else if (.cpb_fl(support)) rep(TRUE, nu) else ymax >= 2
}
.fe_mean_identified <- function(fe_hat, ident) if (any(ident)) mean(fe_hat[ident]) else mean(fe_hat)
## The same rule for fixed effects entered as dummy columns (count_reg and the
## two-part count models): in a first difference the dummies of each `fe`
## variable are held at their means over the rows of the identified units, as
## judged by `ident()` on a unit's responses, so that the diverging
## coefficients of the others do not move the reference unit.
.fe_dummy_ref <- function(xref, X, y, fe, ident) {
  for (v in fe) {
    cols <- colnames(X)[startsWith(colnames(X), paste0("factor(", v, ")"))]
    if (!length(cols)) next
    lev <- as.integer(X[, cols, drop = FALSE] %*% seq_along(cols))        # 0: the reference level
    ok <- vapply(split(y, lev), ident, logical(1))
    rows <- lev %in% as.integer(names(ok)[ok])
    if (any(rows)) xref[cols] <- colMeans(X[rows, cols, drop = FALSE])
  }
  xref
}
.ident_count  <- function(y) any(y > 0)                       # a count equation
.ident_trunc  <- function(y) any(y >= 2)                      # a zero-truncated count equation
.ident_binary <- function(y) any(y > 0) && any(y == 0)        # a participation or inflation equation

## the reference unit effect of a concentrated fixed-effects fit (an object from before 0.2.0: the mean effect)
.fe_ref <- function(object) if (!is.null(object$fe_ref)) object$fe_ref else mean(object$fe)
## The unit effect of each new row of a concentrated fixed-effects fit: a unit the fit estimated (named in the
## `fe_var` column of `newdata`) at its own effect, an unseen unit, or every row when `newdata` has no such column,
## at the reference unit `ref`
.fe_newdata_effect <- function(newdata, fe_var, effects, ref) {
  if (is.null(fe_var) || !fe_var %in% names(newdata)) return(rep(ref, nrow(newdata)))
  e <- unname(effects[as.character(newdata[[fe_var]])]); e[is.na(e)] <- ref
  e
}

## rows whose unit identifier is missing are dropped, as rows with a missing
## response or covariate are; `fe` may name several columns or be NULL
.ud_drop_na_fe <- function(data, fe) {
  fe <- fe[!vapply(fe, is.null, logical(1))]
  fe <- unlist(fe); fe <- fe[is.character(fe) & fe %in% names(data)]
  if (!length(fe)) return(data)
  bad <- Reduce(`|`, lapply(fe, function(v) is.na(data[[v]])))
  if (any(bad)) { data <- data[!bad, , drop = FALSE]; attr(data, "ud_kept") <- !bad }
  data
}
## a vector `offset`, `weights`, or `cluster` given for the original rows of a
## data set that .ud_drop_na_fe() reduced follows the reduction
.ud_align_vec <- function(v, data) {
  kept <- attr(data, "ud_kept")
  if (is.null(v) || is.null(kept) || (is.character(v) && length(v) == 1L) || length(v) != length(kept)) return(v)
  v[kept]
}

## The model matrix of `newdata` in the fit's own coding: the stored factor
## levels and contrasts (a factor's contrasts attribute is lost when
## model.frame() re-levels it through xlev, so a data set whose reference
## level was set through that attribute would otherwise be re-coded), with the
## columns selected by the coefficient names, so that a coding mismatch is
## refused rather than multiplied into the coefficients by position. A newdata
## factor's own contrasts attribute is dropped where the stored coding replaces
## it (model.frame() would otherwise warn that it drops it).
.ud_newdata_matrix <- function(terms_obj, newdata, xlev = NULL, contrasts = NULL, coefnames = NULL) {
  terms_obj <- stats::delete.response(terms_obj)
  miss <- setdiff(all.vars(terms_obj), names(newdata))
  if (length(miss)) stop("'newdata' is missing required variable(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  for (nm in intersect(names(contrasts), names(newdata)))
    if (is.factor(newdata[[nm]]) && !is.null(attr(newdata[[nm]], "contrasts"))) attr(newdata[[nm]], "contrasts") <- NULL
  ## a column holding only missing values reads as logical (data.frame(x = NA)); it is a missing value of
  ## the fitted variable's own type, so its rows predict NA like any other missing covariate
  for (nm in intersect(all.vars(terms_obj), names(newdata)))
    if (is.logical(newdata[[nm]]) && all(is.na(newdata[[nm]])))
      newdata[[nm]] <- if (nm %in% names(xlev)) factor(newdata[[nm]], levels = xlev[[nm]]) else as.numeric(newdata[[nm]])
  mf <- stats::model.frame(terms_obj, newdata, xlev = xlev, na.action = stats::na.pass)
  X <- stats::model.matrix(terms_obj, mf, contrasts.arg = contrasts)
  if (!is.null(coefnames)) {
    lack <- setdiff(coefnames, colnames(X))
    if (length(lack)) stop("'newdata' does not reproduce the fit's design columns (", paste(lack, collapse = ", "),
                           "): each variable must have the type it had in the fitting data, and factors must ",
                           "carry its levels and coding.", call. = FALSE)
    X <- X[, coefnames, drop = FALSE]
  }
  X
}

## Rows of newdata with a missing covariate or offset predict NA, as in predict.glm(). The distribution
## routines take no missing rates, so those rows are evaluated at a placeholder and masked afterwards.
.ud_na_rows <- function(X, off = NULL) {
  bad <- !stats::complete.cases(X)
  if (!is.null(off)) bad <- bad | is.na(rep_len(off, nrow(X)))
  bad
}
.ud_mask <- function(x, bad) { if (any(bad)) x[bad] <- NA; x }

## A response without a single positive count leaves the rate at its lower bound and the dispersion
## parameter without information; the fit still returns, so it says so.
.ud_warn_all_zero <- function(y) {
  if (length(y) && all(y == 0))
    warning("every count is zero: the rate runs to its lower bound, so the coefficients diverge and the ",
            "dispersion parameter is not identified by these data.", call. = FALSE)
  invisible(NULL)
}

## Predicted probabilities of a participation/inflation glm on `newdata`, in the
## glm's own factor coding: predict.lm() re-levels each factor through the
## stored xlevels and applies the stored contrasts, so a newdata factor's own
## contrasts attribute is dropped first (predict.lm() would otherwise warn that
## it drops it).
.ud_glm_predict <- function(object, newdata) {
  if (is.null(newdata)) return(as.numeric(stats::predict(object, type = "response")))
  for (nm in intersect(names(object$xlevels), names(newdata)))
    if (is.factor(newdata[[nm]]) && !is.null(attr(newdata[[nm]], "contrasts"))) attr(newdata[[nm]], "contrasts") <- NULL
  as.numeric(stats::predict(object, newdata = newdata, type = "response"))
}

