## ---------------------------------------------------------------------------
## Simulators for the CPB family
## ---------------------------------------------------------------------------

#' Simulate from the continuous parameter binomial
#'
#' Draws counts from the continuous parameter binomial (CPB) or its
#' zero-truncated variant.
#'
#' @param n Number of draws (recycled against `lambda`).
#' @param lambda Rate parameter; a scalar or a length-`n` vector (see [dcpb()]).
#' @param alpha Shape/dispersion parameter in (0, 1).
#' @param truncated If `TRUE`, draw from the zero-truncated CPB.
#' @param support The CPB's support, `"ceiling"` (King's CPB, the default) or
#'   `"floor"`; see [cpb()].
#' @details King's CPB has the support `0, ..., ceiling(lambda / (1 - alpha))`,
#' which contains 1 for every positive rate. The floor support,
#' `0, ..., floor(lambda / (1 - alpha))`, is `{0}` for a rate whose ceiling is
#' below 1, so its zero-truncated distribution does not exist. With
#' `truncated = TRUE` such a draw is set to 1 and a warning reports how many
#' there were: the stated parameters give that count probability zero, so data
#' simulated this way are not data from the model (in a simulation study, keep
#' every ceiling at 1 or above).
#' @return An integer vector of counts.
#' @examples
#' set.seed(1)
#' table(rcpb(1000, lambda = 3, alpha = 0.5))
#' @export
rcpb <- function(n, lambda, alpha, truncated = FALSE, support = c("ceiling", "floor")) {
  support <- match.arg(support)
  if (length(alpha) != 1L || alpha <= 0 || alpha >= 1) stop("`alpha` must be a single value in (0, 1).")
  if (n == 0) return(integer(0))
  lambda <- rep_len(lambda, n)
  out <- integer(n); nbelow <- 0L
  for (i in seq_len(n)) {
    pr <- .cpb_pmf_core(lambda[i], alpha, support)
    if (anyNA(pr)) { out[i] <- NA_integer_; next }           # refused at the cap
    K  <- length(pr) - 1L
    if (truncated) {
      if (K < 1) { out[i] <- 1L; nbelow <- nbelow + 1L; next }
      pr[1] <- 0; pr <- pr / sum(pr)
    }
    out[i] <- as.integer(sample(0:K, 1L, prob = pr))
  }
  if (nbelow > 0L)
    warning(nbelow, " of ", n, " rates have a ceiling lambda / (1 - alpha) below 1, where the zero-truncated CPB does not ",
            "exist; those draws are set to 1, a count the stated parameters give probability zero.", call. = FALSE)
  .cpb_warn_cap(out)
}

#' Simulate from the hurdle continuous parameter binomial
#'
#' @param n Number of units.
#' @param lambda Intensity rate for participants (the CPB's `lambda`); scalar or length-`n`.
#' @param alpha CPB shape parameter in (0, 1).
#' @param p Participation probability; scalar or length-`n`.
#' @param support The CPB's support, `"ceiling"` (King's CPB, the default) or
#'   `"floor"`; see [cpb()].
#' @return An integer vector: 0 for nonparticipants, a zero-truncated CPB draw
#'   otherwise.
#' @examples
#' set.seed(1)
#' y <- rhurdle_cpb(1000, lambda = 3, alpha = 0.5, p = 0.6)
#' mean(y == 0)
#' @export
rhurdle_cpb <- function(n, lambda, alpha, p, support = c("ceiling", "floor")) {
  support <- match.arg(support)
  d <- stats::rbinom(n, 1, p)
  y <- integer(n)
  idx <- which(d == 1)
  if (length(idx)) y[idx] <- rcpb(length(idx), rep_len(lambda, n)[idx], alpha, truncated = TRUE, support = support)
  y
}

## ---------------------------------------------------------------------------
## Hurdle-CPB estimator
## ---------------------------------------------------------------------------

#' Hurdle continuous parameter binomial regression
#'
#' Fits a participation (hurdle) model joined to a zero-truncated CPB for the
#' positive counts. Because the hurdle log-likelihood factorizes, the two parts
#' are fit separately: a logistic regression of participation on all units, and
#' a zero-truncated CPB ([cpb()]) on the positive counts. The result is the
#' natural model for a bounded, underdispersed participation process --- most
#' units at zero, participants carrying a tight count.
#'
#' @param formula Intensity model formula (`y ~ x`).
#' @param data A data frame.
#' @param participation Optional one-sided formula (`~ z`) for the participation
#'   (hurdle) model; defaults to the intensity model's right-hand side.
#' @param fe Optional column name for unit fixed effects in the intensity model,
#'   absorbed by a concentrated likelihood ([cpb_fe()]). This is how the
#'   within-unit underdispersion is recovered.
#' @param part_fe Optional column name for fixed effects in the participation
#'   model (added as factors). Units with no within-unit variation in
#'   participation are uninformative for a fixed-effects logit, and the dummy
#'   logit's coefficients carry the incidental-parameters bias of order 1/T (in
#'   short panels the slopes are inflated, roughly toward twice their value at
#'   T = 2); the split-panel jackknife (`bias_correct`) corrects only the
#'   intensity. The hurdle's likelihood separates its two parts, so for the
#'   logit link a conditional logit on the same data
#'   (\code{survival::clogit()}) gives participation slopes free of that bias,
#'   and an analytical correction (the \pkg{bife} package, for the logit and
#'   the probit) reduces it to second order; the hurdle's predictions and
#'   quantities of interest use the dummy logit's unit effects. On King's
#'   support the dummy of a unit that never participates, at minus infinity, is
#'   not counted in the degrees of freedom (see [compare_models()], Degrees of
#'   freedom). In new data ([predict.hurdle_cpb()], [score()], [cv_score()]) such
#'   a unit's rows get the participation probability at that limit, zero, and
#'   those of a unit whose every count is positive (`part_one_units`) one; a unit
#'   the fit has not seen, or every row of `newdata` without the `part_fe`
#'   column, enters at the reference unit of the units the data identify. The
#'   participation equation is checked for separation by [zi_cpb()]'s rule
#'   (Separation; component `part_separation`). On the floor support these
#'   dummies are counted in the degrees of freedom, new rows of these units are
#'   evaluated at the dummies' fitted values, and an unseen unit enters at the
#'   mean of all the unit dummies.
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`. Positive participation coefficients raise the probability of
#'   a positive count.
#' @param offset Optional offset on the log-mean scale for the intensity (an
#'   exposure): a numeric vector or the name of a column in `data`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector). The participation equation's standard errors are then
#'   cluster-robust (the sandwich, with the G/(G - 1) factor), and with
#'   `se = "bootstrap"` the intensity's bootstrap resamples whole clusters; see
#'   [cpb()].
#' @param se Inference for the intensity coefficients: `"none"` or `"bootstrap"`.
#' @param B Bootstrap replicates when `se = "bootstrap"`.
#' @param weights Optional frequency weights (a numeric vector or a column name),
#'   applied to both margins; see [cpb()].
#' @param cores Worker processes for the intensity bootstrap; see [cpb()].
#' @param ... Passed to [cpb()] or [cpb_fe()], including `support` (King's CPB,
#'   `"ceiling"`, by default; see [cpb()]).
#' @return An object of class `"hurdle_cpb"` with elements `participation` (a
#'   `glm`), `intensity` (a `cpb` or `cpb_fe`), and bookkeeping.
#' @examples
#' set.seed(1)
#' n <- 800; x <- rnorm(n); z <- rnorm(n)
#' y <- rhurdle_cpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5,
#'                  p = plogis(-0.2 + 0.8 * z))
#' fit <- hurdle_cpb(y ~ x, data = data.frame(y = y, x = x, z = z),
#'                   participation = ~ z)
#' fit
#' @export
hurdle_cpb <- function(formula, data, participation = NULL, fe = NULL, part_fe = NULL,
                       link = c("logit", "probit", "cloglog"), offset = NULL,
                       cluster = NULL, se = c("none", "bootstrap"), B = 500, weights = NULL,
                       cores = 1L, ...) {
  .ud_no_formula_offset(formula, participation)
  se <- match.arg(se); link <- match.arg(link)
  refit_args <- list(fe = fe, part_fe = part_fe, link = link, offset = offset, weights = weights,
                     cluster = if (is.character(cluster) && length(cluster) == 1L) cluster, dots = list(...))
  data <- .ud_drop_na_fe(data, list(fe, part_fe))
  offset <- .ud_align_vec(offset, data); weights <- .ud_align_vec(weights, data); cluster <- .ud_align_vec(cluster, data)
  ## reduce to complete cases on all model variables so the two margins stay aligned
  mv <- unique(c(all.vars(formula), all.vars(if (is.null(participation)) formula[-2L] else participation), fe, part_fe))
  keep <- stats::complete.cases(data[, intersect(mv, names(data)), drop = FALSE])
  cl_vals <- .ud_cluster_values(cluster, data, keep)           # checked before the rows are reduced
  if (!all(keep)) {
    if (!is.null(offset) && !is.character(offset)) offset <- offset[keep]
    if (!is.null(weights) && !is.character(weights)) weights <- weights[keep]
    data <- data[keep, , drop = FALSE]
  }
  y <- stats::model.response(stats::model.frame(formula, data))
  if (!is.numeric(y)) stop("Response must be a numeric count; got ", class(y)[1L], ".")
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  .ud_warn_all_zero(y)
  d <- as.integer(y > 0)
  int_offset <- if (is.null(offset)) NULL                       # offset applies to the intensity
                else if (is.character(offset) && length(offset) == 1L) offset else offset[y > 0]
  ## carry a cluster identifier as a column so both margins can see it aligned
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) {
    data[[".cluster"]] <- cl_vals; cluster <- ".cluster"
  }
  ## frequency weights as a column so both margins see them aligned
  if (!is.null(weights) && !(is.character(weights) && length(weights) == 1L)) {
    if (length(weights) != nrow(data)) stop("'weights' must have one value per row of 'data'.")
    data[[".hw"]] <- as.numeric(weights); weights <- ".hw"
  }
  w_all <- .ud_weights(weights, data, seq_len(nrow(data)))

  ## participation (hurdle) model: logit of d, optionally with fixed effects
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d; pdata[[".w"]] <- if (is.null(w_all)) rep(1, nrow(data)) else w_all
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, weights = .w,
                     family = stats::binomial(link = link))
  pf <- .part_fields(pfit, part_fe, data, if (is.null(cluster)) "model" else "cluster",
                     if (!is.null(cluster)) data[[cluster]], "hurdle_cpb")

  ## intensity model on the positive counts, with optional unit fixed effects (the cluster reaches it only for its
  ## bootstrap)
  pos <- data[y > 0, , drop = FALSE]; icl <- if (se == "bootstrap") cluster
  ifit <- if (is.null(fe)) cpb(formula, data = pos, truncated = TRUE, se = se, B = B,
                               cluster = icl, offset = int_offset, weights = weights, cores = cores, ...)
          else cpb_fe(formula, data = pos, fe = fe, truncated = TRUE, se = se, B = B,
                      cluster = icl, offset = int_offset, weights = weights, cores = cores, ...)

  ## per-observation intensity lambda (fixed-effects-aware) and a reference
  ## profile, for calibration and the first-difference decomposition
  Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
  beta <- ifit$coefficients
  off_all <- if (is.null(offset)) 0 else if (is.character(offset)) as.numeric(data[[offset]]) else as.numeric(offset)
  if (is.null(fe)) {
    lambda_full <- as.numeric(exp(off_all + Xall %*% beta))
    int_xref <- colMeans(Xall); int_fe_ref <- mean(off_all)
  } else {
    Xcov <- Xall[, names(beta), drop = FALSE]
    int_fe_ref <- .fe_ref(ifit) + mean(off_all)
    fe_i <- ifit$fe[as.character(data[[fe]])]; fe_i[is.na(fe_i)] <- .fe_ref(ifit)
    lambda_full <- as.numeric(exp(off_all + fe_i + Xcov %*% beta))
    int_xref <- colMeans(Xcov)
  }
  ## a unit that never participates has its participation dummy at minus infinity; not counted (.ud_zero_levels)
  zl <- if (is.null(part_fe) || .cpb_fl(.cpb_support(ifit))) character(0) else .ud_fe_zero_levels(y, data, part_fe)
  ## one whose every count is positive has it at plus infinity; counted, but predicted at that limit in new data
  ol <- if (is.null(part_fe) || .cpb_fl(.cpb_support(ifit))) character(0) else .ud_fe_one_levels(y, data, part_fe)

  structure(c(list(participation = pfit, intensity = ifit,
                 converged = isTRUE(pfit$converged) && isTRUE(ifit$converged), Zpart = stats::model.matrix(pfit),
                 y = as.integer(y), p_full = as.numeric(stats::fitted(pfit)),
                 lambda_full = lambda_full, n = length(y), n_participate = sum(d),
                 int_beta = beta, int_xref = int_xref, int_fe_ref = int_fe_ref,
                 part_beta = stats::coef(pfit),
                 part_xref = .fe_dummy_ref(colMeans(stats::model.matrix(pfit)), stats::model.matrix(pfit), d, part_fe,
                                           .ident_binary),
                 fe = fe, part_fe = part_fe, link = link, support = .cpb_support(ifit),
                 weights = w_all, nobs_weighted = if (is.null(w_all)) length(y) else sum(w_all),
                 cluster = if (is.character(cluster)) cluster else NULL,
                 loglik = as.numeric(stats::logLik(pfit)) + ifit$loglik,
                 df = length(stats::coef(pfit)) - length(zl) + ifit$df, part_zero_units = zl, part_one_units = ol,
                 se.type = se, formula = formula, part_formula = participation, refit_args = refit_args, call = match.call()), pf),
            class = "hurdle_cpb")
}

#' @method print hurdle_cpb
#' @export
print.hurdle_cpb <- function(x, ...) {
  cat("Hurdle Continuous Parameter Binomial")
  if (!is.null(x$fe)) cat(" (intensity fixed effects on '", x$fe, "')", sep = "")
  cat("\nCall:  ", deparse(x$call), "\n", sep = "")
  cat(sprintf("Units: %d (%d participate, %.0f%%)\n",
              x$n, x$n_participate, 100 * x$n_participate / x$n))
  pc <- stats::coef(x$participation)
  pc <- pc[!.is_fe_col(names(pc), x$part_fe)]   # hide participation fixed-effect dummies
  lk <- if (is.null(x$link)) "logit" else x$link
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n", lk))
  print(round(pc, 4))
  cat("\nIntensity (zero-truncated CPB) coefficients:\n"); print(round(x$intensity$coefficients, 4))
  cat("Intensity alpha (shape):", round(x$intensity$alpha, 4), "\n")
  .cpb_lower_note(x$intensity$alpha)
  for (nt in .sep_note(x$part_separation, "participation")) cat(nt, "\n")
  .ud_zero_df_note(x$part_zero_units, "participation dummy", "participation dummies")
  invisible(x)
}

#' Predict from a hurdle-CPB fit
#'
#' @param object A `"hurdle_cpb"` object.
#' @param newdata Data frame of covariate profiles (the intensity and participation
#'   covariates). A row of a fixed-effects unit the fit has not seen, or every row
#'   when the `fe` or `part_fe` column is absent, enters that equation at the
#'   reference unit, and a unit at a limit (no positive count, or no zero) at
#'   that limit (see `part_fe` in [hurdle_cpb()]).
#' @param type `"response"` for the marginal expected count E(Y) (zeros
#'   included), `"participation"` for P(Y > 0), or `"intensity"` for
#'   E(Y | Y > 0).
#' @param offset Optional offset (log scale) of the intensity for the `newdata`
#'   rows: a numeric vector or a column name in `newdata`; by default the fit's
#'   own offset column is read from `newdata` (see [predict.cpb()]).
#' @param se.fit,interval,level Standard errors (`se.fit = TRUE`) and a
#'   `level` confidence interval (`interval = "confidence"`) of the predicted
#'   quantity at the rows of `newdata`; see [predict.cpb()].
#' @param ... Unused.
#' @return A numeric vector. With `interval = "confidence"`, a matrix with
#'   columns `fit`, `lwr` and `upr`; with `se.fit = TRUE`, a list of `fit` (the
#'   vector or the matrix), `se.fit` and `method` (each row's interval source).
#'   Without `se.fit`, the source is the attribute `method`.
#' @method predict hurdle_cpb
#' @export
predict.hurdle_cpb <- function(object, newdata = NULL,
                               type = c("response", "participation", "intensity"), offset = NULL,
                               se.fit = FALSE, interval = c("none", "confidence"),
                               level = 0.95, ...) {
  type <- match.arg(type)
  if (isTRUE(se.fit) || match.arg(interval) != "none")
    return(.ud_predict_ci(object, newdata, type, se.fit, match.arg(interval), level, offset = offset))
  if (!is.null(newdata)) {
    need <- setdiff(all.vars(stats::delete.response(stats::terms(object$participation))), c(names(newdata), object$part_fe))
    if (length(need)) stop("'newdata' is missing participation variable(s): ", paste(need, collapse = ", "), ".")
  }
  p <- .part_at_limit(.ud_part_prob(object, newdata), object, newdata)
  if (type == "participation") return(as.numeric(p))
  ## newdata = NULL: use the stored FULL-length intensity (lambda_full, one per
  ## observation), NOT predict(intensity) which returns only the positives-only
  ## fitted values (length n_pos) and would recycle-misalign against p (length n).
  ey <- if (is.null(newdata)) .cpb_ztmean(object$lambda_full, object$intensity$alpha, .cpb_support(object))
        else predict(object$intensity, newdata = newdata, type = "response", offset = offset)   # E(Y | Y > 0)
  if (type == "intensity") return(as.numeric(ey))
  as.numeric(p * ey)                                  # unnamed, like fitted() and the siblings
}
