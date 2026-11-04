## ---------------------------------------------------------------------------
## The participation equation of the hurdles (hurdle_cpb, hurdle_gec,
## hurdle_count): a binary glm of P(Y > 0) on every row.
##
## Its covariance is the glm's model-based one, or the sandwich with the meat
## of the count families' (.count_vcov): per-row scores (hurdle_count with
## se = "robust"), or scores summed within clusters with the G/(G - 1) factor
## (a fit given `cluster`), frequency weights entering as duplicated rows. The
## hurdles store it as `part_vcov`, with `part_vcov_type`, and every table,
## interval and quantity of interest reads it from there. A clustered covariance
## asked for without a cluster is NA, as count_reg()'s is (hurdle_count's
## intensity fit says so).
##
## Its separation is checked by the zero-inflated fits' rule (.zi_separation;
## see ?zi_cpb, Separation), which on a binary equation is the profile of the
## glm's own likelihood: a covariate coefficient whose profile 100 standard
## deviations of its covariate further from zero is less than
## qchisq(0.95, 1) / 2 below its value at the estimate is flagged. The intercept
## and the part_fe dummies are not checked (a unit that always or never
## participates has its dummy at infinity, as in any fixed-effects logit).
## ---------------------------------------------------------------------------

## The levels of a two-part fit's fixed-effects column whose rows have no positive count (none when no row has one:
## the fits then warn that every count is zero). Such a unit adds nothing to the likelihood: its dummy in the
## participation or zero equation, and its count effect, sit at their limits, where every one of its rows has
## probability one. As fixest leaves such units out of a fixed-effects fit, the hurdle and zero-inflated fits do not
## count those parameters in their degrees of freedom, whichever equation carries them, so that the two kinds of
## model are charged alike on the same data (CPB fits on the floor support keep the count of versions before 0.2.0).
.ud_zero_levels <- function(y, unit) {
  if (is.null(unit)) return(character(0))
  pos <- tapply(y > 0, as.character(unit), any)
  if (!any(pos)) character(0) else names(pos)[!pos]
}

## The levels of a hurdle's participation fixed-effects column whose rows are all positive. Such a unit's dummy is at
## plus infinity, as in any fixed-effects logit; unlike the dummy of a unit with no positive count it stays counted in
## the degrees of freedom, as the zero-inflated fits count their zero-equation dummy for the unit, which sits at its
## own limit (pi = 0) there. The fit records the units so that new data can be predicted at the limit (.part_at_limit).
.ud_one_levels <- function(y, unit) {
  if (is.null(unit)) return(character(0))
  allpos <- tapply(y > 0, as.character(unit), all)
  names(allpos)[allpos]
}

## The participation probability of a hurdle in new data, at the limit for the units .ud_zero_levels() and
## .ud_one_levels() found. The dummy of a unit with no positive count is at minus infinity, and that of a unit whose
## every count is positive at plus infinity, but the glm stops at a finite value for either, and that value plus the
## covariates can give the unit's new rows any probability; the limit gives them zero or one, whatever the
## covariates. The zero-inflated fits' count effect at its limit gives the first kind of unit the same zero.
.part_at_limit <- function(p, object, newdata) {
  zl <- object$part_zero_units; ol <- object$part_one_units; pf <- object$part_fe
  if (is.null(newdata) || !(length(zl) || length(ol)) || length(pf) != 1L || is.null(newdata[[pf]])) return(p)
  u <- as.character(newdata[[pf]])
  p[u %in% zl] <- 0
  p[u %in% ol] <- 1
  p
}

## The mixtures' analogue for their zero equation (zi_cpb() on King's support, zi_gec(), zi_count()): a zero_fe
## unit with no positive count has its dummy at plus infinity, so its rows in new data get the structural-zero
## probability one, and a unit whose every count is positive has it at minus infinity, so they get zero. With the
## count effect of the first kind of unit at its own limit as well, P(Y = 0) is one either way.
.zero_at_limit <- function(p, object, newdata) {
  zl <- object$zero_fe_zero_units; ol <- object$zero_fe_one_units; zf <- .ud_zero_fe(object)
  if (is.null(newdata) || !(length(zl) || length(ol)) || length(zf) != 1L || is.null(newdata[[zf]])) return(p)
  u <- as.character(newdata[[zf]])
  p[u %in% zl] <- 1
  p[u %in% ol] <- 0
  p
}
.ud_zero_fe <- function(object) if (!is.null(object$zero_fe)) object$zero_fe else object$refit_args$zero_fe

## The rows of new data whose count-equation fixed-effects unit has no positive count in the fitted data (its
## effect at minus infinity: every count of the unit is zero), for the mixtures whose count component carries the
## unit effect: zi_cpb(fe =) on King's support (`zero_units`) and zi_count() with one `fe` column (`fe_zero_units`).
.count_at_limit <- function(object, newdata) {
  if (is.null(newdata)) return(logical(0))
  zu <- if (inherits(object, "zi_cpb")) object$zero_units else object$fe_zero_units
  fv <- object$fe
  if (!length(zu) || length(fv) != 1L || is.null(newdata[[fv]])) return(rep(FALSE, nrow(newdata)))
  as.character(newdata[[fv]]) %in% zu
}

## --- new data: units the fit has not seen ---------------------------------------------------------------------
## An equation whose fixed effects enter as factor dummies (a hurdle's participation glm, a mixture's zero equation,
## the count equations of count_reg(), hurdle_count(), zi_count() and hurdle_gec()) cannot evaluate a dummy for a
## unit the fit did not see. Such a row, or every row when `newdata` lacks the fixed-effects column, enters at the
## reference profile `xref` of that factor: the average over the units whose effect the data identify
## (.fe_dummy_ref()), the reference unit first_difference() uses. `fe` names the fixed-effects columns; with
## `all_factors` an unseen level of any factor enters at its reference profile in the same way (held-out scoring,
## where a held-out block can hold levels its training folds lack). An object without `xref` (saved before 0.2.0)
## uses the mean of the factor's estimated effects, the reference level at zero, as held-out scoring did then.

## the variable behind a model-frame factor name (`factor(iso3)` -> `iso3`)
.ud_xlev_var <- function(fv) if (grepl("^factor\\(.*\\)$", fv)) sub("^factor\\((.*)\\)$", "\\1", fv) else fv

## whether newdata holds a level (or lacks a fixed-effects column) that .ud_dummy_eta() has to place
.ud_has_unseen <- function(xlev, newdata, fe = NULL, all_factors = FALSE) {
  if (is.null(newdata)) return(FALSE)
  for (fv in names(xlev)) {
    v <- .ud_xlev_var(fv)
    if (!(all_factors || v %in% fe)) next
    if (!v %in% names(newdata)) { if (v %in% fe) return(TRUE); next }
    vals <- as.character(newdata[[v]])
    if (any(!is.na(vals) & !(vals %in% xlev[[fv]]))) return(TRUE)
  }
  FALSE
}

.ud_dummy_eta <- function(terms_obj, xlev, coefs, newdata, contrasts = NULL, xref = NULL, fe = NULL,
                          all_factors = FALSE) {
  coefs[is.na(coefs)] <- 0
  extra <- numeric(nrow(newdata))
  for (fv in names(xlev)) {
    v <- .ud_xlev_var(fv)
    if (!(all_factors || v %in% fe)) next
    if (!v %in% names(newdata)) {
      if (!v %in% fe) next                            # the design matrix then reports the missing variable
      unseen <- rep(TRUE, nrow(newdata)); newdata[[v]] <- rep(xlev[[fv]][1L], nrow(newdata))
    } else {
      vals <- as.character(newdata[[v]])
      unseen <- !is.na(vals) & !(vals %in% xlev[[fv]])
      if (!any(unseen)) next
      vals[unseen] <- xlev[[fv]][1L]; newdata[[v]] <- vals   # the reference level, whose dummies are zero
    }
    cols <- intersect(paste0(fv, xlev[[fv]]), names(coefs))
    ref <- if (length(cols) && !is.null(xref) && all(cols %in% names(xref))) sum(xref[cols] * coefs[cols])
           else if (length(cols)) mean(c(0, unname(coefs[cols]))) else 0
    extra[unseen] <- extra[unseen] + ref
  }
  X <- .ud_newdata_matrix(terms_obj, newdata, xlev, contrasts, names(coefs))
  as.numeric(X %*% coefs) + extra
}

## A hurdle's participation probability in new data before its limits (.part_at_limit()): the glm's own prediction,
## or, where newdata holds a unit the fit has not seen (or no fixed-effects column), .ud_dummy_eta() at the
## participation equation's reference unit
.ud_part_prob <- function(object, newdata, all_factors = FALSE) {
  pg <- object$participation
  if (is.null(newdata) || !.ud_has_unseen(pg$xlevels, newdata, object$part_fe, all_factors))
    return(.ud_glm_predict(pg, newdata))
  eta <- .ud_dummy_eta(stats::delete.response(stats::terms(pg)), pg$xlevels, stats::coef(pg), newdata,
                       pg$contrasts, object$part_xref, object$part_fe, all_factors)
  as.numeric(stats::family(pg)$linkinv(eta))
}

## --- several fixed-effects columns -----------------------------------------------------------------------------
## the levels, over every column of `fe`, with no positive count (their dummies at minus infinity)
.ud_fe_zero_levels <- function(y, data, fe) {
  z <- unlist(lapply(fe, function(f) .ud_zero_levels(y, data[[f]])))
  if (is.null(z)) character(0) else z
}
## the levels whose every count is positive, for a single fixed-effects column (predicted at the limit, which
## .part_at_limit() and .zero_at_limit() apply to one column)
.ud_fe_one_levels <- function(y, data, fe) if (length(fe) == 1L) .ud_one_levels(y, data[[fe]]) else character(0)
## the groups of the separation check's count of extreme probabilities: the fixed-effects units, or their
## combinations when there are several columns
.ud_fe_group <- function(data, fe)
  if (is.null(fe)) NULL else if (length(fe) == 1L) data[[fe]] else interaction(data[fe], drop = TRUE)

## a hurdle's degrees of freedom: the stored count, or for an object without one the count of its two parts
.ud_twopart_df <- function(m) {
  if (!is.null(m$df)) return(m$df)
  idf <- if (!is.null(m$intensity$df)) m$intensity$df else length(m$intensity$coefficients) + 1L
  length(stats::coef(m$participation)) + idf
}

## the print methods' note on the units .ud_zero_levels() found (`one`, `many`: the parameters, singular and plural)
.ud_zero_df_note <- function(units, one, many) {
  k <- length(units)
  if (!k) return(invisible(NULL))
  cat(if (k == 1L) sprintf("Note: 1 unit has no positive count. Its %s is at the limit", one)
      else sprintf("Note: %d units have no positive count. Their %s are at the limit", k, many),
      " and not counted in the degrees of freedom.\n", sep = "")
}

.part_fit_vcov <- function(pfit, type = "model", cluster = NULL) {
  V <- tryCatch(stats::vcov(pfit), error = function(e) NULL)
  if (is.null(V) || type == "model") return(V)
  if (type == "cluster" && is.null(cluster)) { V[] <- NA_real_; return(V) }
  ok <- !is.na(stats::coef(pfit)); X <- stats::model.matrix(pfit)[, ok, drop = FALSE]
  mu <- pfit$fitted.values; w <- pfit$prior.weights
  S <- X * ((pfit$y - mu) * pfit$family$mu.eta(pfit$linear.predictors) / pmax(mu * (1 - mu), 1e-300))
  meat <- if (type == "cluster") {
    G <- length(unique(cluster))
    (G / (G - 1)) * crossprod(as.matrix(rowsum(S * w, group = cluster)))
  } else crossprod(S * sqrt(w))
  bread <- V[ok, ok, drop = FALSE]
  V[ok, ok] <- bread %*% meat %*% bread
  V
}

## the participation equation's covariance of a hurdle (an object from before 0.2.0: the glm's)
.part_vcov <- function(object) {
  if (!is.null(object$part_vcov)) return(object$part_vcov)
  tryCatch(stats::vcov(object$participation), error = function(e) NULL)
}

## the separation check of the participation glm: its successes enter .zi_separation() as the "zeros" of a
## mixture whose count component has no zero mass, which leaves exactly the binary likelihood
.part_separation <- function(pfit, part_fe = NULL, group = NULL) {
  g <- stats::coef(pfit); ok <- !is.na(g)
  Z <- stats::model.matrix(pfit)[, ok, drop = FALSE]
  sep <- tryCatch(.zi_separation(Z, g[ok], ifelse(pfit$y > 0.5, 0L, 1L), rep(0, length(pfit$y)),
                                 w = pfit$prior.weights, group = group, fe_cols = .zi_fe_cols(Z, part_fe),
                                 link = pfit$family$link),
                  error = function(e) NULL)
  if (is.null(sep)) return(NULL)
  sep[c("diverging", "profile", "n_extreme", "groups_extreme")]
}

.part_warn_separation <- function(sep, fun) {
  if (is.null(sep) || !length(sep$diverging)) return(invisible(NULL))
  pl <- length(sep$diverging) > 1L
  warning(fun, ": the participation equation is separated or nearly so. The profile likelihood of the coefficient",
          if (pl) "s", " on ", paste(sep$diverging, collapse = ", "), " falls by less than qchisq(0.95, 1) / 2 over ",
          "100 standard deviations of ", if (pl) "their covariates" else "its covariate", " further from zero, so ",
          "the data put no useful bound on ", if (pl) "them" else "it", ", and ", if (pl) "their" else "its",
          " standard error, odds ratio and effect on the participation probability are not meaningful (",
          sep$n_extreme, " observations have a participation probability numerically 0 or 1).", call. = FALSE)
}

## The participation equation of a fitted hurdle: its covariance (model-based, robust, or clustered on `cl`) and
## separation check, as fields of the fit
.part_fields <- function(pfit, part_fe, data, type, cl, fun) {
  sep <- .part_separation(pfit, part_fe, .ud_fe_group(data, part_fe))
  .part_warn_separation(sep, fun)
  list(part_vcov = .part_fit_vcov(pfit, type, cl), part_vcov_type = type,
       part_n_clusters = if (type == "cluster") length(unique(cl)) else NA_integer_, part_separation = sep)
}

## the dummy columns of fixed-effects variables `fe` among the names `nm`
.is_fe_col <- function(nm, fe) {
  out <- rep(FALSE, length(nm))
  for (v in fe) out <- out | startsWith(nm, paste0("factor(", v, ")"))
  out
}

## The participation table of a hurdle's summary: the participation covariance's standard errors, the part_fe
## dummies left out, and the separated terms without a standard error
.part_block <- function(object) {
  b <- object$part_beta; V <- .part_vcov(object)
  se <- if (is.null(V)) NULL else sqrt(pmax(diag(V), 0))[names(b)]
  keep <- !.is_fe_col(names(b), object$part_fe) & !is.na(b)
  .sum_mask(.sum_block(b[keep], if (is.null(se)) NULL else se[keep]), object$part_separation$diverging)
}
.part_note <- function(object) {
  ty <- object$part_vcov_type
  c(if (identical(ty, "cluster") && isTRUE(object$part_n_clusters > 0))
      sprintf("Participation standard errors are cluster-robust (%d clusters).", object$part_n_clusters)
    else if (identical(ty, "robust")) "Participation standard errors are heteroskedasticity-robust.",
    .sep_note(object$part_separation, "participation"))
}

## a summary table's separated rows lose their standard error, z and p
.sum_mask <- function(bl, terms) {
  if (is.null(bl) || !length(terms) || ncol(bl) < 4L) return(bl)
  bl[rownames(bl) %in% terms, 2:4] <- NA
  bl
}
.sep_note <- function(sep, equation) {
  if (is.null(sep) || !length(sep$diverging)) return(NULL)
  paste0("Note: the data put no useful bound on ", paste(sep$diverging, collapse = ", "), " in the ", equation,
         " equation (separation); ", if (length(sep$diverging) > 1L) "their standard errors are" else
           "its standard error is", " not shown.")
}

## the separated terms of a fit (zero-inflation or participation equation), for irr(), zi_test() and compare_models()
.ud_separated <- function(object) c(object$separation$diverging, object$part_separation$diverging)
