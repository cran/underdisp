## ---------------------------------------------------------------------------
## Simulators for the CPB family
## ---------------------------------------------------------------------------

#' Simulate from the continuous parameter binomial
#'
#' Draws counts from the continuous parameter binomial (CPB) or its
#' zero-truncated variant.
#'
#' @param n Number of draws (recycled against `lambda`).
#' @param lambda Mean parameter; a scalar or a length-`n` vector.
#' @param alpha Shape/dispersion parameter in (0, 1).
#' @param truncated If `TRUE`, draw from the zero-truncated CPB.
#' @return An integer vector of counts.
#' @examples
#' set.seed(1)
#' table(rcpb(1000, lambda = 3, alpha = 0.5))
#' @export
rcpb <- function(n, lambda, alpha, truncated = FALSE) {
  if (alpha <= 0 || alpha >= 1) stop("`alpha` must be in (0, 1).")
  lambda <- rep_len(lambda, n)
  out <- integer(n)
  for (i in seq_len(n)) {
    pr <- .cpb_pmf_core(lambda[i], alpha)
    K  <- length(pr) - 1L
    if (truncated) {
      if (K < 1) { out[i] <- 1L; next }
      pr[1] <- 0; pr <- pr / sum(pr)
    }
    out[i] <- as.integer(sample(0:K, 1L, prob = pr))
  }
  out
}

#' Simulate from the hurdle continuous parameter binomial
#'
#' @param n Number of units.
#' @param lambda Intensity mean for participants; scalar or length-`n`.
#' @param alpha CPB shape parameter in (0, 1).
#' @param p Participation probability; scalar or length-`n`.
#' @return An integer vector: 0 for nonparticipants, a zero-truncated CPB draw
#'   otherwise.
#' @examples
#' set.seed(1)
#' y <- rhurdle_cpb(1000, lambda = 3, alpha = 0.5, p = 0.6)
#' mean(y == 0)
#' @export
rhurdle_cpb <- function(n, lambda, alpha, p) {
  d <- stats::rbinom(n, 1, p)
  y <- integer(n)
  idx <- which(d == 1)
  if (length(idx)) y[idx] <- rcpb(length(idx), rep_len(lambda, n)[idx], alpha, truncated = TRUE)
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
#'   participation are uninformative for a fixed-effects logit.
#' @param link Link for the participation model: `"logit"` (default), `"probit"`,
#'   or `"cloglog"`. Positive participation coefficients raise the probability of
#'   a positive count.
#' @param offset Optional offset on the log-mean scale for the intensity (an
#'   exposure): a numeric vector or the name of a column in `data`.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector) for cluster-robust bootstrap inference on the intensity
#'   coefficients; see [cpb()]. Cluster-robust errors for the participation
#'   logit are available directly via
#'   `sandwich::vcovCL(fit$participation, cluster = ...)`.
#' @param se Inference for the intensity coefficients: `"none"` or `"bootstrap"`.
#' @param B Bootstrap replicates when `se = "bootstrap"`.
#' @param ... Passed to [cpb()] or [cpb_fe()].
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
                       cluster = NULL, se = c("none", "bootstrap"), B = 500, ...) {
  se <- match.arg(se); link <- match.arg(link)
  y <- stats::model.response(stats::model.frame(formula, data))
  if (any(y < 0) || any(y != floor(y))) stop("The response must be nonnegative integer counts.")
  d <- as.integer(y > 0)
  int_offset <- if (is.null(offset)) NULL                       # offset applies to the intensity
                else if (is.character(offset) && length(offset) == 1L) offset else offset[y > 0]
  ## carry a cluster identifier as a column so both margins can see it aligned
  if (!is.null(cluster) && !(is.character(cluster) && length(cluster) == 1L)) {
    data[[".cluster"]] <- cluster; cluster <- ".cluster"
  }

  ## participation (hurdle) model: logit of d, optionally with fixed effects
  prhs <- if (is.null(participation)) formula[-2L] else participation
  if (!is.null(part_fe))
    prhs <- stats::reformulate(c(labels(stats::terms(prhs)), paste0("factor(", part_fe, ")")))
  pdata <- data; pdata[[".d"]] <- d
  pfit <- stats::glm(stats::update(prhs, .d ~ .), data = pdata, family = stats::binomial(link = link))

  ## intensity model on the positive counts, with optional unit fixed effects
  pos <- data[y > 0, , drop = FALSE]
  ifit <- if (is.null(fe)) cpb(formula, data = pos, truncated = TRUE, se = se, B = B,
                               cluster = cluster, offset = int_offset, ...)
          else cpb_fe(formula, data = pos, fe = fe, truncated = TRUE, se = se, B = B,
                      cluster = cluster, offset = int_offset, ...)

  ## per-observation intensity lambda (fixed-effects-aware) and a reference
  ## profile, for calibration and the first-difference decomposition
  Xall <- stats::model.matrix(stats::delete.response(stats::terms(formula)), data)
  beta <- ifit$coefficients
  if (is.null(fe)) {
    lambda_full <- as.numeric(exp(Xall %*% beta))
    int_xref <- colMeans(Xall); int_fe_ref <- 0
  } else {
    Xcov <- Xall[, names(beta), drop = FALSE]
    int_fe_ref <- mean(ifit$fe)
    fe_i <- ifit$fe[as.character(data[[fe]])]; fe_i[is.na(fe_i)] <- int_fe_ref
    lambda_full <- as.numeric(exp(fe_i + Xcov %*% beta))
    int_xref <- colMeans(Xcov)
  }

  structure(list(participation = pfit, intensity = ifit, Zpart = stats::model.matrix(pfit),
                 y = as.integer(y), p_full = as.numeric(stats::fitted(pfit)),
                 lambda_full = lambda_full, n = length(y), n_participate = sum(d),
                 int_beta = beta, int_xref = int_xref, int_fe_ref = int_fe_ref,
                 part_beta = stats::coef(pfit), part_xref = colMeans(stats::model.matrix(pfit)),
                 fe = fe, part_fe = part_fe, link = link,
                 cluster = if (is.character(cluster)) cluster else NULL,
                 formula = formula, part_formula = participation, call = match.call()),
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
  pc <- pc[!grepl("^factor\\(", names(pc))]   # hide participation fixed-effect dummies
  lk <- if (is.null(x$link)) "logit" else x$link
  cat(sprintf("\nParticipation (%s link) -- positive coefficients raise P(Y > 0), i.e. participation:\n", lk))
  print(round(pc, 4))
  cat("\nIntensity (zero-truncated CPB) coefficients:\n"); print(round(x$intensity$coefficients, 4))
  cat("Intensity alpha (shape):", round(x$intensity$alpha, 4), "\n")
  invisible(x)
}

#' Predict from a hurdle-CPB fit
#'
#' @param object A `"hurdle_cpb"` object.
#' @param newdata Data frame of covariate profiles (must contain both the
#'   intensity and participation covariates).
#' @param type `"response"` for the marginal expected count E(Y) (zeros
#'   included), `"participation"` for P(Y > 0), or `"intensity"` for
#'   E(Y | Y > 0).
#' @param ... Unused.
#' @return A numeric vector.
#' @method predict hurdle_cpb
#' @export
predict.hurdle_cpb <- function(object, newdata = NULL,
                               type = c("response", "participation", "intensity"), ...) {
  type <- match.arg(type)
  p <- stats::predict(object$participation, newdata = newdata, type = "response")
  if (type == "participation") return(p)
  ## newdata = NULL: use the stored FULL-length intensity (lambda_full, one per
  ## observation), NOT predict(intensity) which returns only the positives-only
  ## fitted values (length n_pos) and would recycle-misalign against p (length n).
  ey <- if (is.null(newdata)) object$lambda_full
        else predict(object$intensity, newdata = newdata, type = "response")
  if (type == "intensity") return(ey)
  p * ey
}
