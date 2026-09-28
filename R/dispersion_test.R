## ---------------------------------------------------------------------------
## Regression-adjusted tests of equidispersion.
##
## Two tests, one interface. (1) The likelihood-ratio test of a fitted family's
## dispersion parameter against its Poisson value, on the fitted model's own
## design: the Poisson is the equidispersed member of every family in the
## package (CPB alpha = 1, GEC delta = 1, COM-Poisson nu = 1, generalized
## Poisson lambda = 0, gamma-count alpha = 1, double Poisson theta = 1, negative
## binomial 1/theta = 0). (2) The auxiliary-regression test of Cameron and
## Trivedi (1990) for a Poisson fit: regress ((y - mu)^2 - y) / mu on mu through
## the origin; the t statistic of the slope tests Var = mu against
## Var = mu + a * mu^2, with a < 0 underdispersion.
##
## Reference distribution of the likelihood-ratio statistic, one rule.
## Asymptotically the statistic is chi-square(1) where the Poisson value is
## interior and the 0.5 chi^2_0 + 0.5 chi^2_1 mixture where it is on the
## boundary (CPB, negative binomial). The CPB qualifies although its support
## moves with its parameters: with theta = 1 - alpha its ceiling lambda / theta
## diverges as theta -> 0, the probability beyond the ceiling vanishes faster
## than any power of theta, and the log-likelihood expands one-sidedly as the
## Poisson's plus theta * [k - lambda/2 - k(k-1)/(2 lambda)] + O(theta^2) per
## observation, whose score is the classical dispersion score
## -[(y - lambda)^2 - y] / (2 lambda), with information 1/2 per observation and
## orthogonal to the regression coefficients (Self and Liang 1987; Andrews 2001);
## the zero-truncated model reaches the same mixture through its efficient score. The asymptotic
## distribution ignores that the null's mean parameters are estimated: to first
## order the signed root moves toward underdispersion by p / sqrt(2 n) standard
## deviations for p mean parameters (unit intercepts included) on n observations
## (the bias Dean and Lawless 1989 remove from the score test). The p-value is a
## parametric bootstrap under the Poisson fitted to the same design whenever the
## Poisson value is on the boundary (the one-sided statistic is skewed toward
## rejection in finite samples: the package's size simulations found the
## asymptotic CPB test over-rejecting with 60 to 200 observations) or the fit has
## unit fixed effects; where the Poisson value is interior the asymptotic
## distribution is used unless that shift pushes a 5% test's first-order size
## above 6%. The
## bootstrap simulates at the Poisson value itself, so it is not the bootstrap
## from a boundary estimate that fails.
## ---------------------------------------------------------------------------

#' Regression-adjusted test of equidispersion
#'
#' Tests whether a fitted count model's dispersion departs from the Poisson,
#' conditional on the covariates. For a fit from a family with a dispersion
#' parameter (`method = "lr"`, the default), the test is the likelihood-ratio
#' test of that parameter against its Poisson value, with the Poisson refit on
#' the same design (fixed effects included); its p-value comes from the
#' asymptotic distribution or from a parametric bootstrap under the fitted
#' Poisson, by the rule in Details. For a Poisson fit, `method = "auxiliary"`
#' runs the Cameron--Trivedi (1990) auxiliary regression, which needs no
#' alternative family.
#'
#' @details
#' The Poisson is nested in every family at one value of the dispersion
#' parameter: `alpha = 1` (CPB), `delta = 1` (GEC), `nu = 1` (both COM-Poisson
#' parameterizations), `lambda = 0` (generalized Poisson), `alpha = 1`
#' (gamma-count), `theta = 1` (double Poisson), and `1/theta = 0` (negative
#' binomial). When that value is interior to the parameter space the
#' likelihood-ratio statistic is asymptotically \eqn{\chi^2_1} under the null; a
#' directional alternative (`"under"` or `"over"`) uses the signed root
#' \eqn{r = \pm\sqrt{LR}}, positive when the estimate lies on the underdispersed
#' side of the Poisson value, which is asymptotically standard normal. When the
#' Poisson value is on the boundary (the CPB, whose `alpha` lives in (0, 1); the
#' negative binomial, whose overdispersion parameter is non-negative) the
#' asymptotic null distribution is the
#' \eqn{\tfrac12\chi^2_0 + \tfrac12\chi^2_1} mixture (Self and Liang 1987;
#' Andrews 2001) and the test is one-sided by construction (`"under"` for the
#' CPB, `"over"` for the negative binomial). The CPB qualifies although its
#' support moves with its parameters: near `alpha = 1` its ceiling
#' `lambda / (1 - alpha)` diverges, the log-likelihood is regular on that side,
#' and its score there is the classical dispersion score
#' \eqn{-[(y - \lambda)^2 - y] / (2\lambda)}; the zero-truncated model reaches
#' the same mixture through its efficient score. Because the Poisson is the CPB's
#' limit as `alpha` approaches 1, the statistic cannot be negative in principle;
#' when the `max.support` guard stops `alpha` short of that limit and the fitted
#' log-likelihood lands just below the Poisson's, the statistic takes its
#' boundary value 0. Any other negative statistic means the optimizer fell
#' short, and no p-value is reported.
#'
#' **Calibration.** Where the Poisson value is on the boundary of the family's
#' parameter space (the CPB, the negative binomial), or the fit has unit fixed
#' effects, the p-value is by default a parametric bootstrap with 199
#' replicates. The one-sided boundary statistic is skewed toward rejection in
#' finite samples, and unit intercepts bias the dispersion estimate toward
#' underdispersion by an amount of order 1/T; in both cases the asymptotic
#' distribution over-rejects. Where the Poisson value is interior, the
#' asymptotic distribution is used unless the estimated mean parameters shift
#' it too far: to first order they move the signed root toward
#' underdispersion by \eqn{p/\sqrt{2n}} standard deviations for \eqn{p} mean
#' parameters on \eqn{n} observations (the weight total), the bias Dean and
#' Lawless (1989) remove from the score test, and the bootstrap takes over when
#' that shift pushes the first-order size of a 5\% test above 6\%. The bootstrap
#' simulates responses from the Poisson (zero-truncated Poisson) fitted to the
#' same design, offset, weights, and unit effects, refits both models to each,
#' and reports the share of simulated statistics at least as extreme as the
#' observed one, counting the observed one, among the replicates whose two fits
#' both succeeded. Simulating at the Poisson value itself keeps the bootstrap
#' valid at the boundary. `B = 0` requests the asymptotic distribution (with
#' unit fixed effects it returns the statistic without a p-value), and a
#' positive `B` sets the number of bootstrap replicates. Every replicate refits
#' both models, so the bootstrap takes about `B` times as long as the original
#' fit; `cores` spreads the replicates over worker processes.
#'
#' @param object A fitted `cpb`, `cpb_fe`, `gec`, `gec_fe`, or `count_reg`
#'   model (for `method = "auxiliary"`, a `count_reg` fit with
#'   `family = "poisson"`).
#' @param alternative `"two.sided"` (default where the null is interior),
#'   `"under"` (underdispersion), or `"over"` (overdispersion).
#' @param method `"lr"` (likelihood ratio against the Poisson refit; default)
#'   or `"auxiliary"` (the Cameron--Trivedi auxiliary regression for a
#'   Poisson fit).
#' @param B Parametric-bootstrap replicates for `method = "lr"`: `NULL` (the
#'   default) applies the calibration rule in Details (199 replicates when the
#'   bootstrap is needed), `0` requests the asymptotic distribution, and a
#'   positive whole number forces the bootstrap.
#' @param cores Worker processes for the bootstrap replicates.
#' @return An object of class `c("dispersion_test", "htest")` with the
#'   statistic, the p-value, the estimated dispersion parameter and its Poisson
#'   value, the two log-likelihoods, and the alternative; for `method = "lr"`
#'   also `calibration` (`"asymptotic"` or `"parametric bootstrap"`), `B` and
#'   `B_ok` (replicates requested and completed), `boot` (the simulated
#'   statistics), `first_order_size` (the asymptotic 5\% test's first-order
#'   size), `asymptotic_ok` (whether the calibration rule admits the asymptotic
#'   distribution), and `note` where a statistic or p-value needs one.
#' @references Andrews, D. W. K. (2001). Testing when a parameter is on the
#'   boundary of the maintained hypothesis. \emph{Econometrica}, 69(3), 683-734.
#'   Cameron, A. C. and
#'   Trivedi, P. K. (1990). Regression-based tests for overdispersion in the
#'   Poisson model. \emph{Journal of Econometrics}, 46(3), 347-364. Dean, C. and
#'   Lawless, J. F. (1989). Tests for detecting overdispersion in Poisson
#'   regression models. \emph{Journal of the American Statistical Association},
#'   84(406), 467-472. Self, S. G. and Liang, K.-Y. (1987). Asymptotic
#'   properties of maximum likelihood estimators and likelihood ratio tests
#'   under nonstandard conditions. \emph{Journal of the American Statistical
#'   Association}, 82(398), 605-610.
#' @examples
#' set.seed(2); n <- 400; x <- rnorm(n)
#' d <- data.frame(y = rgammacount(n, exp(1 + 0.4 * x), alpha = 2), x = x)
#' dispersion_test(count_reg(y ~ x, d, family = "gammacount"), alternative = "under")
#' dispersion_test(count_reg(y ~ x, d, family = "poisson"), method = "auxiliary")
#' @seealso [count_reg()], [dispersion_profile()], [zi_test()]
#' @export
dispersion_test <- function(object, alternative = c("two.sided", "under", "over"),
                            method = c("lr", "auxiliary"), B = NULL, cores = 1L) {
  method <- match.arg(method); alt_missing <- missing(alternative)
  alternative <- match.arg(alternative)
  if (method == "auxiliary") return(.disp_auxiliary(object, alternative))
  if (!is.null(B) && (length(B) != 1L || !is.numeric(B) || is.na(B) || B < 0 || B != round(B)))
    stop("B must be NULL or a non-negative whole number.")
  .disp_lr(object, alternative, alt_missing, B, cores)
}

## the asymptotic distribution is kept while the first-order size of a 5% test
## stays within one percentage point of nominal
.disp_size_tol <- 0.06

## first-order size of the asymptotic 5% test when the signed root sits
## p / sqrt(2 n) standard deviations toward underdispersion
.disp_first_order_size <- function(p, n, alternative) {
  s <- p / sqrt(2 * n)
  switch(alternative,
    under = stats::pnorm(stats::qnorm(0.95) - s, lower.tail = FALSE),
    over = stats::pnorm(-stats::qnorm(0.95) - s),
    two.sided = stats::pnorm(stats::qnorm(0.975) - s, lower.tail = FALSE) + stats::pnorm(-stats::qnorm(0.975) - s))
}

## family label, dispersion parameter, its Poisson value, boundary status, and
## the side of the Poisson value that is underdispersion
.disp_family <- function(object) {
  if (inherits(object, c("cpb", "cpb_fe")))
    return(list(label = "CPB", par = "alpha", est = object$alpha, null = 1, boundary = TRUE,
                forced = "under", under = function(est) est < 1))
  if (inherits(object, c("gec", "gec_fe")))
    return(list(label = "GEC", par = "delta", est = object$delta, null = 1, boundary = FALSE,
                forced = NULL, under = function(est) est < 1))
  if (inherits(object, "count_reg")) {
    fam <- .count_fam(object$family)
    if (!fam$nshape) stop("The fitted family (Poisson) has no dispersion parameter; use method = \"auxiliary\" ",
                          "or fit a family with a dispersion parameter.")
    if (object$family == "negbin")
      return(list(label = fam$label, par = "1/theta", est = 1 / object$theta, null = 0, boundary = TRUE,
                  forced = "over", under = function(est) FALSE))
    nul <- fam$null_shape
    return(list(label = fam$label, par = sub("^(log|atanh)\\((.*)\\)$", "\\2", fam$shape_name),
                est = object$theta, null = nul, boundary = FALSE, forced = NULL,
                under = switch(object$family,
                  compois = , mpcmp = , gammacount = , doublepois = function(est) est > nul,
                  genpois = function(est) est < nul)))
  }
  stop("dispersion_test() takes a cpb, cpb_fe, gec, gec_fe, or count_reg fit.")
}

## unit fixed effects, the number of mean parameters (unit intercepts
## included), and the observation count (the weight total)
.disp_design <- function(object) {
  conc <- inherits(object, c("cpb_fe", "gec_fe"))
  list(fe = conc || (inherits(object, "count_reg") && !is.null(object$fe)),
       p = ncol(object$X) + if (conc) length(unique(object$unit)) else 0L,
       n = if (is.null(object$weights)) length(object$Y) else sum(object$weights))
}

.disp_offset <- function(object)
  if (is.null(object$offset)) rep(0, length(object$Y)) else rep_len(object$offset, length(object$Y))

## the likelihood-ratio statistic against the Poisson nest: a negative value
## that the support guard produced (the fit stopped short of the alpha -> 1
## limit) is the boundary value 0; any other negative value is an optimizer
## shortfall and carries no p-value
.disp_negative_lr <- function(LR, binding, label) {
  if (is.na(LR) || LR >= -1e-8) return(list(LR = if (is.na(LR)) NA_real_ else max(LR, 0), shortfall = FALSE, note = NULL))
  if (binding)
    return(list(LR = 0, shortfall = FALSE,
                note = paste0("the max.support guard stopped alpha short of its Poisson limit, where the fitted ",
                              "log-likelihood falls just below the Poisson's, so the statistic takes its boundary value 0")))
  list(LR = LR, shortfall = TRUE,
       note = paste0("the fitted ", label, " log-likelihood lies below its Poisson nest (LR < 0): the optimizer ",
                     "fell short, and no p-value is reported"))
}

## the CPB statistic with its asymptotic boundary-mixture p-value (summary.cpb)
.disp_boundary_lr <- function(LR, binding) {
  v <- .disp_negative_lr(LR, binding, "CPB")
  p <- if (is.na(v$LR) || v$shortfall) NA_real_ else 0.5 * stats::pchisq(v$LR, 1, lower.tail = FALSE)
  list(LR = v$LR, p = p, note = v$note, shortfall = v$shortfall)
}

## Poisson (zero-truncated Poisson) fitted to the fit's own design, offset,
## weights, and unit effects, for the response y: log-likelihood and rates
.disp_null_fit <- function(object, y = object$Y) {
  if (inherits(object, c("cpb_fe", "gec_fe"))) return(.disp_fe_null(object, y))
  ft <- .count_fit(object$X, y, .count_fam("poisson"), isTRUE(object$truncated),
                   offset = .disp_offset(object), w = object$weights)
  list(loglik = ft$loglik, rate = ft$mu)
}

## Poisson null for the concentrated fixed-effects fits, with the unit effects
## profiled out: in closed form for the Poisson (a_u = log(sum_t w y) - log(sum_t
## w exp(eta))) and by a vectorized Newton iteration on all units at once for the
## zero-truncated Poisson (score sum_t w (y - m), m = lambda / (1 - e^-lambda))
.disp_fe_null <- function(object, y = object$Y) {
  X <- object$X; Y <- y; u <- object$unit
  if (is.null(X) || is.null(u)) stop("The fixed-effects fit does not carry its design; refit with the current package version.")
  w <- .ud_w1(object$weights, length(Y))
  off <- .disp_offset(object)
  trunc <- isTRUE(object$truncated)
  s <- .ud_colscale(X, w); Xs <- sweep(X, 2, s, "/")
  sy <- as.numeric(rowsum(w * Y, u)[, 1])
  conc <- function(b) {
    eta <- off + as.numeric(Xs %*% b)
    a <- log(pmax(sy, 1e-300)) - log(as.numeric(rowsum(w * exp(eta), u)[, 1]))
    a[sy <= 0] <- -30
    if (!trunc) { lam <- exp(a[u] + eta); return(list(ll = sum(w * stats::dpois(Y, lam, log = TRUE)), lam = lam)) }
    ## a unit whose counts are all one has its zero-truncated Poisson supremum
    ## at lambda -> 0 (log-likelihood 0): it is held at the floor and skipped
    ones <- as.numeric(rowsum(w * (Y > 1), u)[, 1]) <= 0
    a[ones] <- -8 - max(eta)
    for (it in seq_len(60L)) {
      lam <- exp(a[u] + eta); em <- exp(-lam); m <- lam / (1 - em)
      g <- as.numeric(rowsum(w * (Y - m), u)[, 1])
      h <- as.numeric(rowsum(w * ((1 - em * (1 + lam)) / (1 - em)^2) * lam, u)[, 1])
      step <- pmin(pmax(g / pmax(h, 1e-10), -2), 2)
      step[ones | !is.finite(step)] <- 0
      a <- pmax(a + step, -8 - max(eta))
      if (max(abs(step)) < 1e-9) break
    }
    lam <- exp(a[u] + eta)
    list(ll = sum(w * (stats::dpois(Y, lam, log = TRUE) - log1p(-exp(-lam)))), lam = lam)
  }
  b0 <- tryCatch({ v <- stats::glm.fit(cbind(1, Xs), Y, weights = w, offset = off, family = stats::poisson())$coefficients[-1]
                   v[!is.finite(v)] <- 0; v }, error = function(e) rep(0, ncol(X)))
  o <- stats::optim(b0, function(b) { v <- conc(b)$ll; if (is.finite(v)) -v else 1e10 },
                    method = "BFGS", control = list(maxit = 1000, reltol = 1e-12))
  list(loglik = -o$value, rate = conc(o$par)$lam)
}

## Frequency weights in the parametric bootstrap. A weight of w stands for w
## independent observations, so a replicate carries w independent draws for that
## row: the rows are expanded and the refits run unweighted. (One draw per row,
## with the weight passed to the refit, would make the w copies identical and
## inflate the replicate statistics by about sum(w^2)/sum(w).) Returns NULL when
## the weights are not whole numbers, where no such expansion exists.
.disp_expand <- function(object) {
  w <- object$weights
  if (is.null(w)) return(list(object = object, index = seq_along(object$Y)))
  if (any(abs(w - round(w)) > 1e-8)) return(NULL)
  idx <- rep(seq_along(w), times = as.integer(round(w)))
  ob <- object
  ob$offset <- .disp_offset(object)[idx]
  ob$X <- object$X[idx, , drop = FALSE]; ob$Y <- object$Y[idx]
  if (!is.null(object$unit)) ob$unit <- object$unit[idx]
  ob["weights"] <- list(NULL)
  list(object = ob, index = idx)
}

## Poisson (zero-truncated Poisson) draws at the given rates
.disp_simulate <- function(rate, truncated) {
  if (!truncated) return(stats::rpois(length(rate), rate))
  y <- stats::qpois(stats::runif(length(rate), exp(-rate), 1), rate)
  y[!is.finite(y) | y < 1] <- 1
  as.integer(y)
}

## a data frame and formula that rebuild the design matrix X column for column
## (the intercept, when X has one, from the formula), with an optional unit
## column, so a new response can be refit exactly as the model was
.disp_frame <- function(X, y, unit = NULL) {
  ic <- match("(Intercept)", colnames(X))
  Z <- if (is.na(ic)) X else X[, -ic, drop = FALSE]
  d <- data.frame(.y = y)
  nm <- if (ncol(Z)) paste0(".v", seq_len(ncol(Z))) else character(0)
  for (j in seq_along(nm)) d[[nm[j]]] <- as.numeric(Z[, j])
  if (!is.null(unit)) d$.unit <- unit
  rhs <- if (length(nm)) paste(nm, collapse = " + ") else "1"
  if (is.null(unit) && is.na(ic)) rhs <- paste("0 +", rhs)
  list(data = d, formula = stats::as.formula(paste(".y ~", rhs)))
}

## the fitted family refit to the response y on the fit's own design, offset,
## weights, unit effects, and support guard: log-likelihood and dispersion estimate
.disp_alt_fit <- function(object, y) {
  off <- .disp_offset(object); w <- object$weights
  if (inherits(object, "count_reg")) {
    ft <- .count_fit(object$X, y, .count_fam(object$family), isTRUE(object$truncated), offset = off, w = w)
    return(c(loglik = ft$loglik, est = if (object$family == "negbin") 1 / ft$theta else ft$theta))
  }
  conc <- inherits(object, c("cpb_fe", "gec_fe"))
  fr <- .disp_frame(object$X, y, if (conc) object$unit)
  if (inherits(object, "cpb_fe")) {
    f <- cpb_fe(fr$formula, fr$data, fe = ".unit", truncated = isTRUE(object$truncated), se = "none",
                offset = off, weights = w, max.support = object$max.support)
    return(c(loglik = f$loglik, est = f$alpha))
  }
  if (inherits(object, "gec_fe")) {
    f <- gec_fe(fr$formula, fr$data, fe = ".unit", se = "none", offset = off, weights = w,
                max.support = object$max.support)
    return(c(loglik = f$loglik, est = f$delta))
  }
  if (inherits(object, "cpb")) {
    f <- cpb(fr$formula, fr$data, truncated = isTRUE(object$truncated), se = "none",
             offset = off, weights = w, max.support = object$max.support)
    return(c(loglik = f$loglik, est = f$alpha))
  }
  f <- gec(fr$formula, fr$data, truncated = isTRUE(object$truncated), se = "none",
           offset = off, weights = w, max.support = object$max.support)
  c(loglik = f$loglik, est = f$delta)
}

## the statistic compared across replicates: the likelihood ratio (at least 0)
## for boundary and two-sided tests, the signed root for a directional test
.disp_signed <- function(LR, est, fam, alternative) {
  if (fam$boundary || alternative == "two.sided") return(LR)
  sqrt(LR) * if (isTRUE(fam$under(est))) 1 else -1
}

.disp_lr <- function(object, alternative, alt_missing, B, cores) {
  fam <- .disp_family(object)
  if (fam$boundary) {
    if (!alt_missing && alternative != fam$forced)
      warning("the Poisson value is on the boundary of this family's parameter space; the test is one-sided (alternative = \"",
              fam$forced, "\").")
    alternative <- fam$forced
  } else if (alt_missing) alternative <- "two.sided"
  des <- .disp_design(object)
  size1 <- .disp_first_order_size(des$p, des$n, alternative)
  asym_ok <- !fam$boundary && !des$fe && size1 <= .disp_size_tol
  ll1 <- object$loglik
  null0 <- .disp_null_fit(object); ll0 <- null0$loglik
  neg <- .disp_negative_lr(2 * (ll1 - ll0), fam$boundary && isTRUE(object$support_binding), fam$label)
  LR <- neg$LR; note <- neg$note
  nboot <- if (neg$shortfall) 0L else if (is.null(B)) (if (asym_ok) 0L else 199L) else as.integer(B)
  stat <- if (neg$shortfall) NA_real_ else .disp_signed(LR, fam$est, fam, alternative)
  p <- NA_real_; boot <- NULL; B_ok <- NA_integer_
  ex <- .disp_expand(object)
  remedy <- "the default (B = NULL) calibrates it by parametric bootstrap"
  if (is.null(ex)) {                       # weights that are not whole numbers: there is no replicate data to simulate
    if (nboot > 0L) note <- c(note, "the weights are not whole numbers, so the parametric bootstrap is skipped")
    nboot <- 0L; remedy <- "the parametric bootstrap that calibrates it needs whole-number weights"
  }
  if (nboot > 0L) {
    rate <- null0$rate[ex$index]; trunc <- isTRUE(object$truncated); obj <- ex$object
    one <- function(b) {
      ys <- .disp_simulate(rate, trunc)
      tryCatch(suppressWarnings({
        n0 <- .disp_null_fit(obj, ys)
        a1 <- .disp_alt_fit(obj, ys)
        if (!is.finite(n0$loglik) || !is.finite(a1[["loglik"]]) || a1[["loglik"]] < -1e9) NA_real_
        else .disp_signed(max(2 * (a1[["loglik"]] - n0$loglik), 0), a1[["est"]], fam, alternative)
      }), error = function(e) NA_real_)
    }
    boot <- unlist(.ud_lapply(seq_len(nboot), one, cores = cores))
    ok <- is.finite(boot); B_ok <- sum(ok)
    if (B_ok < 0.8 * nboot)
      warning("dispersion_test: ", nboot - B_ok, " of ", nboot, " bootstrap replicates failed to fit; the p-value uses the ",
              B_ok, " that did.", call. = FALSE)
    if (B_ok > 0L) {
      hit <- if (!fam$boundary && alternative == "over") boot[ok] <= stat + 1e-10 else boot[ok] >= stat - 1e-10
      p <- (1 + sum(hit)) / (1 + B_ok)
    }
  } else if (!neg$shortfall && des$fe) {
    note <- c(note, paste0("unit fixed effects bias the dispersion estimate toward underdispersion, so the asymptotic ",
                           "distribution does not apply; ",
                           if (is.null(ex)) remedy else "B > 0 gives the parametric-bootstrap p-value"))
  } else if (!neg$shortfall) {
    p <- if (fam$boundary) 0.5 * stats::pchisq(LR, 1, lower.tail = FALSE)
         else switch(alternative,
           two.sided = stats::pchisq(LR, 1, lower.tail = FALSE),
           under = stats::pnorm(stat, lower.tail = FALSE),
           over  = stats::pnorm(stat))
    if (!asym_ok)
      note <- c(note, if (fam$boundary)
        paste0("the p-value is asymptotic; at a boundary null the statistic is skewed toward rejection in finite ",
               "samples, and ", remedy)
        else sprintf(paste0("the p-value is asymptotic; with %d mean parameters on %s observations its first-order ",
                            "size at the 5%% level is %.3f, and %s"),
                     as.integer(des$p), format(des$n), size1, remedy))
  }
  structure(list(statistic = c(LR = LR), parameter = c(df = 1), p.value = p, note = note,
                 estimate = stats::setNames(fam$est, fam$par), null.value = stats::setNames(fam$null, fam$par),
                 alternative = switch(alternative, two.sided = "two.sided", under = "underdispersion", over = "overdispersion"),
                 method = paste0("Likelihood-ratio test of equidispersion (", fam$label, " vs Poisson",
                                 if (fam$boundary) "; boundary null" else "",
                                 if (nboot > 0L) sprintf("; parametric bootstrap under the fitted Poisson, %d of %d replicates",
                                                         as.integer(B_ok), nboot)
                                 else if (fam$boundary) "; asymptotic 1/2 chi2_0 + 1/2 chi2_1 mixture"
                                 else "; asymptotic chi-square", ")"),
                 data.name = paste(deparse(object$call$formula), collapse = ""),
                 loglik = c(alternative = ll1, poisson = ll0), boundary = fam$boundary,
                 calibration = if (nboot > 0L) "parametric bootstrap" else "asymptotic",
                 B = nboot, B_ok = B_ok, boot = boot, first_order_size = size1, asymptotic_ok = asym_ok),
            class = c("dispersion_test", "htest"))
}

.disp_auxiliary <- function(object, alternative) {
  if (!inherits(object, "count_reg") || !identical(object$family, "poisson"))
    stop("method = \"auxiliary\" applies to a count_reg(family = \"poisson\") fit.")
  if (isTRUE(object$truncated)) stop("the auxiliary regression is defined for the untruncated Poisson.")
  y <- object$Y; mu <- object$mu; w <- .ud_w1(object$weights, length(y))
  z <- ((y - mu)^2 - y) / mu                      # E[(y-mu)^2 - y] = a * mu^2 under Var = mu + a mu^2
  fit <- stats::lm(z ~ mu + 0, weights = w)
  cf <- summary(fit)$coefficients
  a <- cf[1, 1]; tstat <- cf[1, 3]
  ## lm() counts rows, not weight: rescale the statistic and the df to the weight total
  n <- length(y); sw <- sum(w)
  if (abs(sw - n) > 1e-8) tstat <- tstat * sqrt((sw - 1) / (n - 1))
  p <- switch(alternative,
    two.sided = 2 * stats::pnorm(-abs(tstat)),
    under = stats::pnorm(tstat),                  # a < 0
    over  = stats::pnorm(tstat, lower.tail = FALSE))
  structure(list(statistic = c(t = tstat), parameter = c(df = sw - 1), p.value = p,
                 estimate = c(a = a), null.value = c(a = 0),
                 alternative = switch(alternative, two.sided = "two.sided", under = "underdispersion", over = "overdispersion"),
                 method = "Cameron-Trivedi auxiliary regression test of equidispersion (Var = mu + a mu^2)",
                 data.name = paste(deparse(object$call$formula), collapse = ""), boundary = FALSE),
            class = c("dispersion_test", "htest"))
}

#' @export
print.dispersion_test <- function(x, digits = 4, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data: ", x$data.name, "\n", sep = "")
  cat(sprintf("%s = %.4f", names(x$statistic), x$statistic),
      if (!is.null(x$loglik)) sprintf("  (logLik: fitted family = %.2f, Poisson = %.2f)", x$loglik[1], x$loglik[2]),
      sprintf(",  p-value = %s\n", format.pval(x$p.value, digits = digits)), sep = "")
  cat("alternative hypothesis: ", x$alternative, "\n", sep = "")
  for (nt in x$note) cat("note: ", nt, "\n", sep = "")
  cat(sprintf("estimate: %s = %.4f  (Poisson value %s)\n", names(x$estimate), x$estimate, format(x$null.value)))
  invisible(x)
}
