## Numerical accuracy at the edges the third review wave probed: the EWP's kinked likelihood and its standard errors,
## the gamma count's log-probabilities and moments, the families' variances, the support caps, the mean-parameterized
## COM-Poisson's root, and units with no positive count under uncentred covariates.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

test_that("EWP standard errors are the inverse expected information, from the family's own scores", {
  set.seed(101); n <- 200; x <- rnorm(n); lam <- exp(1 + 0.4 * x)
  d <- data.frame(y = vapply(lam, function(l) rewp(1, l, 1, 0.6), numeric(1)), x = x)
  m <- suppressWarnings(count_reg(y ~ x, d, family = "ewp3"))
  b1 <- m$theta[1]; b2 <- m$theta[2]; I <- matrix(0, 4, 4)
  for (i in seq_len(n)) {                     # the term derivatives in (eta, beta1, beta2), centred under the fit
    l <- m$mu[i]; k <- 0:80; p <- dewp(k, l, b1, b2); lo <- k <= l
    dt <- cbind(k - l + l * ifelse(lo, -b1, b2), ifelse(lo, -(l - k), 0), ifelse(lo, 0, -(k - l)))
    C <- crossprod(sweep(dt, 2, colSums(dt * p)) * sqrt(p))
    J <- rbind(c(1, x[i], 0, 0), c(0, 0, 1, 0), c(0, 0, 0, 1))
    I <- I + crossprod(J, C %*% J)
  }
  ## to 1%: the maximum of a kinked likelihood puts some centralities exactly at an integer (as a least-absolute-
  ## deviations fit interpolates data points), where a support point's score has two one-sided values; this
  ## reference takes the left one and the package their average
  expect_equal(unname(sqrt(diag(m$vcov_full))), sqrt(diag(solve(I))), tolerance = 1e-2)
  expect_true(all(is.finite(suppressWarnings(count_reg(y ~ x, d, family = "ewp3", se = "robust"))$se.beta)))
})

test_that("the gamma count's log-probabilities do not underflow, and its moments refuse at the cap", {
  expect_equal(dgammacount(0, 400, 2, log = TRUE), log(801) - 800, tolerance = 1e-10)   # Q(2, 800) = 801 e^-800
  expect_equal(dgammacount(0:6, 3, 2), exp(dgammacount(0:6, 3, 2, log = TRUE)), tolerance = 1e-14)
  gm <- underdisp:::.gc_moments(c(5e4, 1e5), 2)
  expect_true(all(is.finite(gm[, 1]))); expect_true(all(is.na(gm[, 2])))
})

test_that("variances are taken about the mean, and sums the cap cuts while rising are refused", {
  p <- 40 / 41                                  # the COM-Poisson at nu = 1e4 and rate 40 is a Bernoulli
  expect_equal(underdisp:::cmp_moments_cpp(log(40), 1e4)[, 2], p * (1 - p), tolerance = 1e-8)
  expect_true(all(is.na(underdisp:::ewp_norm_cpp(log(9e4), -0.5, -0.5))))
  expect_equal(underdisp:::cmp_loglambda_cpp(0.05, 2000), log(0.05 / 0.95), tolerance = 1e-6)
})

test_that("units with no positive count stay at rates near zero under an uncentred covariate", {
  set.seed(3)
  pf <- do.call(rbind, lapply(1:40, function(i) { x <- rnorm(6); yr <- 2001:2006
    data.frame(unit = i, x = x, year = yr, y = rcpb(6, exp(0.6 + rnorm(1, 0, .4) + 0.3 * x + 0.04 * (yr - 2003)), 0.5)) }))
  zr <- pf$unit %in% c(4, 17, 30); pf$y[zr] <- 0L
  for (fn in list(cpb_fe, gec_fe)) {
    fu <- suppressWarnings(fn(y ~ x + year, pf, fe = "unit"))
    fc <- suppressWarnings(fn(y ~ x + I(year - 2003), pf, fe = "unit"))
    expect_false(anyNA(fitted(fu))); expect_lt(max(fitted(fu)[zr]), 1e-10)
    expect_equal(unname(fu$coefficients[["x"]]), unname(fc$coefficients[["x"]]), tolerance = 1e-4)
  }
  zf <- suppressMessages(zi_cpb(y ~ x + year, pf, zero = ~ 1, fe = "unit"))
  expect_lt(max(zf$lambda_full[zr]), 1e-10); expect_false(anyNA(fitted(zf)))
})

test_that("count_reg's distinct-row evaluation gives the full data's fit and covariances", {
  set.seed(6); n <- 300; g <- factor(sample(letters[1:3], n, TRUE)); cl <- sample(1:30, n, TRUE)
  d <- data.frame(y = rpois(n, exp(0.5 + 0.3 * (g == "b"))), g = g, cl = cl)
  f <- count_reg(y ~ g, d, family = "poisson", se = "robust")
  gl <- glm(y ~ g, poisson, d)
  expect_equal(unname(coef(f)), unname(coef(gl)), tolerance = 1e-7)
  X <- model.matrix(gl); r <- d$y - fitted(gl); Bm <- solve(crossprod(X * sqrt(fitted(gl))))
  expect_equal(unname(f$vcov), unname(Bm %*% crossprod(X * r) %*% Bm), tolerance = 1e-5)
  fc <- count_reg(y ~ g, d, family = "poisson", se = "cluster", cluster = "cl")
  S <- rowsum(X * r, d$cl)
  expect_equal(unname(fc$vcov), unname(30 / 29 * Bm %*% crossprod(S) %*% Bm), tolerance = 1e-5)
  expect_length(f$fitted.values, n)
})
