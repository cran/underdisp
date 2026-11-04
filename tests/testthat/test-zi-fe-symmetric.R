## Regression tests for the FE zero-inflated CPB estimators. The seeded EM is the
## vehicle for a comparison symmetric to a hurdle with participation fixed effects:
## it uses the zeros (unlike the two-step positives-only offset), supports per-unit
## inflation fixed effects (zero_fe), and does not collapse the structural-zero
## probability toward a degenerate pi -> 0.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

test_that("EM FE-ZI recovers parameters and does not collapse pi -> 0", {
  set.seed(3); nu <- 25; per <- 16; n <- nu * per
  u <- factor(rep(seq_len(nu), each = per)); x <- rnorm(n); z <- rnorm(n)
  fe <- rnorm(nu, 1.0, 0.4)[as.integer(u)]; lam <- exp(fe + 0.5 * x); pit <- plogis(-0.5 + 0.8 * z)
  y <- rzicpb(n, lam, 0.5, pit); d <- data.frame(y = y, x = x, z = z, u = u)
  em <- zi_cpb(y ~ x, d, zero = ~ z, fe = "u", method = "em", se = "none")
  expect_gt(mean(em$pi_full), 0.05)                       # no degenerate collapse to pi = 0
  expect_equal(unname(em$coefficients["x"]), 0.5, tolerance = 0.25)
  expect_gt(unname(em$zero_coef["z"]), 0.2)               # right sign/magnitude (small-sample noisy)
})

test_that("FE-ZI accepts per-unit inflation fixed effects (zero_fe) for a symmetric comparison", {
  skip_on_cran()
  set.seed(4); nu <- 25; per <- 16; n <- nu * per
  u <- factor(rep(seq_len(nu), each = per)); x <- rnorm(n)
  fe <- rnorm(nu, 1.0, 0.4)[as.integer(u)]; zfe <- rnorm(nu, -0.3, 0.6)[as.integer(u)]
  lam <- exp(fe + 0.5 * x); pit <- plogis(zfe)             # per-unit structural-zero rate
  y <- rzicpb(n, lam, 0.5, pit); d <- data.frame(y = y, x = x, u = u)
  ems <- zi_cpb(y ~ x, d, zero = ~ 1, fe = "u", zero_fe = "u", method = "em", se = "none")
  expect_s3_class(ems, "zi_cpb")
  expect_true(is.finite(AIC(ems)))
  expect_gt(mean(ems$pi_full), 0.05)
})

test_that("EM FE-ZI handles units with no positive counts (all-zero units)", {
  ## Regression: units whose counts are all zero never enter the count fit, so
  ## fe_hat used to come back shorter than the number of levels and the fit
  ## crashed ('names' attribute must be the same length as the vector). The fix
  ## expands fe_hat to every level with the EM's mean-effect fallback and counts
  ## only the genuinely estimated unit effects toward df.
  skip_on_cran()
  set.seed(5); nu <- 20; per <- 12; n <- nu * per
  u <- factor(rep(seq_len(nu), each = per)); x <- rnorm(n)
  fe <- rnorm(nu, 1.0, 0.4)[as.integer(u)]
  pit <- rep(0.25, n); pit[as.integer(u) <= 5] <- 1        # first 5 units: all structural zeros
  y <- rzicpb(n, exp(fe + 0.5 * x), 0.5, pit)
  expect_true(all(tapply(y, u, max)[1:5] == 0))
  d <- data.frame(y = y, x = x, u = u)
  em <- zi_cpb(y ~ x, d, zero = ~ 1, fe = "u", zero_fe = "u", method = "em", se = "none")
  expect_s3_class(em, "zi_cpb")
  expect_length(em$fe_hat, nu)                             # expanded to every level
  expect_true(all(is.finite(em$fe_hat)))
  expect_true(all(is.finite(em$lambda_full)))
  expect_true(is.finite(AIC(em)))
  full_count <- length(em$coefficients) + nu + 1L + length(em$zero_coef)
  expect_lt(em$df, full_count)                             # df excludes never-estimated unit effects
  ## On King's support (above) those units are left out before the EM runs (.zi_zero_units); the EM's own
  ## mean-effect expansion is the floor support's path, as before 0.2.0
  ef <- suppressWarnings(zi_cpb(y ~ x, d, zero = ~ 1, fe = "u", zero_fe = "u", method = "em", se = "none",
                                support = "floor"))
  expect_length(ef$fe_hat, nu)
  expect_true(all(is.finite(ef$fe_hat)))
  expect_false(any(ef$fe_estimated[1:5]))
  expect_lt(ef$df, length(ef$coefficients) + nu + 1L + length(ef$zero_coef))
})

test_that("the EM reports convergence as it happened on King's support", {
  set.seed(6); dz <- data.frame(x = rnorm(300)); dz$y <- rzicpb(300, exp(1 + 0.4 * dz$x), 0.5, 0.3)
  expect_false(suppressWarnings(zi_cpb(y ~ x, dz, zero = ~ 1, method = "em", maxit = 2))$converged)
  expect_true(suppressWarnings(zi_cpb(y ~ x, dz, zero = ~ 1, method = "em", maxit = 2, support = "floor"))$converged)
})

test_that("zi_cpb(fe, ml) is a joint maximum on King's support, whatever the covariates' centring", {
  skip_on_cran()
  set.seed(5)
  pz <- do.call(rbind, lapply(1:40, function(i) { x <- rnorm(8); yr <- 2000 + 1:8
    data.frame(unit = i, x = x, year = yr, y = rzicpb(8, exp(0.5 + rnorm(1, 0, .5) + 0.4 * x + 0.05 * (yr - 2004)), 0.5, 0.3)) }))
  f0 <- suppressWarnings(zi_cpb(y ~ x + yr, transform(pz, yr = year - 2000), zero = ~ 1, fe = "unit"))
  f2 <- suppressWarnings(zi_cpb(y ~ x + yr, transform(pz, yr = year - 2000, x = x + 5), zero = ~ 1, fe = "unit"))
  em <- suppressWarnings(zi_cpb(y ~ x + yr, transform(pz, yr = year - 2000), zero = ~ 1, fe = "unit", method = "em"))
  expect_lt(abs(f2$loglik - f0$loglik), 0.1)                            # the plug-in alone differed by about 3
  expect_gt(f0$loglik, em$loglik - 0.1)                                  # the plug-in alone lay about 7 below
})