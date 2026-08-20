## Regression tests for the FE zero-inflated CPB estimators. The seeded EM is the
## vehicle for a comparison symmetric to a hurdle with participation fixed effects:
## it uses the zeros (unlike the two-step positives-only offset), supports per-unit
## inflation fixed effects (zero_fe), and does not collapse the structural-zero
## probability toward a degenerate pi -> 0.

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
})
