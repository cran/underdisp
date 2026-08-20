## Regression tests for held-out / cross-validated scoring (score(newdata=), cv_score).

test_that("score(fit, newdata = training data) reproduces the in-sample score (quick)", {
  set.seed(1); n <- 200; x <- rnorm(n)
  d <- data.frame(y = rcpb(n, exp(1 + 0.4 * x), 0.5, truncated = FALSE), x = x)
  same <- function(fit, dat) expect_equal(unname(score(fit, newdata = dat)), unname(score(fit)), tolerance = 1e-8)
  same(cpb(y ~ x, d, truncated = FALSE, se = "none"), d)
  same(count_reg(y ~ x, d, "poisson", se = "none"), d)
})

test_that("score(fit, newdata = training data) reproduces the in-sample score", {
  skip_on_cran()
  set.seed(1); n <- 400; x <- rnorm(n); z <- rnorm(n); u <- factor(sample(15, n, replace = TRUE))
  yc <- rcpb(n, exp(1 + 0.4 * x), 0.5, truncated = FALSE); d <- data.frame(y = yc, x = x, z = z, u = u)
  yh <- ifelse(rbinom(n, 1, plogis(0.3 + 0.6 * z)) == 1, 0L, 1L + rpois(n, exp(0.6 + 0.4 * x)))
  dh <- data.frame(y = yh, x = x, z = z, u = u)
  yz <- ifelse(rbinom(n, 1, plogis(-0.4 + 0.7 * z)) == 1, 0L, rzicpb(n, exp(1.1 + 0.4 * x), 0.5, 0))
  dz <- data.frame(y = yz, x = x, z = z, u = u)

  same <- function(fit, dat) expect_equal(unname(score(fit, newdata = dat)), unname(score(fit)), tolerance = 1e-8)
  same(cpb(y ~ x, d, truncated = FALSE, se = "none"), d)
  same(cpb_fe(y ~ x, d, fe = "u", se = "none"), d)
  same(hurdle_cpb(y ~ x, dh, participation = ~ z, se = "none"), dh)
  same(hurdle_cpb(y ~ x, dh, participation = ~ z, fe = "u", se = "none"), dh)
  same(zi_cpb(y ~ x, dz, zero = ~ z, se = "none"), dz)
  same(zi_cpb(y ~ x, dz, zero = ~ z, fe = "u", se = "none"), dz)
  same(count_reg(y ~ x, d, "poisson", se = "none"), d)
  same(hurdle_count(y ~ x, dh, "poisson", participation = ~ z), dh)
  same(zi_count(y ~ x, dz, "poisson", zero = ~ z), dz)
  same(gec(y ~ x, d, se = "none"), d)
  same(gec_fe(y ~ x, d, fe = "u", se = "none"), d)   # gec_fe inherits gec: must dispatch to the specific branch
  same(hurdle_gec(y ~ x, dh, participation = ~ z, se = "none"), dh)
  same(zi_gec(y ~ x, dz, zero = ~ z, se = "none"), dz)
})

test_that("held-out scoring tolerates fixed-effect units unseen in training", {
  skip_on_cran()
  ## Regression: units present only in a held-out fold (dissolved or newly created
  ## states in contiguous year-block CV) used to crash the participation/inflation
  ## dummy design and drop the whole fold. Unseen units now score at the average
  ## estimated unit effect in the binary equation (and the mean fixed effect in
  ## the intensity, as before).
  set.seed(3); n <- 480; x <- rnorm(n); z <- rnorm(n)
  u <- factor(rep(sprintf("u%02d", 1:24), each = 20))
  pu <- plogis(0.5 + 0.6 * z + rnorm(24, 0, 0.7)[as.integer(u)])
  yh <- ifelse(rbinom(n, 1, pu) == 1, 1L + rpois(n, exp(0.6 + 0.4 * x)), 0L)
  dh <- data.frame(y = yh, x = x, z = z, u = u)
  yz <- ifelse(rbinom(n, 1, plogis(-0.6 + rnorm(24, 0, 0.7)[as.integer(u)])) == 1, 0L,
               rzicpb(n, exp(1.1 + 0.4 * x), 0.5, 0))
  dz <- data.frame(y = yz, x = x, z = z, u = u)

  tr_h <- dh[dh$u != "u01", ]; te_h <- dh[dh$u == "u01", ]
  h  <- hurdle_cpb(y ~ x, tr_h, participation = ~ z, fe = "u", part_fe = "u", se = "none")
  sh <- score(h, newdata = te_h)
  expect_true(all(is.finite(sh)))
  ## with every level seen, the tolerant path must reproduce the in-sample score
  expect_equal(unname(score(h, newdata = tr_h)), unname(score(h)), tolerance = 1e-8)

  tr_z <- dz[dz$u != "u01", ]; te_z <- dz[dz$u == "u01", ]
  zf <- zi_cpb(y ~ x, tr_z, zero = ~ z, fe = "u", zero_fe = "u", method = "em", se = "none")
  sz <- score(zf, newdata = te_z)
  expect_true(all(is.finite(sz)))
  expect_equal(unname(score(zf, newdata = tr_z)), unname(score(zf)), tolerance = 1e-8)
})

test_that("cv_score is genuinely out of sample (>= in-sample) and returns per-obs scores", {
  skip_on_cran()
  set.seed(2); n <- 500; x <- rnorm(n); z <- rnorm(n); u <- factor(sample(20, n, replace = TRUE))
  yh <- ifelse(rbinom(n, 1, plogis(0.3 + 0.6 * z)) == 1, 0L, 1L + rpois(n, exp(0.6 + 0.4 * x)))
  dh <- data.frame(y = yh, x = x, z = z, u = u)
  ins <- score(hurdle_cpb(y ~ x, dh, participation = ~ z, fe = "u", se = "none"))["logscore"]
  cv  <- cv_score(function(tr) hurdle_cpb(y ~ x, tr, participation = ~ z, fe = "u", se = "none"), dh, k = 5)
  expect_s3_class(cv, "cv_score")
  expect_length(cv$logscore_i, n)
  expect_gte(cv$logscore, ins - 1e-8)                 # held-out cannot beat in-sample in expectation
})
