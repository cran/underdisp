cal_data <- function(seed, n = 120) {
  set.seed(seed); x <- rnorm(n)
  data.frame(y = rcpb(n, exp(1.5 + 0.4 * x), 0.5, truncated = TRUE), x = x)
}

test_that("calibrate_alpha() stores the bootstrap signed roots and the accessors use them", {
  skip_on_cran()
  fit <- cpb(y ~ x, cal_data(41), se = "none")
  first <- alpha_confint(fit)
  expect_match(attr(first, "method"), "first-order")
  set.seed(5); cal <- calibrate_alpha(fit, B = 39, seed = 1); u1 <- runif(1)
  set.seed(5); u2 <- runif(1)
  expect_equal(u1, u2)                                         # a seeded calibration leaves the caller's stream alone
  expect_identical(cal$alpha.cal$B, 39L)
  expect_gte(cal$alpha.cal$B_ok, 35)
  expect_equal(cal$coefficients, fit$coefficients); expect_equal(cal$alpha, fit$alpha)
  ci <- alpha_confint(cal, level = 0.90)
  expect_match(attr(ci, "method"), "calibrated")
  expect_true(ci[["lower"]] <= fit$alpha && fit$alpha <= ci[["upper"]])
  expect_true(all(is.finite(ci)))
  expect_identical(calibrate_alpha(fit, B = 39, seed = 1)$alpha.cal$r, cal$alpha.cal$r)
  expect_error(alpha_confint(cal, level = 0.99), "larger B")    # 39 replicates cannot serve a 99% interval
  expect_output(print(summary(cal)), "calibrated profile")
  expect_output(print(summary(fit)), "first-order profile")
  ic <- implied_ceiling(cal, newdata = data.frame(x = 0), level = 0.90)
  expect_true(ic$lower <= ic$ceiling && ic$ceiling <= ic$upper)
})

test_that("the calibration does not depend on the number of worker processes", {
  skip_on_cran()
  fit <- cpb(y ~ x, cal_data(42, n = 80), se = "none")
  r1 <- calibrate_alpha(fit, B = 19, cores = 1, seed = 3)$alpha.cal$r
  r2 <- calibrate_alpha(fit, B = 19, cores = 2, seed = 3)$alpha.cal$r
  expect_equal(r1, r2)
})

test_that("whole-number weights are frequency weights in the calibration", {
  skip_on_cran()
  d <- cal_data(43, n = 60); d$w <- 2
  fw <- cpb(y ~ x, d, se = "none", weights = "w")
  fe <- cpb(y ~ x, d[rep(seq_len(nrow(d)), each = 2), ], se = "none")
  rw <- calibrate_alpha(fw, B = 19, seed = 8)$alpha.cal$r
  re <- calibrate_alpha(fe, B = 19, seed = 8)$alpha.cal$r
  expect_equal(rw, re, tolerance = 1e-5)
})

test_that("calibrate_alpha() refuses the cases it cannot serve", {
  d <- cal_data(44)
  fit <- cpb(y ~ x, d, se = "none")
  expect_error(calibrate_alpha(fit, B = 5), "'B'")
  old <- fit; old["profile"] <- list(NULL)
  expect_error(calibrate_alpha(old, B = 19), "earlier version")
  guard <- fit; guard$support_binding <- TRUE
  expect_error(calibrate_alpha(guard, B = 19), "max.support")
  d$w <- runif(nrow(d), 0.5, 1.5)
  expect_error(calibrate_alpha(cpb(y ~ x, d, se = "none", weights = "w"), B = 19), "whole numbers")
  expect_error(calibrate_alpha(lm(y ~ x, d)), "cpb\\(\\) fits")
})
