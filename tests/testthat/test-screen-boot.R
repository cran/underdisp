test_that("parametric-bootstrap threshold is sane and detects underdispersion", {
  set.seed(42)
  n <- 400; x <- rnorm(n)
  N <- pmax(round(exp(1.2 + 0.4 * x) / 0.5), 1); y <- rbinom(n, N, 0.5)
  s <- ud_screen(y ~ x, data = data.frame(y = y, x = x), run_cpb = FALSE,
                 run_gp = FALSE, run_comp = FALSE,
                 ztp_threshold = "bootstrap", ztp_boot_B = 59)
  expect_identical(s$ztp_threshold_method, "bootstrap")
  expect_true(is.finite(s$ztp_threshold))
  expect_lt(s$ztp_threshold, 1)
  expect_gt(s$ztp_threshold, 0.5)
  expect_identical(s$verdict_atrisk, "UNDERDISPERSED")
})

test_that("bootstrap threshold does not flag an equidispersed null", {
  set.seed(7)
  n <- 400; x <- rnorm(n); y <- rpois(n, exp(1.0 + 0.3 * x))
  s <- ud_screen(y ~ x, data = data.frame(y = y, x = x), run_cpb = FALSE,
                 run_gp = FALSE, run_comp = FALSE,
                 ztp_threshold = "bootstrap", ztp_boot_B = 59)
  expect_false(identical(s$verdict_atrisk, "UNDERDISPERSED"))
})

test_that("calibrated default is unchanged by the new arguments", {
  set.seed(1); x <- rnorm(500)
  N <- pmax(round(exp(1.4 + 0.4 * x) / 0.6), 1); y <- rbinom(500, N, 0.6)
  d <- data.frame(y = y, x = x)
  a <- ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE)
  b <- ud_screen(y ~ x, data = d, run_cpb = FALSE, run_gp = FALSE, run_comp = FALSE,
                 ztp_threshold = "calibrated")
  expect_equal(a$pearson_ztp, b$pearson_ztp)
  expect_equal(a$ztp_threshold, b$ztp_threshold)
  expect_identical(a$ztp_threshold_method, "calibrated")
})
