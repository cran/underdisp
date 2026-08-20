test_that("cpb_fe recovers coefficients with high-dimensional fixed effects", {
  set.seed(3)
  d <- do.call(rbind, lapply(1:80, function(i) {
    x <- rnorm(15)
    N <- pmax(round(exp(rnorm(1, 0, 0.5) + 0.5 * x) / 0.5), 1)
    data.frame(unit = i, x = x, y = rbinom(15, N, 0.5))
  }))
  fit <- cpb_fe(y ~ x, data = d, fe = "unit")
  expect_s3_class(fit, "cpb_fe")
  # recovers a positive slope of the right sign and rough magnitude (finite-T)
  expect_true(fit$coefficients[["x"]] > 0.25 && fit$coefficients[["x"]] < 0.8)
  expect_true(fit$alpha > 0 && fit$alpha < 1)
  expect_identical(fit$n_units, 80L)
  expect_length(fit$fe, 80L)
})

test_that("cpb_fe requires a covariate and a valid fe column", {
  d <- data.frame(unit = rep(1:5, each = 4), x = rnorm(20), y = rpois(20, 3) + 1L)
  expect_error(cpb_fe(y ~ 1, data = d, fe = "unit"))
  expect_error(cpb_fe(y ~ x, data = d, fe = "nope"))
})
