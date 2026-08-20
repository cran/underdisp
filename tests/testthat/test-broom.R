test_that("broom tidiers work on a cpb fit", {
  skip_if_not_installed("broom")
  set.seed(1); x <- rnorm(400)
  N <- pmax(round(exp(1.6 + 0.4 * x) / 0.5), 1); y <- rbinom(400, N, 0.5)
  fit <- cpb(y ~ x, data = data.frame(y = y, x = x), truncated = FALSE, se = "none")

  td <- broom::tidy(fit)
  expect_true(all(c("term", "estimate", "std.error", "statistic", "p.value") %in% names(td)))
  expect_equal(nrow(td), 2L)

  gl <- broom::glance(fit)
  expect_true(all(c("alpha", "logLik", "AIC", "df", "nobs") %in% names(gl)))
  expect_equal(gl$nobs, length(y))

  ag <- broom::augment(fit)
  expect_true(all(c(".fitted", ".resid", ".ceiling") %in% names(ag)))
  expect_equal(nrow(ag), length(y))
})
