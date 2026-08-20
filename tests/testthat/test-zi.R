test_that("zi_cpb recovers intensity and zero-inflation parameters", {
  skip_on_cran()  # EM is a little slow
  set.seed(1)
  n <- 700; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z))
  fit <- zi_cpb(y ~ x, data = data.frame(y = y, x = x, z = z), zero = ~ z)
  expect_s3_class(fit, "zi_cpb")

  expect_true(fit$coefficients[["x"]] > 0.3 && fit$coefficients[["x"]] < 0.7)
  expect_true(fit$alpha > 0.3 && fit$alpha < 0.7)
  expect_true(fit$zero_coef[["z"]] > 0.3)      # positive z effect on structural zeros

  nd <- data.frame(x = 0, z = 0)
  pz <- predict(fit, nd, type = "zero")
  m  <- predict(fit, nd, type = "response")
  expect_true(pz > 0 && pz < 1 && m > 0)
})

test_that("rzicpb produces excess zeros", {
  set.seed(2)
  y <- rzicpb(500, lambda = 4, alpha = 0.5, pi = 0.4)
  expect_true(mean(y == 0) > 0.3)
})
