test_that("zi_test rejects incorrect model types", {
  skip_on_cran()
  set.seed(1); n <- 300; x <- rnorm(n)
  y <- rzicpb(n, exp(1 + 0.4 * x), 0.5, 0.3); d <- data.frame(y = y, x = x)
  cf <- cpb(y ~ x, data = d, truncated = FALSE, se = "none")
  zf <- zi_cpb(y ~ x, data = d, zero = ~ 1)
  expect_s3_class(zi_test(cf, zf), "zi_test")            # valid: nest + zero-inflated
  expect_s3_class(zi_test(zf, cf), "zi_test")            # order-independent
  expect_error(zi_test(cf, cf), "zero-inflated")         # two non-inflated models
  expect_error(zi_test(zf, zf), "Both")                  # two inflated models
})

test_that("zi_test still works from a formula", {
  skip_on_cran()
  set.seed(2); n <- 400; x <- rnorm(n)
  y <- rzicpb(n, exp(1.2 + 0.4 * x), 0.5, 0.3)
  tt <- zi_test(y ~ x, data = data.frame(y = y, x = x))
  expect_s3_class(tt, "zi_test")
  expect_true(tt$LR >= 0)
})

test_that("compare_models validates inputs and ranks models", {
  skip_on_cran()
  set.seed(3); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- rhurdle_cpb(n, exp(1 + 0.5 * x), 0.5, plogis(-0.2 + 0.8 * z))
  d <- data.frame(y = y, x = x, z = z)
  h  <- hurdle_cpb(y ~ x, data = d, participation = ~ z)
  zf <- zi_cpb(y ~ x, data = d, zero = ~ z)
  cm <- compare_models(hurdle = h, zi = zf)
  expect_true(all(c("df", "logLik", "AIC", "BIC", "logscore", "rps") %in% colnames(cm)))
  expect_equal(nrow(cm), 2L)
  expect_error(compare_models(h), "at least two")
  expect_error(compare_models(h, lm(y ~ x, d)), "fitted models from this package")
})
