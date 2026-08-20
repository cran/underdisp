test_that("cpb recovers parameters and computes quantities of interest", {
  set.seed(1)
  n <- 300; x <- rnorm(n)
  N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1)
  y <- rbinom(n, N, 0.5)                       # underdispersed, var/mean approx 0.5
  d <- data.frame(y = y, x = x); d <- d[d$y > 0, ]

  fit <- cpb(y ~ x, data = d, se = "none")
  expect_s3_class(fit, "cpb")
  expect_true(fit$alpha > 0 && fit$alpha < 1)                 # underdispersion detected
  expect_equal(unname(coef(fit)[["x"]]), 0.5, tolerance = 0.2)
  expect_equal(as.numeric(logLik(fit)), fit$loglik)
  expect_identical(nobs(fit), nrow(d))

  # predictions and the implied ceiling
  expect_length(predict(fit), nrow(d))
  expect_true(all(fit$ceiling >= d$y))                        # ceiling respected in-sample
  ic <- implied_ceiling(fit, newdata = data.frame(x = 0))
  expect_true(is.finite(ic$ceiling))
  p5 <- predict(fit, newdata = data.frame(x = 0), type = "prob", at = 5)
  expect_true(p5 >= 0 && p5 <= 1)
})

test_that("cpb bootstrap yields confidence intervals", {
  skip_on_cran()
  set.seed(2)
  n <- 200; x <- rnorm(n)
  N <- pmax(round(exp(1.6 + 0.5 * x) / 0.5), 1); y <- rbinom(n, N, 0.5)
  d <- data.frame(y = y, x = x); d <- d[d$y > 0, ]
  fit <- suppressWarnings(cpb(y ~ x, data = d, se = "bootstrap", B = 60))
  ci <- confint(fit)
  expect_true("alpha" %in% rownames(ci))
  expect_true(all(ci[, 1] <= ci[, 2]))
  fd <- first_difference(fit, "x", from = -1, to = 1, quantity = "mean")
  expect_true(fd[["lower"]] <= fd[["upper"]])
})

test_that("cpb validates its input", {
  expect_error(cpb(y ~ x, data = data.frame(y = c(1, -1), x = 1:2), se = "none"))
  expect_error(cpb(y ~ x, data = data.frame(y = c(0, 1), x = 1:2),
                   truncated = TRUE, se = "none"))
})
