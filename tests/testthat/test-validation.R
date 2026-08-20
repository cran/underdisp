test_that("estimators reject malformed responses with clear messages", {
  set.seed(1); d <- data.frame(y = rpois(120, 2) + 1L, x = rnorm(120))
  expect_error(cpb(I(y - 5) ~ x, d, truncated = FALSE, se = "none"), "non-negative")
  expect_error(cpb(I(y + 0.5) ~ x, d, truncated = FALSE, se = "none"), "non-negative")
  d0 <- data.frame(y = rpois(120, 1), x = rnorm(120))          # contains zeros
  expect_error(cpb(y ~ x, d0, truncated = TRUE, se = "none"), "requires all Y")
  expect_error(hurdle_cpb(I(y - 5) ~ x, d, participation = ~ x), "nonnegative")
})

test_that("cluster and fixed-effects arguments validate their columns", {
  set.seed(1); d <- data.frame(y = rpois(120, 2) + 1L, x = rnorm(120))
  expect_error(cpb(y ~ x, d, truncated = TRUE, se = "bootstrap", B = 3, cluster = "nope"), "cluster")
  expect_error(cpb_fe(y ~ x, d, fe = "nope"), "fe")
  expect_error(cpb_fe(y ~ 1, d, fe = "x"), "covariate")
})

test_that("predict requires newdata to contain the model variables", {
  set.seed(1); x <- rnorm(120); d <- data.frame(y = rpois(120, 2) + 1L, x = x)
  f <- cpb(y ~ x, d, truncated = FALSE, se = "none")
  expect_error(predict(f, newdata = data.frame(w = 1)), "missing required variable")
  expect_length(predict(f, newdata = data.frame(x = c(-1, 0, 1))), 3L)   # correct newdata works
})

test_that("estimators survive unusual but valid inputs", {
  set.seed(1)
  d <- data.frame(y = rcpb(250, exp(1 + 0.5 * rnorm(250)), 0.5, truncated = FALSE),
                  x = rnorm(250), fx = factor(sample(letters[1:3], 250, TRUE)))
  expect_s3_class(cpb(y ~ 1, d, truncated = FALSE, se = "none"), "cpb")     # intercept-only
  expect_s3_class(cpb(y ~ fx, d, truncated = FALSE, se = "none"), "cpb")    # factor covariate
  expect_s3_class(zi_cpb(y ~ x, d, zero = ~ 1), "zi_cpb")
  d2 <- d; d2$x[1:6] <- NA
  expect_s3_class(suppressWarnings(cpb(y ~ x, d2, truncated = FALSE, se = "none")), "cpb")  # NA rows dropped
})
