## Offset support for the ZI mixtures. The decisive check is reparameterization
## invariance: a CONSTANT log-offset c must shift the count intercept by
## exactly -c, leave slopes/alpha/zero coefficients essentially unchanged, and
## leave the maximized log-likelihood identical.

test_that("zi_cpb offset is a reparameterization (constant-offset invariance)", {
  skip_on_cran()
  set.seed(21); n <- 300; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, exp(1.2 + 0.3 * x), 0.5, pi = plogis(-0.6 + 0.8 * z))
  d <- data.frame(y = y, x = x, z = z)
  f0 <- zi_cpb(y ~ x, d, zero = ~ z)
  f1 <- zi_cpb(y ~ x, d, zero = ~ z, offset = rep(log(2), n))
  expect_equal(f1$loglik, f0$loglik, tolerance = 1e-3)
  expect_equal(unname(f1$coefficients["(Intercept)"]),
               unname(f0$coefficients["(Intercept)"]) - log(2), tolerance = 0.02)
  expect_equal(unname(f1$coefficients["x"]), unname(f0$coefficients["x"]), tolerance = 0.02)
  expect_equal(f1$alpha, f0$alpha, tolerance = 0.03)
  ## lambda_full includes the offset, so fitted marginals agree
  expect_equal(f1$lambda_full, f0$lambda_full, tolerance = 0.05)
})

test_that("zi_gec offset is a reparameterization", {
  skip_on_cran()
  set.seed(22); n <- 300; x <- rnorm(n); z <- rnorm(n)
  y <- rzicpb(n, exp(1.2 + 0.3 * x), 0.5, pi = plogis(-0.6 + 0.8 * z))
  d <- data.frame(y = y, x = x, z = z)
  f0 <- zi_gec(y ~ x, d, zero = ~ z)
  f1 <- zi_gec(y ~ x, d, zero = ~ z, offset = rep(log(2), n))
  expect_equal(f1$loglik, f0$loglik, tolerance = 1e-3)
  expect_equal(unname(f1$coefficients["(Intercept)"]),
               unname(f0$coefficients["(Intercept)"]) - log(2), tolerance = 0.02)
  expect_equal(f1$delta, f0$delta, tolerance = 0.05)
})

test_that("offset guards: fe and em paths refuse loudly", {
  skip_on_cran()
  d <- data.frame(y = c(0L, 1L, 2L, 0L, 3L, 1L), x = rnorm(6),
                  u = rep(1:2, each = 3))
  expect_error(zi_cpb(y ~ x, d, zero = ~ 1, fe = "u", offset = rep(0.1, 6)),
               "not yet supported")
  expect_error(zi_cpb(y ~ x, d, zero = ~ 1, method = "em", offset = rep(0.1, 6)),
               "method = \"ml\" only")
})
