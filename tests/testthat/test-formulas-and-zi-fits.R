## Formulas and new data the count classes must take as given, and the zero-inflated fits' reported quantities:
## transformed responses and '.' with fixed effects, character covariates in new data, the observed-data
## log-likelihood of the EM fits, the degeneracy guard, zi_test()'s nests, weighted nobs, and the support guard.

test_that("transformed responses and '.' survive the fixed-effects formula", {
  set.seed(7); n <- 300; x <- rnorm(n); z <- rnorm(n); u <- factor(sample(sprintf("u%02d", 1:15), n, TRUE))
  d <- data.frame(y = rpois(n, exp(1 + 0.3 * x)), x, z, u)
  a <- count_reg(pmin(y, 3) ~ x, d, family = "poisson", fe = "u")
  b <- count_reg(pmin(y, 3) ~ x + factor(u), d, family = "poisson")
  expect_equal(a$coefficients[["x"]], b$coefficients[["x"]])
  expect_identical(max(a$Y), 3)
  zc <- zi_count(pmin(y, 3) ~ x, d, family = "poisson", zero = ~ z)
  expect_identical(max(zc$y), 3L)
  cd <- count_reg(y ~ ., d[, c("y", "x", "u")], family = "poisson", fe = "u")
  expect_false("u" %in% names(cd$coefficients))
  expect_true("x" %in% names(cd$coefficients))
  ## zi_count: the inflation equation's default right-hand side is the count formula's covariates, so a '.'
  ## leaves the unit identifier out of it as well
  d$y0 <- ifelse(z > 0.5, 0L, d$y)
  zd <- zi_count(y0 ~ ., d[, c("y0", "x", "u")], family = "poisson", fe = "u", se = "none")
  zx <- zi_count(y0 ~ x, d, family = "poisson", fe = "u", se = "none")
  expect_identical(names(zd$zero.coefficients), c("(Intercept)", "x"))
  expect_equal(zd$coefficients, zx$coefficients)
  expect_equal(zd$loglik, zx$loglik)
})

test_that("a character covariate predicts on new data holding one of its values", {
  set.seed(10)
  dch <- data.frame(y = rpois(200, 3), x = rnorm(200), g = sample(c("a", "b", "c"), 200, TRUE), stringsAsFactors = FALSE)
  cr <- count_reg(y ~ x + g, dch, family = "poisson")
  gl <- glm(y ~ x + g, poisson, dch)
  expect_equal(unname(predict(cr, data.frame(x = 0, g = "b"))),
               unname(predict(gl, data.frame(x = 0, g = "b"), type = "response")), tolerance = 1e-6)
  expect_true(is.finite(predict(cpb(y ~ x + g, dch, truncated = FALSE, se = "none"), data.frame(x = 0, g = "b"))))
})

test_that("the zero-inflated EM fits report the observed-data log-likelihood", {
  skip_on_cran()
  set.seed(13); pz <- do.call(rbind, lapply(1:20, function(i) data.frame(u = i, x = rnorm(10))))
  pz$y <- ifelse(runif(200) < .3, 0L, rcpb(200, exp(2.7 + .1 * pz$x + rep(rnorm(20, 0, .2), each = 10)), .06))
  ze <- suppressWarnings(zi_cpb(y ~ x, pz, zero = ~ 1, fe = "u", method = "em"))
  f0 <- dcpb(0, ze$lambda_full, ze$alpha); py <- dcpb(ze$y, ze$lambda_full, ze$alpha)
  ll <- sum(ifelse(ze$y == 0, log(ze$pi_full + (1 - ze$pi_full) * f0), log1p(-ze$pi_full) + log(py)))
  expect_equal(ze$loglik, ll, tolerance = 1e-6)
  zp <- suppressWarnings(zi_cpb(y ~ x, pz, zero = ~ 1, method = "em"))
  f0 <- dcpb(0, zp$lambda_full, zp$alpha); py <- dcpb(zp$y, zp$lambda_full, zp$alpha)
  ll <- sum(ifelse(zp$y == 0, log(zp$pi_full + (1 - zp$pi_full) * f0), log1p(-zp$pi_full) + log(py)))
  expect_equal(zp$loglik, ll, tolerance = 1e-6)
  expect_identical(zp$max.support, 500L)
})

test_that("predict.zi_cpb returns NA for a row with a missing covariate", {
  set.seed(1); n <- 300; x <- rnorm(n); z <- rnorm(n)
  d <- data.frame(y = rzicpb(n, lambda = exp(1.3 + 0.5 * x), alpha = 0.5, pi = plogis(-0.5 + 0.8 * z)), x, z)
  zz <- zi_cpb(y ~ x, d, zero = ~ z)
  p <- predict(zz, newdata = data.frame(x = c(0, NA), z = c(0, 0)))
  expect_true(is.finite(p[1]) && is.na(p[2]))
})

test_that("zi_test() refuses nests on other data, truncated nests and other supports", {
  skip_on_cran()
  set.seed(14); x <- rnorm(300); d <- data.frame(y = rcpb(300, exp(.4 + .3 * x), .3), x = x)
  zk <- suppressWarnings(zi_cpb(y ~ x, d, zero = ~ 1))
  expect_error(zi_test(zk, cpb(y ~ x, d, truncated = FALSE, support = "floor")), "different CPB supports")
  expect_error(zi_test(zk, cpb(y ~ x, d[d$y > 0, ])), "different numbers of observations")
  expect_s3_class(zi_test(zk, cpb(y ~ x, d, truncated = FALSE)), "zi_test")
})

test_that("weighted fits count the weight total in glance() and compare_models()' BIC", {
  set.seed(7); d <- data.frame(y = rpois(300, 3), x = rnorm(300)); w <- sample(1:3, 300, TRUE)
  cw <- count_reg(y ~ x, d, family = "poisson", weights = w)
  expect_equal(underdisp:::glance.count_reg(cw)$nobs, sum(w))
  cm <- compare_models(a = cw, b = count_reg(y ~ x, d, family = "negbin", weights = w))
  expect_equal(cm["a", "BIC"], -2 * cw$loglik + cw$df * log(sum(w)))
})

test_that("max.support stops one short of the integer maximum", {
  d <- data.frame(y = rpois(50, 3), x = rnorm(50))
  expect_error(cpb(y ~ x, d, truncated = FALSE, max.support = .Machine$integer.max), "2147483646")
})
