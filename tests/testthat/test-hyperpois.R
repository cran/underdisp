## The hyper-Poisson family: pmf identities, the mean parameterization, oracle parity
## with DGLMExtPois, recovery, and inheritance of the two-part machinery.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

test_that("the hyper-Poisson pmf sums to one, has the Poisson at gamma = 1, and has mean mu", {
  k <- 0:120
  for (g in c(0.2, 0.6, 1.5, 4)) for (mu in c(0.4, 3.2, 15)) {
    p <- dhyperpois(k, mu, g)
    expect_equal(sum(p), 1, tolerance = 1e-10)
    expect_equal(sum(k * p), mu, tolerance = 1e-8)
    v <- sum(k^2 * p) - mu^2
    if (g < 1) expect_lt(v, mu) else expect_gt(v, mu)
    expect_equal(underdisp:::.count_families$hyperpois$varfun(mu, g), v, tolerance = 1e-8)
  }
  expect_equal(dhyperpois(k, 3.2, 1), dpois(k, 3.2), tolerance = 1e-12)
})

test_that("the compiled normalizer and moments match the R terms", {
  lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  for (ll in c(-2, 0.5, 2.7)) for (g in c(0.3, 1, 3)) {
    lt <- underdisp:::.hp_lterm(0:2000, ll, g); p <- exp(lt - lse(lt))
    r <- underdisp:::hp_norm_cpp(ll, g)
    expect_equal(r[1, 1], lse(lt), tolerance = 1e-10)
    expect_equal(r[1, 2], sum((0:2000) * p), tolerance = 1e-8)
    expect_equal(r[1, 3], sum((0:2000)^2 * p) - sum((0:2000) * p)^2, tolerance = 1e-8)
  }
})

test_that("the pmf matches DGLMExtPois at the solved rate", {
  skip_if_not_installed("DGLMExtPois")
  k <- 0:30
  for (g in c(0.3, 0.8, 2)) for (mu in c(0.7, 4)) {
    lam <- exp(underdisp:::hp_loglambda_cpp(mu, g))
    expect_equal(dhyperpois(k, mu, g), DGLMExtPois::dhP(k, gamma = g, lambda = lam), tolerance = 1e-9)
  }
})

test_that("d/p/q/r are mutually consistent and handle zero-length input", {
  p <- phyperpois(0:12, 3.2, 0.5)
  expect_equal(p, cumsum(dhyperpois(0:12, 3.2, 0.5)), tolerance = 1e-10)
  expect_equal(qhyperpois(p[1:8] - 1e-9, 3.2, 0.5), 0:7)
  expect_identical(dhyperpois(numeric(0), 3.2, 0.5), numeric(0))
  expect_identical(qhyperpois(numeric(0), 3.2, 0.5), numeric(0))
  expect_identical(rhyperpois(0, 3.2, 0.5), numeric(0))
  expect_error(dhyperpois(0:3, 3.2, c(0.5, 0.5)), "single value")
  set.seed(1); s <- rhyperpois(4000, 3.2, 0.5)
  expect_equal(mean(s), 3.2, tolerance = 0.03)
  expect_lt(var(s), 3.2)
})

test_that("count_reg recovers the dispersion and the coefficients and reproduces DGLMExtPois", {
  set.seed(12); n <- 500; x <- rnorm(n); mu <- exp(1 + 0.5 * x)
  d <- data.frame(y = rhyperpois(n, mu, gamma = 0.4), x = x)
  m <- count_reg(y ~ x, d, family = "hyperpois")
  expect_equal(m$theta, 0.4, tolerance = 0.35); expect_equal(unname(coef(m)[2]), 0.5, tolerance = 0.1)
  expect_equal(mean(fitted(m)), mean(d$y), tolerance = 0.03)
  expect_equal(count_reg(y ~ x, d, family = "hyper_poisson")$loglik, m$loglik)   # alias
  dt <- dispersion_test(m)
  expect_true(dt$p.value < 0.05)
  skip_if_not_installed("DGLMExtPois")
  ref <- suppressWarnings(DGLMExtPois::glm.hP(y ~ x, y ~ 1, data = d))   # logver: minus the log-likelihood
  expect_gte(m$loglik, -ref$logver - 0.01)
  expect_equal(m$theta, exp(unname(ref$deltas)), tolerance = 0.05)
})

test_that("the hyper-Poisson inherits the hurdle, zero-inflated, scoring, and simulation machinery", {
  set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1, rhyperpois(n, exp(1 + 0.4 * x), 0.5) + 1L, 0L)
  d <- data.frame(y = y, x = x, z = z)
  h <- hurdle_count(y ~ x, d, family = "hyperpois", participation = ~ z)
  expect_true(is.finite(h$loglik))
  expect_equal(predict(h, newdata = d), fitted(h), tolerance = 1e-8)
  zz <- zi_count(y ~ x, d, family = "hyperpois", zero = ~ z)
  expect_true(is.finite(zz$loglik))
  m <- count_reg(y ~ x, d, family = "hyperpois")
  s <- simulate(m, nsim = 5, seed = 1)
  expect_equal(dim(s), c(n, 5L)); expect_true(all(s >= 0))
  expect_equal(predict(m, newdata = d), fitted(m), tolerance = 1e-8)
  expect_s3_class(first_difference(m, "x", 0, 1), "ud_fd")
  expect_true(all(is.finite(confint(m))))
  expect_true(all(is.finite(score(m))))
  mt <- count_reg(y ~ x, d[d$y > 0, ], family = "hyperpois", truncated = TRUE)
  expect_true(is.finite(mt$loglik))
})

test_that("the normalizer is exact at large gamma, and means the support cap cannot hold are refused", {
  hn <- underdisp:::hp_norm_cpp; hll <- underdisp:::hp_loglambda_cpp
  for (g in c(1e4, 1e8, 1e15)) expect_equal(hn(log(0.6 * g), g)[1, 2], 0.6 / 0.4, tolerance = 1e-3)   # the geometric limit
  for (g in c(1e-6, 0.3, 1, 7, 1e4, 1e8)) {
    mu <- c(0.2, 2, 30); expect_equal(hn(hll(mu, g), g)[, 2], mu, tolerance = 1e-8)
  }
  expect_true(is.na(hll(1e9, 0.5)))
  expect_lt(system.time(hll(rep(1e9, 50), 0.5))[["elapsed"]], 5)
})
