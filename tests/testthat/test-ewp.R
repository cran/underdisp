## The three-parameter exponentially weighted Poisson: pmf identities, oracle parity
## with the ewp package, recovery, the two-parameter dispersion test, and inheritance
## of the two-part machinery.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

test_that("the EWP pmf sums to one, is the Poisson at beta1 = beta2 = 0, and is underdispersed for positive betas", {
  k <- 0:150
  for (b in list(c(0.5, 0.5), c(1.5, 0.3), c(-0.3, 0.8), c(2, 2))) for (lam in c(0.6, 3.2, 14)) {
    p <- dewp(k, lam, b[1], b[2])
    expect_equal(sum(p), 1, tolerance = 1e-10)
    m <- sum(k * p); v <- sum(k^2 * p) - m^2
    expect_equal(underdisp:::.count_families$ewp3$meanfun(lam, b), m, tolerance = 1e-8)
    expect_equal(underdisp:::.count_families$ewp3$varfun(lam, b), v, tolerance = 1e-8)
    if (all(b > 0)) expect_lt(v, m)
  }
  expect_equal(dewp(k, 3.2, 0, 0), dpois(k, 3.2), tolerance = 1e-12)
})

test_that("the compiled normalizer matches the R terms, including a thickened upper tail", {
  lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  for (ll in c(-1.5, 1.2, 2.8)) for (b in list(c(0.5, 0.5), c(1, -0.5), c(-0.4, 2))) {
    lt <- underdisp:::.ewp_lterm(0:3000, ll, b); p <- exp(lt - lse(lt))
    r <- underdisp:::ewp_norm_cpp(ll, b[1], b[2])
    expect_equal(r[1, 1], lse(lt), tolerance = 1e-10)
    expect_equal(r[1, 2], sum((0:3000) * p), tolerance = 1e-8)
  }
})

test_that("the pmf matches the ewp package", {
  skip_if_not_installed("ewp")
  k <- 0:20
  for (b in list(c(0.5, 0.5), c(1.5, 0.3), c(-0.3, 0.8))) for (lam in c(0.6, 3.2, 8))
    expect_equal(dewp(k, lam, b[1], b[2]), ewp::dewp3(k, lam, b[1], b[2], sum_limit = 150), tolerance = 1e-9)
})

test_that("d/p/q/r are mutually consistent, take two dispersion values, and handle zero-length input", {
  p <- pewp(0:12, 3.2, 1, 1)
  expect_equal(p, cumsum(dewp(0:12, 3.2, 1, 1)), tolerance = 1e-10)
  expect_equal(qewp(p[1:8] - 1e-9, 3.2, 1, 1), 0:7)
  expect_identical(dewp(numeric(0), 3.2, 1, 1), numeric(0))
  expect_identical(rewp(0, 3.2, 1, 1), numeric(0))
  expect_error(underdisp:::.fam_d(underdisp:::.count_families[["ewp3"]], 0:3, 3.2, 1, FALSE), "2 single values")
  set.seed(1); s <- rewp(4000, 3.2, 1, 1)
  m <- sum((0:60) * dewp(0:60, 3.2, 1, 1))
  expect_equal(mean(s), m, tolerance = 0.03)
  expect_true(any(s == 0) || dewp(0, 3.2, 1, 1) < 1e-3)        # zeros are in the support
})

test_that("count_reg recovers the EWP's parameters, tests both against the Poisson, and reproduces ewp_reg", {
  set.seed(21); n <- 600; x <- rnorm(n); lam <- exp(1 + 0.4 * x)
  d <- data.frame(y = vapply(lam, function(l) rewp(1, l, 1, 0.6), numeric(1)), x = x)
  m <- count_reg(y ~ x, d, family = "ewp", se = "analytic")
  expect_length(m$theta, 2L)
  expect_equal(m$theta, c(1, 0.6), tolerance = 0.5)
  expect_equal(unname(coef(m)[2]), 0.4, tolerance = 0.15)
  expect_equal(mean(fitted(m)), mean(d$y), tolerance = 0.03)
  ci <- confint(m); expect_true(all(c("beta1", "beta2") %in% rownames(ci)))
  dt <- dispersion_test(m, B = 49)
  expect_equal(unname(dt$parameter), 2)
  expect_true(dt$p.value < 0.05)
  expect_error(dispersion_test(m, alternative = "under"), "two-sided")
  skip_if_not_installed("ewp")
  ref <- suppressWarnings(ewp::ewp_reg(y ~ x, data = d, verbose = FALSE, sum_limit = 150))
  expect_gte(m$loglik, as.numeric(logLik(ref)) - 0.01)
})

test_that("the EWP inherits the hurdle, zero-inflated, scoring, and simulation machinery", {
  set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1,
              vapply(exp(1 + 0.4 * x), function(l) rewp(1, l, 1, 1), numeric(1)) + 1L, 0L)
  d <- data.frame(y = y, x = x, z = z)
  h <- hurdle_count(y ~ x, d, family = "ewp3", participation = ~ z)
  expect_true(is.finite(h$loglik))
  expect_equal(predict(h, newdata = d), fitted(h), tolerance = 1e-8)
  zz <- zi_count(y ~ x, d, family = "ewp3", zero = ~ z)
  expect_true(is.finite(zz$loglik)); expect_length(zz$theta, 2L)
  m <- count_reg(y ~ x, d, family = "ewp3")
  s <- simulate(m, nsim = 5, seed = 1)
  expect_equal(dim(s), c(n, 5L)); expect_true(all(s >= 0))
  expect_equal(predict(m, newdata = d), fitted(m), tolerance = 1e-8)
  expect_s3_class(first_difference(m, "x", 0, 1), "ud_fd")
  expect_true(all(is.finite(score(m))))
  expect_output(print(m), "beta1, beta2")
})

test_that("the two-parameter EWP is the three-parameter one with equal tails", {
  k <- 0:60
  expect_equal(dewp(k, 3.2, 0.8), dewp(k, 3.2, 0.8, 0.8))
  f2 <- underdisp:::.count_families$ewp2; f3 <- underdisp:::.count_families$ewp3
  expect_equal(f2$logpmf(k, rep(3.2, length(k)), 0.8), f3$logpmf(k, rep(3.2, length(k)), c(0.8, 0.8)))
  set.seed(22); n <- 500; x <- rnorm(n); lam <- exp(1 + 0.4 * x)
  d <- data.frame(y = vapply(lam, function(l) rewp(1, l, 1), numeric(1)), x = x)
  m2 <- count_reg(y ~ x, d, family = "ewp2"); m3 <- count_reg(y ~ x, d, family = "ewp3")
  expect_equal(m2$theta, 1, tolerance = 0.3)
  expect_gte(m3$loglik, m2$loglik - 1e-6)                     # nested: EWP2 is EWP3 at beta1 = beta2
  dt <- dispersion_test(m2)
  expect_equal(unname(dt$parameter), 1); expect_true(dt$p.value < 0.05)
})

test_that("dispersion_test() of the EWP3's two tails bootstraps by default; B = 0 gives the chi-square(2)", {
  set.seed(1); d <- data.frame(y = rpois(200, 3))
  f <- suppressWarnings(count_reg(y ~ 1, d, family = "ewp3", se = "none"))
  dt <- dispersion_test(f, B = 0)
  expect_false(dt$asymptotic_ok)
  expect_identical(dt$calibration, "asymptotic")
  expect_equal(dt$p.value, stats::pchisq(unname(dt$statistic), 2, lower.tail = FALSE))
})
