## The four families added in 0.1.1: pmf identities, oracle parity, recovery,
## and inheritance of the two-part machinery.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

test_that("each family's pmf sums to one and has the Poisson as its equidispersed member", {
  k <- 0:80
  expect_equal(sum(dgammacount(k, 3.2, 2.5)), 1, tolerance = 1e-10)
  expect_equal(sum(ddoublepois(k, 3.2, 2)), 1, tolerance = 1e-10)
  expect_equal(sum(dgenpois(k, 3.2, -0.3)), 1, tolerance = 1e-10)
  expect_equal(sum(dgenpois(k, 3.2, 0.3)), 1, tolerance = 1e-8)
  expect_equal(sum(dcompois(k, nu = 1.8, mu = 3.2)), 1, tolerance = 1e-10)
  expect_equal(dgammacount(k, 3.2, 1), dpois(k, 3.2), tolerance = 1e-12)
  expect_equal(ddoublepois(k, 3.2, 1), dpois(k, 3.2), tolerance = 1e-12)
  expect_equal(dgenpois(k, 3.2, 0), dpois(k, 3.2), tolerance = 1e-12)
  expect_equal(dcompois(k, nu = 1, mu = 3.2), dpois(k, 3.2), tolerance = 1e-10)
})

test_that("the pmfs match independent implementations", {
  k <- 0:40
  skip_if_not_installed("rmutil")
  expect_equal(dgammacount(k, 3.2, 2.5), rmutil::dgammacount(k, m = 3.2, s = 2.5), tolerance = 1e-12)
  skip_if_not_installed("gamlss.dist")
  for (th in c(0.5, 2, 4))
    expect_equal(ddoublepois(k, 3.2, th), gamlss.dist::dDPO(k, mu = 3.2, sigma = 1 / th), tolerance = 1e-12)
  expect_equal(dgenpois(k, 3.2, 0.3), VGAM::dgenpois0(k, theta = 3.2 * 0.7, lambda = 0.3), tolerance = 1e-12)
})

test_that("the mean parameterization of the COM-Poisson recovers its mean and variance", {
  ll <- underdisp:::cmp_loglambda_cpp(c(0.5, 3.2, 20), 1.8)
  m <- underdisp:::cmp_moments_cpp(ll, 1.8)
  expect_equal(m[, 1], c(0.5, 3.2, 20), tolerance = 1e-8)
  expect_true(all(m[, 2] < m[, 1]))                      # underdispersed at nu = 1.8
  expect_equal(dcompois(0:60, nu = 1.8, mu = 3.2), dcompois(0:60, lambda = exp(ll[2]), nu = 1.8), tolerance = 1e-12)
})

test_that("d/p/q are mutually consistent and handle zero-length input", {
  for (fam in list(list(d = dgammacount, p = pgammacount, q = qgammacount, r = rgammacount, th = 2),
                   list(d = ddoublepois, p = pdoublepois, q = qdoublepois, r = rdoublepois, th = 2),
                   list(d = dgenpois, p = pgenpois, q = qgenpois, r = rgenpois, th = -0.3))) {
    p <- fam$p(0:12, 3.2, fam$th)
    expect_equal(p, cumsum(fam$d(0:12, 3.2, fam$th)), tolerance = 1e-10)
    expect_equal(fam$q(p[1:8] - 1e-9, 3.2, fam$th), 0:7)
    expect_identical(fam$d(numeric(0), 3.2, fam$th), numeric(0))
    expect_identical(fam$q(numeric(0), 3.2, fam$th), numeric(0))
    expect_identical(fam$r(0, 3.2, fam$th), numeric(0))
    expect_error(fam$d(0:3, 3.2, c(fam$th, fam$th)), "single value")
    set.seed(1); s <- fam$r(4000, 3.2, fam$th)
    m_pmf <- sum((0:60) * fam$d(0:60, 3.2, fam$th))                 # the pmf's own mean (mu is the rate)
    expect_equal(mean(s), m_pmf, tolerance = 0.06)
  }
  expect_identical(dcpb(numeric(0), 3, 0.5), numeric(0))
  expect_identical(dgec(numeric(0), 3, 0.7), numeric(0))
  expect_identical(qgec(numeric(0), 3, 0.7), numeric(0))
  expect_identical(dcompois(numeric(0), 3, 1.5), numeric(0))
})

test_that("count_reg recovers each family's dispersion and reports the exact mean", {
  set.seed(11); n <- 500; x <- rnorm(n); mu <- exp(1 + 0.5 * x)
  d <- data.frame(y = rgammacount(n, mu, alpha = 2), x = x)
  m <- count_reg(y ~ x, d, family = "gammacount")
  expect_equal(m$theta, 2, tolerance = 0.25); expect_equal(unname(coef(m)[2]), 0.5, tolerance = 0.15)
  expect_equal(mean(fitted(m)), mean(d$y), tolerance = 0.03)
  d <- data.frame(y = rdoublepois(n, mu, theta = 2), x = x)
  m <- count_reg(y ~ x, d, family = "doublepois")
  expect_equal(m$theta, 2, tolerance = 0.3)
  d <- data.frame(y = rgenpois(n, mu, lambda = -0.3), x = x)
  m <- count_reg(y ~ x, d, family = "genpois")
  expect_equal(m$theta, -0.3, tolerance = 0.4)
  d <- data.frame(y = rcompois(n, nu = 1.8, mu = mu), x = x)
  m <- count_reg(y ~ x, d, family = "mpcmp")
  expect_equal(m$theta, 1.8, tolerance = 0.25); expect_equal(unname(coef(m)[2]), 0.5, tolerance = 0.15)
  expect_equal(mean(fitted(m)), mean(d$y), tolerance = 0.03)
  expect_equal(count_reg(y ~ x, d, family = "compois_mean")$loglik, m$loglik)   # alias
})

test_that("the generalized-Poisson and mean-parameterized COM-Poisson fits reproduce glmmTMB", {
  skip_on_cran(); skip_if_not_installed("glmmTMB")
  set.seed(3); n <- 400; x <- rnorm(n); mu <- exp(1 + 0.4 * x)
  d <- data.frame(y = rcompois(n, nu = 1.8, mu = mu), x = x)
  ours <- count_reg(y ~ x, d, family = "mpcmp")
  ref <- suppressWarnings(glmmTMB::glmmTMB(y ~ x, data = d, family = glmmTMB::compois()))
  expect_gte(ours$loglik, as.numeric(logLik(ref)) - 0.01)
  expect_equal(ours$theta, 1 / sigma(ref), tolerance = 0.05)
  d2 <- data.frame(y = rgenpois(n, mu, -0.3), x = x)
  ours2 <- count_reg(y ~ x, d2, family = "genpois")
  ref2 <- suppressWarnings(glmmTMB::glmmTMB(y ~ x, data = d2, family = glmmTMB::genpois()))
  expect_gte(ours2$loglik, as.numeric(logLik(ref2)) - 0.01)
  expect_equal(ours2$theta, 1 - 1 / sqrt(sigma(ref2)), tolerance = 0.05)
})

test_that("the new families inherit the hurdle, zero-inflated, scoring, and simulation machinery", {
  set.seed(5); n <- 400; x <- rnorm(n); z <- rnorm(n)
  y <- ifelse(rbinom(n, 1, plogis(0.3 + 0.8 * z)) == 1, rgammacount(n, exp(1 + 0.4 * x), 2) + 1L, 0L)
  d <- data.frame(y = y, x = x, z = z)
  h <- hurdle_count(y ~ x, d, family = "gammacount", participation = ~ z)
  expect_true(is.finite(h$loglik))
  expect_equal(predict(h, newdata = d), fitted(h), tolerance = 1e-8)
  expect_true(all(is.finite(score(h))))
  zz <- zi_count(y ~ x, d, family = "genpois", zero = ~ z)
  expect_true(is.finite(zz$loglik)); expect_true(all(is.finite(confint(zz))))
  for (fam in c("gammacount", "doublepois", "genpois", "mpcmp")) {
    m <- count_reg(y ~ x, d, family = fam)
    s <- simulate(m, nsim = 5, seed = 1)
    expect_equal(dim(s), c(n, 5L)); expect_true(all(s >= 0))
    expect_equal(predict(m, newdata = d), fitted(m), tolerance = 1e-8)
    expect_s3_class(first_difference(m, "x", 0, 1), "ud_fd")
    expect_true(all(is.finite(confint(m))))
    expect_true(all(is.finite(score(m))))
  }
  cm <- compare_models(gc = count_reg(y ~ x, d, family = "gammacount"), po = count_reg(y ~ x, d, family = "poisson"))
  expect_equal(nrow(cm), 2L)
})

test_that("the compiled double-Poisson and generalized-Poisson normalizers match the R terms", {
  lse <- function(x) { m <- max(x); m + log(sum(exp(x - m))) }
  for (mu in c(0.3, 2.5, 14)) for (th in c(0.5, 2, 6)) {
    K <- underdisp:::.dp_kmax(mu, th); lt <- underdisp:::.dp_lterm(0:K, mu, th); p <- exp(lt - lse(lt))
    r <- underdisp:::dp_norm_cpp(mu, th, 100000L)
    expect_equal(r[1, 1], lse(lt), tolerance = 1e-10); expect_equal(r[1, 2], sum((0:K) * p), tolerance = 1e-8)
  }
  for (mu in c(0.3, 2.5, 14)) for (la in c(-0.4, -0.05, 0.3)) {
    K <- underdisp:::.gp_kmax(mu, la); lt <- underdisp:::.gp_lterm(0:K, mu, la); ok <- is.finite(lt)
    r <- underdisp:::gp_norm_cpp(mu, la, 100000L)
    expect_true(is.finite(r[1, 1]) && abs(r[1, 1] - lse(lt[ok])) < 1e-9,
                info = sprintf("mu %g lambda %g: C++ logZ %.12g, R logZ %.12g, K %d, finite terms %d", mu, la, r[1, 1], lse(lt[ok]), K, sum(ok)))
    if (la < 0) { p <- exp(lt[ok] - lse(lt[ok])); expect_equal(r[1, 2], sum((0:K)[ok] * p), tolerance = 1e-8) }
  }
})
