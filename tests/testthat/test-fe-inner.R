## The concentrated fixed-effects likelihood: each unit's intercept is solved to
## the supremum of a saw-tooth objective (src/cpb_fe.cpp, src/gec_fe.cpp). These
## tests pin the inner search against an R-side brute force and check that the
## reported log-likelihood is attained at the reported intercepts.

testthat::skip_on_cran()   # validation battery: runs in the package's CI (NOT_CRAN = true), not on CRAN

brute_unit <- function(yy, oo, alpha, width = 4, by = 0.002) {
  if (sum(yy) == 0) return(0)                 # all-zero unit: the supremum is P(0) -> 1 as the rate -> 0
  ll <- function(a) { v <- sum(dcpb(yy, exp(a + oo), alpha, log = TRUE)); if (is.finite(v)) v else -1e300 }
  ap <- log(sum(yy) / sum(exp(oo)))
  g <- seq(ap - width, ap + width, by = by); vals <- vapply(g, ll, numeric(1))
  k <- which.max(vals)
  o <- stats::optimize(ll, c(g[max(1, k - 1)], g[min(length(g), k + 1)]), maximum = TRUE, tol = 1e-10)
  max(o$objective, vals[k])
}

test_that("cpb_fe's inner search is never below a brute-force profile and its loglik is attained", {
  set.seed(11); n <- 300; x <- rnorm(n)
  d <- data.frame(y = rcpb(n, exp(1.0 + 0.5 * x), 0.5), x = x, u = factor(sample(1:12, n, TRUE)))
  f <- cpb_fe(y ~ x, d, fe = "u")
  d <- d[order(d$u), ]
  ## (a) the reported loglik is the CPB log-likelihood at the reported parameters
  lam <- exp(f$fe[as.character(d$u)] + f$coefficients[["x"]] * d$x)
  expect_equal(sum(dcpb(d$y, lam, f$alpha, log = TRUE)), f$loglik, tolerance = 1e-8)
  ## (b) per unit, the compiled inner search at the fitted (beta, alpha) is at
  ##     least the brute-force supremum (a 0.002 grid refined by optimize())
  tot <- 0
  for (u in levels(d$u)) {
    idx <- d$u == u
    tot <- tot + brute_unit(d$y[idx], f$coefficients[["x"]] * d$x[idx], f$alpha)
  }
  expect_gte(f$loglik, tot - 1e-4)   # the supremum is a left limit: attained to O(1e-8) in the intercept
  ## (c) the fitted intercepts are feasible: every count is below its ceiling
  expect_true(all(d$y <= floor(lam / (1 - f$alpha) + 1e-9)))
})

test_that("cpb_fe's inner search handles small ceilings, where the teeth are wide and tall", {
  set.seed(5); n <- 400; x <- rnorm(n); u <- factor(sample(1:16, n, TRUE)); ue <- rnorm(16, 0, 0.6)[u]
  d <- data.frame(y = rcpb(n, exp(0.2 + 0.3 * x + ue), 0.4), x = x, u = u)
  f <- cpb_fe(y ~ x, d, fe = "u")
  d <- d[order(d$u), ]
  lam <- exp(f$fe[as.character(d$u)] + f$coefficients[["x"]] * d$x)
  expect_equal(sum(dcpb(d$y, lam, f$alpha, log = TRUE)), f$loglik, tolerance = 1e-8)
  tot <- 0
  for (u in levels(d$u)) {
    idx <- d$u == u
    tot <- tot + brute_unit(d$y[idx], f$coefficients[["x"]] * d$x[idx], f$alpha, width = 5, by = 0.001)
  }
  expect_gte(f$loglik, tot - 1e-4)   # the supremum is a left limit: attained to O(1e-8) in the intercept
})

test_that("gec_fe's inner search matches a brute-force profile for an underdispersed Katz fit", {
  set.seed(4); n <- 300; x <- rnorm(n)
  d <- data.frame(y = rcpb(n, exp(0.8 + 0.4 * x), 0.5), x = x, u = factor(sample(1:10, n, TRUE)))
  f <- gec_fe(y ~ x, d, fe = "u")
  expect_lt(f$delta, 1)
  d <- d[order(d$u), ]
  mu <- exp(f$fe[as.character(d$u)] + f$coefficients[["x"]] * d$x)
  expect_equal(sum(dgec(d$y, mu, f$delta, log = TRUE)), f$loglik, tolerance = 1e-8)
  tot <- 0
  for (u in levels(d$u)) {
    idx <- d$u == u; yy <- d$y[idx]; oo <- f$coefficients[["x"]] * d$x[idx]
    ll <- function(a) { v <- sum(dgec(yy, exp(a + oo), f$delta, log = TRUE)); if (is.finite(v)) v else -1e300 }
    ap <- log(sum(yy) / sum(exp(oo))); g <- seq(ap - 4, ap + 4, by = 0.002); vals <- vapply(g, ll, numeric(1))
    k <- which.max(vals)
    o <- stats::optimize(ll, c(g[max(1, k - 1)], g[min(length(g), k + 1)]), maximum = TRUE, tol = 1e-10)
    tot <- tot + max(o$objective, vals[k])
  }
  expect_gte(f$loglik, tot - 1e-4)   # the supremum is a left limit: attained to O(1e-8) in the intercept
})

test_that("the concentrated CPB profile in alpha is smooth at the fit (no saw-tooth in the outer objective)", {
  set.seed(11); n <- 300; x <- rnorm(n)
  d <- data.frame(y = rcpb(n, exp(1.0 + 0.5 * x), 0.5), x = x, u = factor(sample(1:12, n, TRUE)))
  f <- cpb_fe(y ~ x, d, fe = "u")
  d <- d[order(d$u), ]
  X <- matrix(d$x, ncol = 1); Y <- as.integer(d$y); uf <- factor(d$u); nu <- nlevels(uf)
  ustart <- as.integer(c(0, cumsum(tabulate(as.integer(uf), nu))))
  off <- rep(0, n); wv <- rep(1, n); ms <- f$max.support
  prof <- vapply(f$alpha + seq(-0.02, 0.02, by = 0.004), function(a)
    -underdisp:::cpb_fe_nll_cpp(c(f$coefficients[["x"]], qlogis(a)), X, Y, off, wv, ustart, nu, ms, FALSE, 30L, numeric(0))[[1L]],
    numeric(1))
  ## the profile peaks at the fit and its second differences are all negative (concave, no jags)
  expect_equal(which.max(prof), 6L)
  expect_true(all(diff(prof, differences = 2) < 0))
})
